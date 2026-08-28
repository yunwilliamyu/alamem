use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::{Arc, Mutex};
use needletail::{parse_fastx_file};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use clap::{Parser};
use memchr;
use libc;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
use alamem::{YIndex,
    process_query_sequence,
    filter_overlapping_hits,
    encode_sequence,
    reverse_complement,
    AlignConfig, CONFIG,
    ThreadBuffers
};


fn resolve_input(path: &str) -> Vec<String> {
    if path.ends_with(".txt") {
        let file: File = File::open(path).expect("Failed to open input list");
        let reader: BufReader<File> = BufReader::new(file);
        reader.lines().filter_map(|l: Result<String, std::io::Error>| l.ok()).collect()
    } else {
        vec![path.to_string()]
    }
}

pub struct FlatChunk {
    pub memory: Vec<u8>,
    pub bounds: Vec<(usize, usize, usize, usize)>,
}

impl FlatChunk {
    pub fn new(byte_capacity: usize, item_capacity: usize) -> Self {
        Self {
            memory: Vec::with_capacity(byte_capacity),
            bounds: Vec::with_capacity(item_capacity),
        }
    }
}

#[derive(Parser, Debug)]
#[command(name = "medival_search")]
struct MedivalSearchCli {
    /// The streamed database. Can be .fasta, .fa, .gz, or .txt (listing paths to fastas)
    pub database: String,
    /// The indexed query file(s). Can be .fasta, .fa, .gz, or .txt (listing paths to fastas)
    pub y_files: String,
    /// TSV file with all the hits
    pub output_path: String,

    /// Thread limit (set to 0 for no limit)
    #[arg(short = 't', long, default_value_t = 0)]
    pub threads: usize,

    #[command(flatten)]
    pub align_config: AlignConfig,
}

fn compute_chunk(
    mut chunk: FlatChunk,
    text_buffer: &mut Vec<u8>,
    thread_bufs: &mut ThreadBuffers,
    y_index: &YIndex,
    writer: &Arc<Mutex<BufWriter<File>>>,
    pb_main: &ProgressBar,
    pb_seqs: &ProgressBar,
    is_single_file: bool,
    recycle_tx: &crossbeam::channel::Sender<FlatChunk>,
) {
    let chunk_len = chunk.bounds.len();
    let chunk_bytes = chunk.memory.len();

    for &(id_start, id_end, seq_start, seq_end) in &chunk.bounds {
        let id = &chunk.memory[id_start..id_end];
        let raw_seq = &chunk.memory[seq_start..seq_end];

        thread_bufs.encoded_fwd.clear();
        thread_bufs.encoded_rev.clear();
        encode_sequence(&raw_seq, &mut thread_bufs.encoded_fwd);
        reverse_complement(&thread_bufs.encoded_fwd, &mut thread_bufs.encoded_rev);

        process_query_sequence(y_index, thread_bufs);
        filter_overlapping_hits(&mut thread_bufs.hits);

        if thread_bufs.hits.is_empty() { continue; }

        let q_id_str = unsafe { std::str::from_utf8_unchecked(id) };
        for hit in &thread_bufs.hits {
            let t_id_str = unsafe { std::str::from_utf8_unchecked(&y_index.seq_names[hit.t_id]) };
            let _ = writeln!(
                text_buffer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}",
                q_id_str, t_id_str, hit.q_size, hit.q_start, hit.q_end,
                hit.t_size, hit.t_start, hit.t_end, hit.strand, hit.ani, hit.score
            );
        }

        if text_buffer.len() > 32 * 1024 {
            writer.lock().unwrap().write_all(text_buffer).unwrap();
            text_buffer.clear();
        }
    }

    if is_single_file { pb_main.inc(chunk_bytes as u64); }
    pb_seqs.inc(chunk_len as u64);

    chunk.memory.clear();
    chunk.bounds.clear();
    recycle_tx.send(chunk).ok();
}


fn main() {
    let cli = MedivalSearchCli::parse();
    eprintln!("{:#?}", cli);

    let x_files = resolve_input(&cli.database);
    let total_x: usize = x_files.len();
    let y_files = resolve_input(&cli.y_files);

    let total_threads = if cli.threads == 0 {
        std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4) // Default to 4 threads if can't figure it out
    } else {
        cli.threads
    };

    rayon::ThreadPoolBuilder::new()
        .stack_size(8 * 1024 * 1024)
        .num_threads(total_threads)
        .build_global()
        .unwrap();

    let config = cli.align_config.clone();
    CONFIG.set(cli.align_config).expect("Config already initialized");

    let output_path = &cli.output_path;

    let out_file: File = File::create(output_path).unwrap_or_else(|_| panic!("Failed to create output file: {}", output_path));
    let writer: Arc<Mutex<BufWriter<File>>> = Arc::new(Mutex::new(BufWriter::new(out_file)));

    {
        let mut handle = writer.lock().unwrap();
        // Because I messed up naming of variables, so the column labels don't correspond to variable names
        // TODO: refactor so variable names match column labels
        // (refactor needed is in the variable names for the compute_chunk text_buffer)
        // For now, use temporary cludge of just renaming the columns.
        //writeln!(handle, "Query\tTarget\tQ_Size\tQ_Start\tQ_End\tT_Size\tT_Start\tT_End\tStrand\tANI\tScore").unwrap();
        writeln!(handle, "Reference\tQuery\tR_Size\tR_Start\tR_End\tQ_Size\tQ_Start\tQ_End\tStrand\tANI\tScore").unwrap();
    }

    let y_index: YIndex = YIndex::build(&y_files);
    eprintln!("Y Indexed. {} distinct {}-mers.", y_index.kmer_map.len(), config.kmer_size);


    eprintln!("Streaming database ({} files total)...", total_x);

    let m = MultiProgress::new();

    // 1. Branch progress bar setup dynamically based on file count
    let (pb_main, is_single_file) = if total_x == 1 {
        let single_file_bytes = std::fs::metadata(&x_files[0])
        .map(|m| m.len())
        .unwrap_or(0);

        let pb = m.add(ProgressBar::new(single_file_bytes));
        pb.set_style(
            ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bytes}/{total_bytes} ({percent}%) ETA: {eta} [{wide_bar:.cyan/blue}]\n{msg}")
            .unwrap()
            .progress_chars("=>-")
        );
        (pb, true)
    } else {
        let pb = m.add(ProgressBar::new(total_x as u64));
        pb.set_style(
            ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {pos}/{len} files ({percent}%) ETA: {eta} [{wide_bar:.cyan/blue}] \n{msg}")
            .unwrap()
            .progress_chars("=>-")
        );
        (pb, false)
    };
    pb_main.enable_steady_tick(std::time::Duration::from_millis(200));

    // 2. Sequence Spinner (Always shows running elapsed time)
    let pb_seqs = m.add(ProgressBar::new_spinner());
    pb_seqs.set_style(
        ProgressStyle::default_spinner()
        .template("{spinner:.green} {human_pos} seqs processed ({per_sec}) | {msg}")
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
    );
    pb_seqs.enable_steady_tick(std::time::Duration::from_millis(200));

    let (recycle_tx, recycle_rx) = crossbeam::channel::unbounded::<FlatChunk>();

    const MAX_SEQ: usize = 100_000;
    const MAX_BYTES: usize = 4 * 1024 * 1024;
    const CAPACITY_BYTES: usize = MAX_BYTES * 2;

    let global_remainder = Arc::new(Mutex::new((0usize, Vec::<FlatChunk>::new())));

    x_files.into_par_iter().for_each(|x_file| {
        let pb_main_clone = pb_main.clone();
        if !is_single_file { pb_main_clone.set_message(x_file.clone()); }

        let mut reader = parse_fastx_file(&x_file).unwrap_or_else(|_| panic!("Failed to open: {}", x_file));

        let get_chunk = || {
            recycle_rx.try_recv().unwrap_or_else(|_| {
                let mut c = FlatChunk::new(CAPACITY_BYTES, MAX_SEQ);
                c.memory.reserve_exact(CAPACITY_BYTES);
                c.bounds.reserve_exact(MAX_SEQ);
                c
            })
        };

        let mut chunk = get_chunk();
        let chunk_iterator = std::iter::from_fn(|| {
            while let Some(record) = reader.next() {
                let rec = record.expect("Invalid record");
                let id_slice = rec.id();

                let id_len = memchr::memchr2(b' ', b'\t', id_slice).unwrap_or(id_slice.len());
                let id_word = &id_slice[..id_len];

                let id_start = chunk.memory.len();
                chunk.memory.extend_from_slice(id_word);
                let id_end = chunk.memory.len();

                let seq_slice = rec.seq();
                let seq_start = chunk.memory.len();
                chunk.memory.extend_from_slice(&seq_slice);
                let seq_end = chunk.memory.len();

                chunk.bounds.push((id_start, id_end, seq_start, seq_end));

                if chunk.bounds.len() >= MAX_SEQ || chunk.memory.len() >= MAX_BYTES {
                    let full_chunk = std::mem::replace(&mut chunk, get_chunk());
                    return Some(full_chunk);
                }
            }
            None
        });

        chunk_iterator.par_bridge().for_each_init(
            || {
                (
                    Vec::<u8>::with_capacity(64 * 1024),
                    ThreadBuffers::new(),
                )
            },
            |(text_buffer, thread_bufs), chunk| {
                //let (text_buffer, encoded_fwd, encoded_rev, thread_bufs) = &mut *bufs;
                compute_chunk(chunk, text_buffer, thread_bufs, &y_index, &writer, &pb_main_clone, &pb_seqs, is_single_file, &recycle_tx);

                if !text_buffer.is_empty() {
                    writer.lock().unwrap().write_all(text_buffer).unwrap();
                    text_buffer.clear();
                }
            }
        );

        if !chunk.bounds.is_empty() {
            let bytes = chunk.memory.len();
            let mut lock = global_remainder.lock().unwrap();
            lock.0 += bytes;
            lock.1.push(chunk);

            if lock.0 >= MAX_BYTES {
                let chunks_to_process = std::mem::take(&mut lock.1);
                lock.0 = 0;
                drop(lock);

                let mut text_buffer = Vec::<u8>::with_capacity(64 * 1024);
                let mut thread_bufs = ThreadBuffers::new();

                for c in chunks_to_process {
                    compute_chunk(c, &mut text_buffer, &mut thread_bufs, &y_index, &writer, &pb_main_clone, &pb_seqs, is_single_file, &recycle_tx);
                }
                if !text_buffer.is_empty() {
                    writer.lock().unwrap().write_all(&text_buffer).unwrap();
                }
            }
        }

        if !is_single_file { pb_main_clone.inc(1); }
    });

    let final_chunks = std::mem::take(&mut global_remainder.lock().unwrap().1);
    if !final_chunks.is_empty() {
        let mut text_buffer = Vec::<u8>::with_capacity(64 * 1024);
        let mut thread_bufs = ThreadBuffers::new();

        for c in final_chunks {
            compute_chunk(c, &mut text_buffer, &mut thread_bufs, &y_index, &writer, &pb_main, &pb_seqs, is_single_file, &recycle_tx);
        }
        if !text_buffer.is_empty() { writer.lock().unwrap().write_all(&text_buffer).unwrap(); }
    }

    writer.lock().unwrap().flush().unwrap();
    pb_main.finish_with_message(format!("All files read."));
    pb_seqs.finish_with_message(format!("Processing complete! Total sequences: {}", pb_seqs.position()));

    //#[allow(unreachable_code)]
    unsafe { libc::_exit(0); }
}
