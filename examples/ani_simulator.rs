use alamem::{
    filter_overlapping_hits, process_query_sequence, reverse_complement, AlignConfig, ThreadBuffers,
    YIndex, CONFIG,
};
use std::time::Instant;
use clap::Parser;

struct SimpleRng {
    seed: u64,
}

impl SimpleRng {
    fn next_u64(&mut self) -> u64 {
        self.seed = self.seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.seed
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[derive(Parser, Debug)]
#[command(name = "ani_simulator")]
struct AniSimCli {
    pub target_fasta: String,
    pub num_reads: usize,
    pub read_len: usize,
    pub target_ani: f64,

    #[command(flatten)]
    pub align_config: AlignConfig,
}

fn main() {
    let cli = AniSimCli::parse();
    println!("{:#?}", cli);
    CONFIG.set(cli.align_config).expect("Config already initialized");

    let target_file = cli.target_fasta;
    let num_reads = cli.num_reads;
    let read_len = cli.read_len;
    let target_ani = cli.target_ani;

    let start_idx = Instant::now();
    let y_index = YIndex::build(&[target_file.clone()]);
    println!("Index built in {:.2?}", start_idx.elapsed());

    if y_index.seq_data.is_empty() {
        eprintln!("Error: Target file contains no valid sequences.");
        std::process::exit(1);
    }

    let mut rng = SimpleRng { seed: 987654321 };

    let mut total_error: f64 = 0.0;
    let mut total_abs_error: f64 = 0.0;
    let mut total_sq_error: f64 = 0.0; // Track sum of squared errors for Standard Error
    let mut successfully_mapped = 0;
    let mut total_len: usize = 0;

    let error_rate = 1.0 - (target_ani / 100.0);

    println!("\nSimulating {} reads (Len: {}, Target ANI: {}%)...", num_reads, read_len, target_ani);

    let mut thread_bufs = ThreadBuffers::new();
    for i in 0..num_reads {
        // 1. Pick a random sequence and position
        let seq_idx = (rng.next_u64() as usize) % y_index.seq_data.len();
        let seq = &y_index.seq_data[seq_idx];
        if seq.len() < read_len { continue; }

        let pos = (rng.next_u64() as usize) % (seq.len() - read_len + 1);
        let original_subseq = &seq[pos..pos + read_len];

        // 2. Mutate sequence and build a cumulative history for Local ANI tracking
        //let mut mutated_fwd = Vec::with_capacity(read_len + (read_len / 10));
        thread_bufs.encoded_fwd.clear();
        thread_bufs.encoded_rev.clear();

        let mut cum_matches = Vec::with_capacity(read_len + (read_len / 10) + 1);
        let mut cum_ref = Vec::with_capacity(read_len + (read_len / 10) + 1);

        cum_matches.push(0);
        cum_ref.push(0);

        let mut current_matches = 0;
        let mut current_ref = 0;

        for &base in original_subseq {
            if base > 3 { // N or other ambiguous base
                thread_bufs.encoded_fwd.push(base);
                current_ref += 1;
                cum_matches.push(current_matches);
                cum_ref.push(current_ref);
                continue;
            }

            if rng.next_f64() < error_rate {
                let mut_type = rng.next_f64();
                if mut_type < 0.70 {
                    // SNP (Substitution)
                    let mut new_base = (rng.next_u64() % 4) as u8;
                    if new_base == base { new_base = (base + 1) % 4; }
                    thread_bufs.encoded_fwd.push(new_base);
                    current_ref += 1;
                    cum_matches.push(current_matches);
                    cum_ref.push(current_ref);
                } else if mut_type < 0.85 {
                    // Deletion
                    current_ref += 1;
                } else {
                    // Insertion
                    thread_bufs.encoded_fwd.push((rng.next_u64() % 4) as u8);
                    cum_matches.push(current_matches);
                    cum_ref.push(current_ref);

                    thread_bufs.encoded_fwd.push(base);
                    current_matches += 1;
                    current_ref += 1;
                    cum_matches.push(current_matches);
                    cum_ref.push(current_ref);
                }
            } else {
                // Match
                thread_bufs.encoded_fwd.push(base);
                current_matches += 1;
                current_ref += 1;
                cum_matches.push(current_matches);
                cum_ref.push(current_ref);
            }
        }

        // 3. Prepare reverse complement
        reverse_complement(&thread_bufs.encoded_fwd, &mut thread_bufs.encoded_rev);

        // 4. Align
        process_query_sequence(&y_index, &mut thread_bufs);
        filter_overlapping_hits(&mut thread_bufs.hits);

        // 5. Evaluate Accuracy against True Local Origin
        let mut longest_hit_info = None;
        let mut max_hit_len = 0;

        for hit in &thread_bufs.hits {
            if hit.t_id == seq_idx {
                let fwd_q_start = hit.q_start;
                let fwd_q_end = hit.q_end;

                let overlap_start = hit.t_start.max(pos);
                let overlap_end = hit.t_end.min(pos + current_ref);

                if overlap_end > overlap_start {
                    let hit_len = overlap_end - overlap_start;

                    if hit_len > max_hit_len {
                        max_hit_len = hit_len;

                        let local_matches = cum_matches[fwd_q_end] - cum_matches[fwd_q_start];
                        let local_ref_consumed = cum_ref[fwd_q_end] - cum_ref[fwd_q_start];
                        let local_q_span = fwd_q_end - fwd_q_start;

                        let local_max_span = local_q_span.max(local_ref_consumed);
                        let true_local_ani = (local_matches as f64 / local_max_span as f64) * 100.0;

                        longest_hit_info = Some((hit.ani, true_local_ani, hit_len));
                    }
                }
            }
        }

        // Apply statistics using only the longest hit found for this read
        if let Some((est_ani, true_local_ani, hit_len)) = longest_hit_info {
            successfully_mapped += 1;

            let error = est_ani - true_local_ani;
            total_error += error;
            total_abs_error += error.abs();
            total_sq_error += error * error;
            total_len += hit_len;

            // Print the first few for visual inspection
            if successfully_mapped <= 5 {
                println!(
                    "Read {:>4} | True Local ANI: {:.2}% | Est ANI: {:.2}% | Error: {:.2}% | Len: {:.0}",
                    i + 1, true_local_ani, est_ani, error, hit_len
                );
            }
        }
    }

    // Final Statistics
    println!("\n=== ANI Estimation Results ===");
    println!("Target ANI specified: {}%", target_ani);
    println!("Successfully mapped:  {}/{} reads", successfully_mapped, num_reads);

    if successfully_mapped > 0 {
        let n = successfully_mapped as f64;
        let mean_absolute_error = total_abs_error / n;
        let mean_bias = total_error / n;

        // Calculate Standard Error of the Mean: SE = s / sqrt(N)
        let standard_error = if successfully_mapped > 1 {
            let sample_variance = (total_sq_error - (total_error * total_error) / n) / (n - 1.0);
            (sample_variance.max(0.0) / n).sqrt()
        } else {
            0.0
        };

        println!("Mean Absolute Error:  {:.3}%", mean_absolute_error);
        println!("          Mean Bias:  {:.3}%", mean_bias);
        println!("     Standard Error:  {:.3}%", standard_error);
        println!("        Mean Length:  {:.0}", total_len as f64 / n);
    }
}
