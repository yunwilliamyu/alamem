use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

#[derive(Debug, Clone)]
struct HitRecord {
    query: String,
    target: String,
    _q_size: usize,
    t_size: usize,
    q_start: usize,
    q_end: usize,
    t_start: usize,
    t_end: usize,
    strand: char,
    ani: f64,
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <input_hits.tsv> [output_regions.tsv]", args[0]);
        std::process::exit(1);
    }

    let input_path = &args[1];
    let file = File::open(input_path).expect("Failed to open input file");
    let reader = BufReader::new(file);

    let mut target_hits: HashMap<String, Vec<HitRecord>> = HashMap::new();

    // 1. Read input hits and apply initial coverage/length filters
    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();

        // Skip headers or empty lines
        if trimmed.is_empty() || trimmed.starts_with("Query") || trimmed.starts_with('#') {
            continue;
        }

        if let Some(hit) = parse_hit_line(trimmed) {
            target_hits.entry(hit.target.clone()).or_default().push(hit);
        }
    }

    // Setup output writer
    let writer_box: Box<dyn Write> = if args.len() >= 3 {
        Box::new(BufWriter::new(File::create(&args[2])?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };
    let mut writer = writer_box;

    writeln!(
        writer,
        "Target\tT_Size\tT_Start\tT_End\tRegion_Size\tTop1_Query\tTop1_Q_Start\tTop1_Q_End\tTop1_Strand\tTop1_ANI\tTop2_Query\tTop2_Q_Start\tTop2_Q_End\tTop2_Strand\tTop2_ANI"
    )?;

    let mut total_regions_found = 0;
    let mut sorted_targets: Vec<_> = target_hits.into_iter().collect();
    sorted_targets.sort_unstable_by(|a, b| a.0.cmp(&b.0));

    // 2. Process targets independently
    for (t_name, mut hits) in sorted_targets {
        // Sort hits by target start to optimize pairwise iteration
        hits.sort_unstable_by_key(|h| h.t_start);

        let mut putative_intersections: Vec<(usize, usize)> = Vec::new();

        // STEP A: Compute all pairwise intersections from distinct query species
        for i in 0..hits.len() {
            for j in (i + 1)..hits.len() {
                // Early exit: hit J starts after hit I ends
                if hits[j].t_start >= hits[i].t_end {
                    break;
                }

                // Require distinct query species
                if hits[i].query != hits[j].query {
                    let start = hits[i].t_start.max(hits[j].t_start);
                    let end = hits[i].t_end.min(hits[j].t_end);

                    if end > start {
                        let size = end - start;
                        if size >= 32 {
                            putative_intersections.push((start, end));
                        }
                    }
                }
            }
        }

        if putative_intersections.is_empty() {
            continue;
        }

        // STEP B: Deduplicate candidates
        putative_intersections.sort_unstable_by(|a, b| {
            a.0.cmp(&b.0).then_with(|| b.1.cmp(&a.1))
        });
        putative_intersections.dedup();

        // STEP C: Containment Filter (remove regions entirely inside another)
        let mut filtered_regions: Vec<(usize, usize)> = Vec::new();

        for &(cand_start, cand_end) in &putative_intersections {
            let mut is_contained = false;

            for &(other_start, other_end) in &putative_intersections {
                if cand_start == other_start && cand_end == other_end {
                    continue;
                }

                // If 'other' completely encompasses 'cand', discard 'cand'
                if other_start <= cand_start && other_end >= cand_end {
                    is_contained = true;
                    break;
                }
            }

            if !is_contained {
                filtered_regions.push((cand_start, cand_end));
            }
        }

        // STEP D: Output surviving maximal putative intersections with top-2 ANI species info
        let t_size = hits.first().map(|h| h.t_size).unwrap_or(0);

        for (r_start, r_end) in filtered_regions {
            let size = r_end - r_start;

            // Gather all hits fully spanning this intersection region
            let mut covering_hits: Vec<&HitRecord> = hits
            .iter()
            .filter(|h| h.t_start <= r_start && h.t_end >= r_end)
            .collect();

            // Sort hits by ANI descending
            covering_hits.sort_unstable_by(|a, b| {
                b.ani.partial_cmp(&a.ani).unwrap_or(std::cmp::Ordering::Equal)
            });

            // Extract highest-ANI hits for distinct query species
            let mut top_hits: Vec<&HitRecord> = Vec::new();
            let mut seen_queries = HashSet::new();

            for hit in covering_hits {
                if seen_queries.insert(&hit.query) {
                    top_hits.push(hit);
                    if top_hits.len() == 2 {
                        break;
                    }
                }
            }

            let format_hit = |h_opt: Option<&&HitRecord>| -> (String, String, String, String, String) {
                match h_opt {
                    Some(h) => (
                        h.query.clone(),
                                h.q_start.to_string(),
                                h.q_end.to_string(),
                                h.strand.to_string(),
                                format!("{:.2}", h.ani),
                    ),
                    None => (
                        "NA".to_string(),
                             "NA".to_string(),
                             "NA".to_string(),
                             "NA".to_string(),
                             "NA".to_string(),
                    ),
                }
            };

            let (h1_q, h1_qs, h1_qe, h1_str, h1_ani) = format_hit(top_hits.get(0));
            let (h2_q, h2_qs, h2_qe, h2_str, h2_ani) = format_hit(top_hits.get(1));

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                t_name,
                t_size,
                r_start,
                r_end,
                size,
                h1_q,
                h1_qs,
                h1_qe,
                h1_str,
                h1_ani,
                h2_q,
                h2_qs,
                h2_qe,
                h2_str,
                h2_ani
            )?;
            total_regions_found += 1;
        }
    }

    eprintln!(
        "Done. Output {} discrete non-contained regions.",
        total_regions_found
    );
    Ok(())
}

fn parse_hit_line(line: &str) -> Option<HitRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    // TSV format: Query\tTarget\tQ_Size\tQ_Start\tQ_End\tT_Size\tT_Start\tT_End\tStrand\tANI\tScore
    if fields.len() < 10 {
        return None;
    }

    let query = fields[0].to_string();
    let target = fields[1].to_string();
    let q_size: usize = fields[2].parse().ok().unwrap_or(0);
    let q_start: usize = fields[3].parse().ok().unwrap_or(0);
    let q_end: usize = fields[4].parse().ok().unwrap_or(0);
    let t_size: usize = fields[5].parse().ok()?;
    let mut t_start: usize = fields[6].parse().ok()?;
    let mut t_end: usize = fields[7].parse().ok()?;
    let strand = fields[8].chars().next().unwrap_or('+');
    let ani: f64 = fields[9].parse().ok().unwrap_or(0.0);

    if t_start > t_end {
        std::mem::swap(&mut t_start, &mut t_end);
    }

    let hit_len = t_end - t_start;

    // Filter out tiny noise (< 32 bp)
    if hit_len < 32 {
        return None;
    }

    // Filter out hits covering > 50% of the target sequence
    if t_size > 0 && hit_len > (t_size / 2) {
        return None;
    }

    // Filter out hits covering > 50% of the query sequence (if size available)
    if q_size > 0 && hit_len > (q_size / 2) {
        return None;
    }

    Some(HitRecord {
        query,
         target,
         _q_size: q_size,
         t_size,
         q_start,
         q_end,
         t_start,
         t_end,
         strand,
         ani,
    })
}
