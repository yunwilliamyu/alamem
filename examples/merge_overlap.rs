use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();

    // Read from file if provided as argument, otherwise read from standard input
    let reader: Box<dyn BufRead> = if args.len() > 1 {
        Box::new(BufReader::new(File::open(&args[1])?))
    } else {
        Box::new(BufReader::new(io::stdin()))
    };

    let mut intervals: Vec<(usize, usize)> = Vec::new();

    // 1. Parse input coordinates
    for line_res in reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();

        // Skip empty lines or comment lines
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() >= 2 {
            // Attempt to parse start and end as integers
            match (fields[0].parse::<usize>(), fields[1].parse::<usize>()) {
                (Ok(mut start), Ok(mut end)) => {
                    if start > end {
                        std::mem::swap(&mut start, &mut end);
                    }
                    intervals.push((start, end));
                }
                _ => {
                    // Non-numeric fields (e.g. "Start End" header) are skipped
                    continue;
                }
            }
        }
    }

    if intervals.is_empty() {
        return Ok(());
    }

    // 2. Sort intervals by start position
    intervals.sort_unstable_by_key(|&(start, _)| start);

    // 3. Merge overlapping and adjacent intervals
    let mut merged: Vec<(usize, usize)> = Vec::with_capacity(intervals.len());

    for (start, end) in intervals {
        if let Some(last) = merged.last_mut() {
            // Extend previous interval if current interval overlaps or touches it
            if start <= last.1 {
                last.1 = last.1.max(end);
            } else {
                merged.push((start, end));
            }
        } else {
            merged.push((start, end));
        }
    }

    // 4. Output merged intervals
    let mut stdout = io::BufWriter::new(io::stdout());
    for (start, end) in merged {
        writeln!(stdout, "{}\t{}", start, end)?;
    }

    Ok(())
}
