// src/lib-old-working.rs
use rustc_hash::FxHashMap;
use needletail::{parse_fastx_file};
use clap::Parser;
use std::sync::OnceLock;
use rayon::prelude::*;


/// Validator function
fn parse_odd_kmer(s: &str) -> Result<usize, String> {
    let k: usize = s
    .parse()
    .map_err(|_| format!("'{}' is not a valid number", s))?;
    if !(9..=15).contains(&k) {
        return Err(format!("k-mer size must be between 9 and 15 (got {})", k));
    }
    if k % 2 == 0 {
        return Err(format!("k-mer size must be an odd number (got {})", k));
    }
    Ok(k)
}
#[derive(Parser, Debug, Clone)]
pub struct AlignConfig {
    /// K-mer size (odd integer between 9 and 15)
    #[arg(short = 'k', long, default_value_t = 13, value_parser = parse_odd_kmer)]
    pub kmer_size: usize,

    /// Stride for checking anchors [default: kmer_size]
    #[arg(short = 's', long)]
    pub stride: Option<usize>,

    /// ANI filter for output hits
    #[arg(
    long,
    default_value_t = 90.0,
    value_parser = |s: &str| -> Result<f64, String> {
        let val = s.parse::<f64>().map_err(|e| e.to_string())?;
        if (85.0..=100.0).contains(&val) {
            Ok(val)
        } else {
            Err(String::from("Value must be between 85.0 and 100.0"))
        }
    }
    )]
    pub min_ani: f64,

    /// Minimum length of output hits
    #[arg(short = 'l', long, default_value_t = 40)]
    pub min_len: usize,

    /// Filter out common seeds from index (shouldn't trigger)
    #[arg(short = 'm', long, default_value_t = 1000000)]
    pub max_seed_multiplicity: usize,

    /// Maximum gap for chaining
    #[arg(long, default_value_t = 100)]
    pub max_chain_gap: usize,
}

// The global, thread-safe container for your runtime configuration
pub static CONFIG: OnceLock<AlignConfig> = OnceLock::new();


/// Safe, cross-platform entry point
pub fn vectorize_base_lookup(seq: &[u8], out: &mut [u8]) {
    // 1. Try x86-64 AVX2
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                vectorize_base_lookup_avx2(seq, out);
                return;
            }
        }
    }

    // 2. Try ARM NEON (Apple Silicon, AWS Graviton)
    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            unsafe {
                vectorize_base_lookup_neon(seq, out);
                return;
            }
        }
    }

    // 3. Universal Scalar Fallback (WASM, older CPUs, RISC-V)
    vectorize_base_lookup_scalar(seq, out);
}

/// The universal fallback loop you already wrote
fn vectorize_base_lookup_scalar(seq: &[u8], out: &mut [u8]) {
    for j in 0..seq.len() {
        let b = seq[j] & 0x0F;
        out[j] = match b {
            1 => 0,       // A
            3 => 1,       // C
            4 | 5 => 3,   // T or U
            7 => 2,       // G
            _ => 128,     // N
        };
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn vectorize_base_lookup_avx2(seq: &[u8], out: &mut [u8]) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let mask = _mm256_set1_epi8(0x0F);
    let lut = _mm256_setr_epi8(
        -128, 0, -128, 1, 3, 3, -128, 2, -128, -128, -128, -128, -128, -128, -128, -128,
        -128, 0, -128, 1, 3, 3, -128, 2, -128, -128, -128, -128, -128, -128, -128, -128,
    );

    let mut i = 0;
    while i + 32 <= seq.len() {
        let ascii_chars = _mm256_loadu_si256(seq.as_ptr().add(i) as *const _);
        let nibbles = _mm256_and_si256(ascii_chars, mask);
        let translated = _mm256_shuffle_epi8(lut, nibbles);
        _mm256_storeu_si256(out.as_mut_ptr().add(i) as *mut _, translated);
        i += 32;
    }

    // Process remainder
    vectorize_base_lookup_scalar(&seq[i..], &mut out[i..]);
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn vectorize_base_lookup_neon(seq: &[u8], out: &mut [u8]) {
    use std::arch::aarch64::*;

    // 1. Bitmask for lower 4 bits
    let mask = vdupq_n_u8(0x0F);

    // 2. In-Register LUT (NEON uses unsigned u8, so 128 works perfectly)
    let lut_data: [u8; 16] = [
        128, 0, 128, 1, 3, 3, 128, 2, 128, 128, 128, 128, 128, 128, 128, 128
    ];
    let lut = vld1q_u8(lut_data.as_ptr());

    let mut i = 0;

    // Process 16 bases at a time
    while i + 16 <= seq.len() {
        let ascii_chars = vld1q_u8(seq.as_ptr().add(i));
        let nibbles = vandq_u8(ascii_chars, mask);

        // vqtbl1q_u8 is the ARM equivalent of _mm256_shuffle_epi8
        let translated = vqtbl1q_u8(lut, nibbles);

        vst1q_u8(out.as_mut_ptr().add(i), translated);
        i += 16;
    }

    // Process remainder
    vectorize_base_lookup_scalar(&seq[i..], &mut out[i..]);
}


/// Encodes an ASCII DNA sequence into 0, 1, 2, 3, or 128/131 (invalid).
/// Uses SIMD where available and falls back to scalar.
pub fn encode_sequence(seq: &[u8], b: &mut Vec<u8>) {
    if seq.is_empty() {
        b.clear();
        return;
    }

    b.clear();
    b.reserve(seq.len());
    // Set length safely since we overwrite it entirely immediately
    unsafe { b.set_len(seq.len()); }

    vectorize_base_lookup(seq, b);
}

/// Only works to reverse complement already encoded in 0, 1, 2, 3, 128/131, string
#[inline]
pub fn reverse_complement(encoded_fwd: &[u8], encoded_rev: &mut Vec<u8>) {
    encoded_rev.clear();
    encoded_rev.extend(encoded_fwd.iter().rev().map(|&byte| byte ^ 3));
}

#[derive(Clone, Copy, Debug)]
pub struct Hit {
    //    q_id: String,
    pub t_id: usize, // Index into y_index.seq_names
    pub q_size: usize,
    pub q_start: usize,
    pub q_end: usize,
    pub t_size: usize,
    pub t_start: usize,
    pub t_end: usize,
    pub strand: char,
    pub ani: f64,
    pub score: i32,
}


#[derive(Clone, Copy, Debug)]
pub struct TargetPos {
    seq_id: u32,
    pos: u32, // Use MSB to represent strand: 0 = fwd, 1 = rev. Rest are position. (max 2.14 billion bases, so won't work for lungfish)
}

#[derive(Clone, Copy, Debug)]
pub struct Anchor {
    t_id: u32,
    q_pos: u32,
    t_pos: u32,
    diag: i64,
}

#[derive(Clone, Copy, Debug)]
pub struct MEM {
    pub q_start: i32,
    pub t_start: i32,
    pub q_end: i32,
    pub t_end: i32,
    pub len: i32,
    pub diag: i64,
}


// Pooling buffers to completely eliminate millions of inner-loop allocations
#[repr(align(64))]
pub struct ThreadBuffers {
    pub encoded_fwd: Vec<u8>,
    pub encoded_rev: Vec<u8>,
    pub anchors_fwd: Vec<Anchor>,
    pub anchors_rev: Vec<Anchor>,
    pub mems: Vec<MEM>,
    pub chain: Vec<MEM>,
    pub scores: Vec<i32>,
    pub preds: Vec<usize>,
    pub order: Vec<usize>,
    pub used: Vec<bool>,
    pub hits: Vec<Hit>,
}

impl ThreadBuffers {
    pub fn new() -> Self {
        Self {
            encoded_fwd: Vec::with_capacity(1024 * 1024),
            encoded_rev: Vec::with_capacity(1024 * 1024),
            anchors_fwd: Vec::with_capacity(1024),
            anchors_rev: Vec::with_capacity(1024),
            mems: Vec::with_capacity(1024),
            chain: Vec::with_capacity(1024),
            scores: Vec::with_capacity(1024),
            preds: Vec::with_capacity(1024),
            order: Vec::with_capacity(1024),
            used: Vec::with_capacity(1024),
            hits: Vec::with_capacity(128),
        }
    }
}


pub struct YIndex {
    pub kmer_map: FxHashMap<u32, (usize, u32)>,
    pub flat_kmers: Vec<TargetPos>,
    pub presence_bits: Vec<u64>,
    pub seq_data: Vec<Vec<u8>>,
    pub seq_names: Vec<Vec<u8>>,
}

impl YIndex {
    pub fn build(y_files: &[String]) -> Self {
        let config = CONFIG.get().expect("Config not initialized");

        let mut seq_data: Vec<Vec<u8>> = Vec::new();
        let mut seq_names: Vec<Vec<u8>> = Vec::new();
        let mut bins: Vec<Vec<(u32, TargetPos)>> = vec![Vec::new(); 256];
        let mut seq_id: u32 = 0;

        let mask: u32 = ((1_u64 << (2 * config.kmer_size)) - 1) as u32;
        let k_shift = 2 * (config.kmer_size - 1);

        let bin_shift = (2 * config.kmer_size) - 8;

        // Pre-allocate buffer to prevent heap thrashing inside the loop
        let mut encoded_seq = Vec::new();

        for file in y_files {
            let mut reader = parse_fastx_file(file).expect("Failed to open Y file");

            while let Some(record) = reader.next() {
                let seq = record.expect("Invalid sequence");
                let raw_seq = seq.seq();

                encoded_seq.clear();
                encoded_seq.reserve(raw_seq.len());
                encode_sequence(&raw_seq, &mut encoded_seq);

                let id_slice = seq.id();
                let id_len = memchr::memchr2(b' ', b'\t', id_slice).unwrap_or(id_slice.len());
                seq_names.push(id_slice[..id_len].to_vec());

                if encoded_seq.len() >= config.kmer_size {
                    let mut h_fwd: u32 = 0;
                    let mut h_rev: u32 = 0;
                    let mut valid: usize = 0;

                    for (i, &val) in encoded_seq.iter().enumerate() {
                        if val > 3 {
                            h_fwd = 0;
                            h_rev = 0;
                            valid = 0;
                            continue;
                        }

                        h_fwd = ((h_fwd << 2) | val as u32) & mask;
                        let rev_val = (val as u32) ^ 3;
                        h_rev = (h_rev >> 2) | (rev_val << k_shift);
                        valid += 1;

                        if valid >= config.kmer_size {
                            let canonical = h_fwd.min(h_rev);
                            let is_rev = h_fwd != canonical;

                            let pos: u32 = (i + 1 - config.kmer_size) as u32;
                            let encoded_pos = if is_rev { pos | 0x8000_0000 } else { pos };

                            let bin_idx = (canonical >> bin_shift) as usize;

                            bins[bin_idx].push((
                                canonical,
                                TargetPos { seq_id, pos: encoded_pos },
                            ));
                        }
                    }
                }
                seq_data.push(encoded_seq.clone());
                seq_id += 1;
            }
        }

        let filter_size_u64s = 65_536;
        let filter_mask = 0x3FFFFF;

        // 2. Compaction and Map construction
        //for mut current_bin in bins.into_iter() {
        let processed_bins: Vec<_> = bins.into_par_iter().map(|mut current_bin| {
            //current_bin.sort_unstable_by_key(|&(k, target_pos)| (k, target_pos.seq_id, target_pos.pos));
            current_bin.sort_unstable_by_key(|&(k, target_pos)| {
                let raw_pos = target_pos.pos & 0x7FFF_FFFF;
                let strand_bit = target_pos.pos & 0x8000_0000;
                (k, target_pos.seq_id, raw_pos, strand_bit)
            });
            /*current_bin.sort_unstable_by_key(|&(k, target_pos)| {
                (k, target_pos.seq_id, target_pos.pos & 0x7FFF_FFFF)
            });*/
            //current_bin.sort_unstable_by_key(|&(k, _)| k);

            let mut local_flat = Vec::new();
            let mut local_map_entries = Vec::new();
            let mut local_presence = vec![0_u64; filter_size_u64s];
            let mut kept = 0;
            let mut masked = 0;

            let mut i = 0;
            while i < current_bin.len() {
                let key = current_bin[i].0;

                let mut j = i + 1;
                while j < current_bin.len() && current_bin[j].0 == key {
                    j += 1;
                }
                let count = j - i;

                if count <= config.max_seed_multiplicity as usize {
                    let local_start = local_flat.len() as usize;

                    // Extract payload
                    for read_idx in i..j {
                        local_flat.push(current_bin[read_idx].1);
                    }

                    local_map_entries.push((key, local_start, count as u32));

                    // Update Bloom filter
                    let mapped_h = (key ^ (key >> 11) ^ (key >> 22)) & filter_mask;
                    local_presence[(mapped_h as usize) >> 6] |= 1_u64 << (mapped_h & 63);

                    kept += 1;
                } else {
                    masked += 1;
                }

                i = j;
            }
            (local_flat, local_map_entries, local_presence, kept, masked)
        }).collect();

        let mut presence_bits = vec![0_u64; filter_size_u64s];
        let total_flat_len: usize = processed_bins.iter().map(|b| b.0.len()).sum();
        let total_kept_keys: usize = processed_bins.iter().map(|b| b.3).sum();
        let mut flat_kmers: Vec<TargetPos> = Vec::with_capacity(total_flat_len);
        let mut kmer_map: FxHashMap<u32, (usize, u32)> = FxHashMap::with_capacity_and_hasher(total_kept_keys, Default::default());

        let mut global_kept = 0;
        let mut global_masked = 0;
        let mut global_offset: usize = 0;

        for (local_flat, local_map_entries, local_presence, kept, masked) in processed_bins {
            flat_kmers.extend(local_flat);

            for (key, local_start, count) in local_map_entries {
                kmer_map.insert(key, (global_offset + local_start, count));
            }

            // The blazing fast SIMD Bitwise OR
            for (global_word, local_word) in presence_bits.iter_mut().zip(local_presence) {
                *global_word |= local_word;
            }

            global_kept += kept;
            global_masked += masked;
            global_offset = flat_kmers.len();
        }

        if global_masked > 0 {
            eprintln!(
                "Masked {} overrepresented k-mers (threshold: {}). Kept {} unique valid k-mers.",
                      global_masked, config.max_seed_multiplicity, global_kept
            );
        }

        if global_masked > global_kept {
            panic!(
                "Fatal Error: Too many k-mers were masked! (Masked: {}, Kept: {}). \
The index would be unsearchable. Please increase `max_seed_multiplicity` \
or increase `kmer_size`.",
global_masked, global_kept
            );
        }

        YIndex {
            kmer_map,
            flat_kmers,
            presence_bits,
            seq_data,
            seq_names,
        }
    }

    /// Fast lookup returning exactly the slice of targets, utilizing the Bloom filter.
    #[inline(always)]
    pub fn get_kmer(&self, query_kmer: u32) -> &[TargetPos] {
        let filter_mask = 0x3FFFFF;
        let mapped_h = (query_kmer ^ (query_kmer >> 11) ^ (query_kmer >> 22)) & filter_mask;

        let bit_idx = (mapped_h as usize) >> 6;
        let bit_mask = 1_u64 << (mapped_h & 63);

        // Fast Path: Bloom filter rejection
        if (self.presence_bits[bit_idx] & bit_mask) == 0 {
            return &[];
        }

        // Slow Path: Hash map resolution
        if let Some(&(start, count)) = self.kmer_map.get(&query_kmer) {
            //let start = start as usize;
            let count = count as usize;
            &self.flat_kmers[start..(start + count)]
        } else {
            &[]
        }
    }
}


#[inline(always)]
fn rev_comp_u32(h_fwd: u32, k: usize) -> u32 {
    let mut x = h_fwd;
    x = ((x >> 2) & 0x3333_3333) | ((x & 0x3333_3333) << 2);
    x = ((x >> 4) & 0x0F0F_0F0F) | ((x & 0x0F0F_0F0F) << 4);
    x = x.swap_bytes();
    (!x) >> (32 - (2 * k))
}
/// Validates all k bases in a single pass using two overlapping 64-bit reads.
/// Safe and sound for any K between 8 and 16.
#[inline(always)]
fn is_invalid_k8_16<const K: usize>(chunk: &[u8]) -> bool {
    debug_assert!(K >= 8 && K <= 16);
    unsafe {
        let ptr = chunk.as_ptr();
        let b_head = std::ptr::read_unaligned(ptr as *const u64);
        let b_tail = std::ptr::read_unaligned(ptr.add(K - 8) as *const u64);
        ((b_head | b_tail) & 0xFCFC_FCFC_FCFC_FCFC) != 0
    }
}

/// Builds the forward hash with parallel tree-splitting, unrolled at compile-time.
#[inline(always)]
fn hash_k_fast<const K: usize>(chunk: &[u8]) -> u32 {
    let mid = (K + 1) / 2;
    let mut h1 = 0u32;
    for &v in &chunk[..mid] {
        h1 = (h1 << 2) | (v as u32);
    }
    let mut h2 = 0u32;
    for &v in &chunk[mid..K] {
        h2 = (h2 << 2) | (v as u32);
    }
    (h1 << ((K - mid) * 2)) | h2
}

/// Identifies anchors by evaluating KMER_SIZE chunks of the query sequence at stride STRIDE
pub fn get_anchors(
    q_seq: &[u8],
    y_index: &YIndex,
    bufs: &mut ThreadBuffers,
) {
    let config = CONFIG.get().expect("Config not initialized");
    let q_len = q_seq.len();
    let k = config.kmer_size;
    if q_len < k { return; }

    let stride = config.stride.unwrap_or(k);

    // Single generic runner: Rust monomorphizes this for each specialized closure
    #[inline(always)]
    fn run<F>(
        q_seq: &[u8],
        y_index: &YIndex,
        bufs: &mut ThreadBuffers,
        k: usize,
        stride: usize,
        mut get_hash: F,
    ) where
    F: FnMut(&[u8]) -> Option<u32>,
    {
        let q_len = q_seq.len();
        let mut i = 0;

        while i <= q_len - k {

            let chunk = unsafe { q_seq.get_unchecked(i..i + k) };

            let Some(h_fwd) = get_hash(chunk) else {
                i += stride;
                continue;
            };

            let h_rev = rev_comp_u32(h_fwd, k);
            let canonical = h_fwd.min(h_rev);

            let t_positions = y_index.get_kmer(canonical);

            if t_positions.is_empty() {
                i += stride;
                continue;
            }

            let q_is_rev = h_fwd != canonical;
            let q_pos = i as u32;
            let rc_pos = (q_len - k - i) as u32;

            for t_pos in t_positions {
                let t_is_rev = (t_pos.pos & 0x8000_0000) != 0;
                let actual_t_pos = t_pos.pos & 0x7FFF_FFFF;

                let (target_buf, active_q_pos) = if q_is_rev == t_is_rev {
                    (&mut bufs.anchors_fwd, q_pos)
                } else {
                    (&mut bufs.anchors_rev, rc_pos)
                };

                target_buf.push(Anchor {
                    t_id: t_pos.seq_id,
                    q_pos: active_q_pos,
                    t_pos: actual_t_pos,
                    diag: actual_t_pos as i64 - active_q_pos as i64,
                });
            }


            i += stride;
        }
    }

    match (k, stride) {
        // HOT PATHS
        (9, _)   => run(q_seq, y_index, bufs, 9, stride,   |c| if is_invalid_k8_16::<9>(c)  { None } else { Some(hash_k_fast::<9>(c)) }),
        (11, _) => run(q_seq, y_index, bufs, 11, stride, |c| if is_invalid_k8_16::<11>(c) { None } else { Some(hash_k_fast::<11>(c)) }),
        (13, _) => run(q_seq, y_index, bufs, 13, stride, |c| if is_invalid_k8_16::<13>(c) { None } else { Some(hash_k_fast::<13>(c)) }),
        (15, _)  => run(q_seq, y_index, bufs, 15, stride,  |c| if is_invalid_k8_16::<15>(c) { None }  else { Some(hash_k_fast::<15>(c)) }),

        // GENERIC FALLBACK
        _ => run(q_seq, y_index, bufs, k, stride, |c| {
            let mut h = 0u32;
            for &val in c {
                if val > 3 { return None; }
                h = (h << 2) | (val as u32);
            }
            Some(h)
        }),
    }
}

pub fn align_strand(
    q_seq: &[u8],
    strand: char,
    y_index: &YIndex,
    bufs: &mut ThreadBuffers,
) {
    let config = CONFIG.get().expect("Config not initialized");
    let q_size: usize = q_seq.len();

    let is_fwd = strand == '+';

    let mut anchors = if is_fwd {
        std::mem::take(&mut bufs.anchors_fwd)
    } else {
        std::mem::take(&mut bufs.anchors_rev)
    };

    if anchors.is_empty() {
        return;
    }

    anchors.sort_unstable_by_key(|a| (a.t_id, a.diag, a.q_pos));

    let mut start = 0;
    while start < anchors.len() {
        let t_id = anchors[start].t_id as usize;
        let mut end = start + 1;

        while end < anchors.len() && anchors[end].t_id as usize == t_id {
            end += 1;
        }

        let t_anchors = &anchors[start..end];
        start = end; // Go ahead and set start to end for next run since we have consumed it
        let t_seq: &[u8] = &y_index.seq_data[t_id];
        let t_size: usize = t_seq.len();

        bufs.mems.clear();
        let mut i: usize = 0;

        while i < t_anchors.len() {
            let q_seed_start = t_anchors[i].q_pos as usize;
            let t_seed_start = t_anchors[i].t_pos as usize;

            debug_assert!(q_seed_start + config.kmer_size <= q_size);
            debug_assert!(t_seed_start + config.kmer_size <= t_size);

            let left_match_len = q_seq[..q_seed_start].iter().rev()
            .zip(t_seq[..t_seed_start].iter().rev())
            .take_while(|(&q, &t)| q == t && q <= 3)
            .count();

            let right_match_len = q_seq[q_seed_start + config.kmer_size..].iter()
            .zip(t_seq[t_seed_start + config.kmer_size..].iter())
            .take_while(|(&q, &t)| q == t && q <= 3)
            .count();

            let q_start = (q_seed_start - left_match_len) as i32;
            let t_start = (t_seed_start - left_match_len) as i32;
            let len = (left_match_len + config.kmer_size + right_match_len) as i32;
            let q_end = q_start + len;
            let t_end = t_start + len;
            let diag = t_start as i64 - q_start as i64;
            let q_r = q_seed_start + config.kmer_size + right_match_len;

            bufs.mems.push(MEM {
                q_start,
                t_start,
                q_end,
                t_end,
                len,
                diag,
            });

            // 3. Skip subsequent anchors on same diag covered by extension
            let diag = t_anchors[i].t_pos as i64 - t_anchors[i].q_pos as i64;
            let mut j = i+1;
            while j < t_anchors.len() {
                let next_diag = t_anchors[j].t_pos as i64 - t_anchors[j].q_pos as i64;
                if next_diag != diag { break; }

                if (t_anchors[j].q_pos as usize + config.kmer_size) <= q_r {
                    j += 1;
                } else {
                    break;
                }
            }

            i = j;
        }

        //bufs.mems.sort_unstable_by_key(|h: &MEM| (h.t_start, h.q_start));
        bufs.mems.sort_unstable_by_key(|h: &MEM| (
            h.t_start + h.q_start,
            h.t_start.max(h.q_start),
            h.t_start.min(h.q_start)
        ));

        bufs.mems.dedup_by(|a: &mut MEM, b: &mut MEM| a.q_start == b.q_start && a.t_start == b.t_start);

        if bufs.mems.is_empty() { continue; }

        let n: usize = bufs.mems.len();
        bufs.scores.clear(); bufs.scores.extend(bufs.mems.iter().map(|h| h.len as i32));
        bufs.preds.clear(); bufs.preds.resize(n, usize::MAX);

        let max_gap = config.max_chain_gap as i32;
        let expected_ani = config.min_ani / 100.0;
        // Fast-math setup: Convert ANI percentage to a /256 fixed-point multiplier
        let ani_fixed_mult = (expected_ani * 256.0) as i32;

        let mut window_start = 0;
        for i_chain in 0..n {
            let h_i = bufs.mems[i_chain];

            while window_start < i_chain {
                let prev_mem = &bufs.mems[window_start];
                let sum_i = h_i.t_start + h_i.q_start;
                let sum_prev_end = prev_mem.t_end + prev_mem.q_end;

                if sum_i.saturating_sub(sum_prev_end) > (max_gap * 2) {
                    window_start += 1;
                } else {
                    break;
                }
            }

            for j_chain in (window_start..i_chain).rev() {
                let h_j = bufs.mems[j_chain];

                // Enforce collinearity
                if h_i.q_start < h_j.q_start || h_i.t_start < h_j.t_start { continue; }

                let true_q_gap = (h_i.q_start - h_j.q_end).max(0);
                let true_t_gap = (h_i.t_start - h_j.t_end).max(0);

                // Enforce max_gap
                if true_q_gap > max_gap || true_t_gap > max_gap {
                    continue;
                }

                let q_added = h_i.q_end - h_j.q_end.max(h_i.q_start);
                let t_added = h_i.t_end - h_j.t_end.max(h_i.t_start);
                let added_match = q_added.min(t_added).max(0);

                if added_match == 0 {
                    continue;
                }

                let diag_diff = (h_i.diag - h_j.diag).abs().min(max_gap as i64) as i32;
                let unaligned_gap = if diag_diff <= 2 && (true_q_gap == 0 || true_t_gap == 0) {
                    0
                } else {
                    diag_diff
                };

                let aligned_gap = true_q_gap.min(true_t_gap);

                // Fixed-point integer math (>> 8 is roughly equivalent to / 256)
                let gap_matches = if aligned_gap > 2 {
                    //((aligned_gap - 2) as f64 * expected_ani) as i32
                    ((aligned_gap - 2) * ani_fixed_mult) >> 8
                } else {
                    0
                };

                let chain_score = bufs.scores[j_chain] + added_match + gap_matches - unaligned_gap;
                if chain_score > bufs.scores[i_chain] {
                    bufs.scores[i_chain] = chain_score;
                    bufs.preds[i_chain] = j_chain;
                }

            }
        }

        bufs.order.clear(); bufs.order.extend(0..n);
        bufs.order.sort_unstable_by_key(|&idx| std::cmp::Reverse(bufs.scores[idx]));
        bufs.used.clear(); bufs.used.resize(n, false);

        for &idx in &bufs.order {
            if bufs.used[idx] || bufs.scores[idx] < config.min_len as i32 { continue; }

            let mut curr: usize = idx;
            bufs.chain.clear();

            while curr != usize::MAX {
                bufs.chain.push(bufs.mems[curr].clone());
                bufs.used[curr] = true;
                curr = bufs.preds[curr];
            }

            bufs.chain.reverse();

            let first = bufs.chain[0];
            let last = bufs.chain[bufs.chain.len() - 1];

            // Calculate average mem length to get regional ANI of chain
            let total_mem_len: usize = bufs.chain.iter().map(|m| m.len as usize).sum();
            let avg_mem_len: f64 = total_mem_len as f64 / bufs.chain.len() as f64;
            let true_mean_len: f64 = (avg_mem_len - config.kmer_size as f64).max(1.0);
            //let regional_expected_ani: f64 = 1.0 - (1.0 / true_mean_len);
            let regional_expected_ani: f64 = true_mean_len / (true_mean_len + 1.0);

            let final_q_start: usize = first.q_start as usize;
            let final_q_end: usize = last.q_end as usize;
            let final_t_start: usize = first.t_start as usize;
            let final_t_end: usize = last.t_end as usize;

            let mut total_matches: f64 = 0.0;
            let mut total_align_span: usize = 0;

            let mut last_q_end = first.q_start;
            let mut last_t_end = first.t_start;

            for (i, h) in bufs.chain.iter().enumerate() {
                let q_gap = (h.q_start - last_q_end).max(0);
                let t_gap = (h.t_start - last_t_end).max(0);
                if i > 0 {
                    // calculate gap sizes between previous MEM and current MEM
                    let aligned_gap = q_gap.min(t_gap);

                    // We assume that any non indel region has the background ANI
                    // But unaligned_gap = q_gap.abs_diff(t_gap) is multiplied by 0
                    // for matches (implicitly), since those can't ever match.
                    if aligned_gap > 2 {
                        let internal_gap = aligned_gap - 2;
                        total_matches += (internal_gap as f64) * regional_expected_ani;
                    }
                }

                total_align_span += q_gap.max(t_gap) as usize;

                let q_added_ani = h.q_end.saturating_sub(last_q_end.max(h.q_start));
                let t_added_ani = h.t_end.saturating_sub(last_t_end.max(h.t_start));

                total_matches += q_added_ani.min(t_added_ani) as f64;

                total_align_span += q_added_ani.max(t_added_ani) as usize;

                last_q_end = last_q_end.max(h.q_end);
                last_t_end = last_t_end.max(h.t_end);
            }

            let q_span = final_q_end - final_q_start;
            let t_span = final_t_end - final_t_start;
            let max_span = q_span.max(t_span);

            let actual_ani: f64 = if max_span > 0 {
                (total_matches / total_align_span as f64) * 100.0
            } else {
                0.0
            };

            if actual_ani >= config.min_ani && max_span >= config.min_len {
                bufs.hits.push(Hit {
                    //q_id: String::from_utf8_lossy(q_id).into_owned(), // Don't need to keep track
                    t_id: t_id,
                    q_size: q_size,
                    q_start: final_q_start,
                    q_end: final_q_end,
                    t_size: t_size,
                    t_start: final_t_start,
                    t_end: final_t_end,
                    strand,
                    ani: actual_ani,
                    score: bufs.scores[idx],
                });
            }
        }
    }
    anchors.clear();
    if is_fwd {
        bufs.anchors_fwd = anchors;
    } else {
        bufs.anchors_rev = anchors;
    }
}


pub fn process_query_sequence(
    y_index: &YIndex,
    bufs: &mut ThreadBuffers,
) {
    bufs.hits.clear();

    let q_seq_fwd = std::mem::take(&mut bufs.encoded_fwd);
    let q_seq_rev = std::mem::take(&mut bufs.encoded_rev);

    // Process forward strand
    get_anchors(&q_seq_fwd, y_index, bufs);
    align_strand(&q_seq_fwd, '+', y_index, bufs);
    let fwd_hits_count = bufs.hits.len();

    align_strand(&q_seq_rev, '-', y_index, bufs);

    let q_len: usize = q_seq_fwd.len();
    for hit in bufs.hits[fwd_hits_count..].iter_mut() {
        let orig_q_start: usize = hit.q_start;
        let orig_q_end: usize = hit.q_end;
        hit.q_start = q_len.saturating_sub(orig_q_end);
        hit.q_end = q_len.saturating_sub(orig_q_start);
    }

    bufs.encoded_fwd = q_seq_fwd;
    bufs.encoded_rev = q_seq_rev;
}

// Filters in place
pub fn filter_overlapping_hits(hits: &mut Vec<Hit>) {
    // Symmetric tiebreaker
    hits.sort_unstable_by(|a, b| {
        b.score.partial_cmp(&a.score)
        .unwrap_or(std::cmp::Ordering::Equal)
        .then_with(|| {
            let a_len = (a.q_end - a.q_start) + (a.t_end - a.t_start);
            let b_len = (b.q_end - b.q_start) + (b.t_end - b.t_start);
            b_len.cmp(&a_len)
        })
    });

    let mut filtered_count = 0;

    for i in 0..hits.len() {
        let hit = hits[i].clone();
        let mut overlaps = false;

        // Calculate lengths for both Query and Target
        let hit_q_len = hit.q_end - hit.q_start;
        let hit_t_len = hit.t_end - hit.t_start;

        for j in 0..filtered_count {
            let kept = &hits[j];

            // 1. Must map to the same target sequence to have a valid target overlap
            if hit.t_id == kept.t_id {
                let q_overlap_start = hit.q_start.max(kept.q_start);
                let q_overlap_end = hit.q_end.min(kept.q_end);

                let t_overlap_start = hit.t_start.max(kept.t_start);
                let t_overlap_end = hit.t_end.min(kept.t_end);

                // 2. Check if physical overlap exists on BOTH axes
                if q_overlap_start < q_overlap_end && t_overlap_start < t_overlap_end {
                    let q_overlap_len = q_overlap_end - q_overlap_start;
                    let t_overlap_len = t_overlap_end - t_overlap_start;

                    let q_ratio = q_overlap_len as f64 / hit_q_len as f64;
                    let t_ratio = t_overlap_len as f64 / hit_t_len as f64;

                    // 3. SYMMETRIZED CHECK: Overlap must be > 50% on BOTH Query and Target
                    if q_ratio > 0.5 && t_ratio > 0.5 {
                        overlaps = true;
                        break;
                    }
                }
            }
        }

        if !overlaps {
            hits[filtered_count] = hit;
            filtered_count += 1;
        }
    }

    hits.truncate(filtered_count);
}


// Unit tests initially generated by Gemini
#[cfg(test)]
mod encoding_tests {
    use super::*;

    #[test]
    fn test_parse_odd_kmer() {
        // Validates correct parsing of odd k-mers within the 9..=15 range
        assert_eq!(parse_odd_kmer("9"), Ok(9));
        assert_eq!(parse_odd_kmer("13"), Ok(13));
        assert_eq!(parse_odd_kmer("15"), Ok(15));

        // Ensures even numbers are rejected gracefully
        assert!(parse_odd_kmer("10").is_err());
        assert!(parse_odd_kmer("14").is_err());

        // Ensures out-of-bound odd numbers are rejected
        assert!(parse_odd_kmer("7").is_err());
        assert!(parse_odd_kmer("17").is_err());

        // Ensures completely invalid strings fail safely
        assert!(parse_odd_kmer("abc").is_err());
    }

    #[test]
    fn test_encode_sequence() {
        // Tests the core scalar/vectorized translation logic.
        // The bitwise & 0x0F lookup gracefully handles upper and lowercase:
        // A/a = 0, C/c = 1, G/g = 2, T/t/U/u = 3, N/Other = 128
        let input = b"ACGTUacgtN";
        let mut output = Vec::new();

        encode_sequence(input, &mut output);

        let expected = vec![
            0, 1, 2, 3, 3, // Upper ACGTU
            0, 1, 2, 3,    // Lower acgt
            128            // Invalid/N
        ];
        assert_eq!(output, expected);
    }

    #[test]
    fn test_reverse_complement() {
        // Tests the in-place reverse complement generation.
        // Original: A(0), C(1), T(3), G(2)
        let encoded_fwd = vec![0, 1, 3, 2];
        let mut encoded_rev = Vec::new();

        reverse_complement(&encoded_fwd, &mut encoded_rev);

        // Reversal yields GTCA -> [2, 3, 1, 0]
        // Bitwise XOR 3 complements -> CAGT -> [1, 0, 2, 3]
        let expected = vec![1, 0, 2, 3];
        assert_eq!(encoded_rev, expected);
    }

    #[test]
    fn test_rev_comp_u32() {
        // Test packing and reverse-complementing a k-mer (K=3)
        // Forward: AAC -> A=0, C=1. Encoded: [0, 0, 1]
        // Packed u32 for K=3: (0 << 4) | (0 << 2) | 1 = 1
        let h_fwd = 1;
        let k = 3;

        // Reverse complement of AAC is GTT
        // G=2, T=3. Encoded: [2, 3, 3]
        // Packed: (2 << 4) | (3 << 2) | 3 = 32 + 12 + 3 = 47
        let h_rev = rev_comp_u32(h_fwd, k);

        assert_eq!(h_rev, 47);
    }

    #[test]
    fn test_is_invalid_k8_16() {
        // Tests the 64-bit overlapping bitwise chunk evaluation
        // Valid bases are 0,1,2,3 (meaning bits 2-7 are always 0)
        let valid_chunk: [u8; 9] = [0, 1, 2, 3, 0, 1, 2, 3, 0];

        // This should pass because no bytes trigger the 0xFCFC... mask
        assert!(!is_invalid_k8_16::<9>(&valid_chunk));

        // Invalid chunk containing an 'N' (128)
        let invalid_chunk: [u8; 9] = [0, 1, 2, 128, 0, 1, 2, 3, 0];

        // 128 (0x80) hits the 0xFC mask, triggering the invalid flag
        assert!(is_invalid_k8_16::<9>(&invalid_chunk));
    }

    #[test]
    fn test_hash_k_fast() {
        // Tests the compile-time unrolled tree-split hasher
        // Sequence of K=9: 8 A's (0) ending in C (1)
        let chunk_a: [u8; 9] = [0, 0, 0, 0, 0, 0, 0, 0, 1];

        // Should simply be 1 at the end of the shift
        let hash_a = hash_k_fast::<9>(&chunk_a);
        assert_eq!(hash_a, 1);

        // Sequence of K=9: Padded with 5 A's, ending in ACGT
        // ACGT -> [0, 1, 2, 3]
        let chunk_b: [u8; 9] = [0, 0, 0, 0, 0, 0, 1, 2, 3];

        // Expected bit math: (1 << 4) | (2 << 2) | 3 = 16 + 8 + 3 = 27
        let hash_b = hash_k_fast::<9>(&chunk_b);
        assert_eq!(hash_b, 27);
    }
}

#[cfg(test)]
mod pipeline_tests {
    use super::*;

    // Safely initialize global CONFIG once across all parallel test threads
    fn init_test_config() {
        let _ = CONFIG.set(AlignConfig {
            kmer_size: 9,
            stride: Some(9),
                           min_ani: 85.0,
                           min_len: 9,
                           max_seed_multiplicity: 1000,
                           max_chain_gap: 100,
        });
    }

    // Helper to build a minimal mock YIndex directly in memory
    fn create_mock_index(target_encoded: Vec<u8>, kmer: u32, is_rev: bool) -> YIndex {
        let mut kmer_map = FxHashMap::default();
        kmer_map.insert(kmer, (0, 1));

        let pos_val = if is_rev { 0x8000_0000 } else { 0 };

        YIndex {
            kmer_map,
            flat_kmers: vec![TargetPos { seq_id: 0, pos: pos_val }],
            presence_bits: vec![u64::MAX; 65536], // Pass all Bloom filter checks
            seq_data: vec![target_encoded],
            seq_names: vec![b"target1".to_vec()],
        }
    }

    #[test]
    fn test_filter_overlapping_hits() {
        // Tests filtering of hits on the same target with > 50% overlap on both query and target.
        let mut hits = vec![
            // Hit 0: Higher score (100), spans Q: 0..50, T: 0..50
            Hit {
                t_id: 0,
                q_size: 100,
                q_start: 0,
                q_end: 50,
                t_size: 100,
                t_start: 0,
                t_end: 50,
                strand: '+',
                ani: 95.0,
                score: 100,
            },
            // Hit 1: Lower score (50), heavily overlaps Hit 0 on both Q (10..40) and T (10..40) -> Should be filtered
            Hit {
                t_id: 0,
                q_size: 100,
                q_start: 10,
                q_end: 40,
                t_size: 100,
                t_start: 10,
                t_end: 40,
                strand: '+',
                ani: 90.0,
                score: 50,
            },
            // Hit 2: Lower score (40), but maps to a DIFFERENT target (t_id = 1) -> Should be retained
            Hit {
                t_id: 1,
                q_size: 100,
                q_start: 0,
                q_end: 50,
                t_size: 100,
                t_start: 0,
                t_end: 50,
                strand: '+',
                ani: 92.0,
                score: 40,
            },
        ];

        filter_overlapping_hits(&mut hits);

        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].score, 100);
        assert_eq!(hits[1].t_id, 1);
    }

    #[test]
    fn test_get_anchors_forward_match() {
        init_test_config();
        let mut bufs = ThreadBuffers::new();

        // 9 bases of A's: packed u32 hash = 0
        let q_seq = vec![0, 0, 0, 0, 0, 0, 0, 0, 0];
        let y_index = create_mock_index(q_seq.clone(), 0, false);

        get_anchors(&q_seq, &y_index, &mut bufs);

        // Standard forward match should produce a forward anchor on t_id 0
        assert_eq!(bufs.anchors_fwd.len(), 1);
        assert_eq!(bufs.anchors_rev.len(), 0);
        assert_eq!(bufs.anchors_fwd[0].t_id, 0);
        assert_eq!(bufs.anchors_fwd[0].q_pos, 0);
        assert_eq!(bufs.anchors_fwd[0].t_pos, 0);
    }

    #[test]
    fn test_align_strand_mem_extension() {
        init_test_config();
        let mut bufs = ThreadBuffers::new();

        // Query & Target: 12 bases matching completely (AAAAAAAAAAAA -> [0; 12])
        let q_seq = vec![0; 12];
        let t_seq = vec![0; 12];
        let y_index = create_mock_index(t_seq, 0, false);

        // Manually place a 9-mer seed anchor at Q:0, T:0
        bufs.anchors_fwd.push(Anchor {
            t_id: 0,
            q_pos: 0,
            t_pos: 0,
            diag: 0,
        });

        align_strand(&q_seq, '+', &y_index, &mut bufs);

        // Alignment should extend 9-mer seed into full 12-base match hit
        assert_eq!(bufs.hits.len(), 1);
        let hit = &bufs.hits[0];
        assert_eq!(hit.q_start, 0);
        assert_eq!(hit.q_end, 12);
        assert_eq!(hit.t_start, 0);
        assert_eq!(hit.t_end, 12);
        assert!((hit.ani - 100.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_process_query_sequence_reverse_strand_coordinates() {
        init_test_config();
        let mut bufs = ThreadBuffers::new();

        // Query (Forward): 10 bases [A, A, A, A, A, A, A, A, A, C]
        bufs.encoded_fwd = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
        reverse_complement(&bufs.encoded_fwd, &mut bufs.encoded_rev);

        // Target (9 T's): 2-bit hash is 0x3FFFF, but canonical hash min(fwd, rev) is 0.
        // Target is reverse relative to canonical 0, so is_rev = true (bit 31 set).
        let t_seq = vec![3; 9];
        let canonical_hash = 0;
        let y_index = create_mock_index(t_seq, canonical_hash, true);

        process_query_sequence(&y_index, &mut bufs);

        // Verify reverse strand coordinate conversion back to forward query space
        assert_eq!(bufs.hits.len(), 1);
        let hit = &bufs.hits[0];
        assert_eq!(hit.strand, '-');

        // Reverse match spans rev_q: 1..10. Re-mapped to forward Q space: (10 - 10) .. (10 - 1) => 0..9
        assert_eq!(hit.q_start, 0);
        assert_eq!(hit.q_end, 9);
    }
}
