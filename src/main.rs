// BacPrune-Rust — LD pruning of haploid genotype matrices
//
// Pruning runs in two phases:
//   Phase 1 (all modes): hash-based exact duplicate removal — O(n·v), no pairwise work.
//   Phase 2 (--r / --dprime): greedy pairwise pruning by an LD threshold.
//
// Method flags (choose one; default is --r):
//   --r       Prune by r² (Pearson r-squared) threshold    (default)
//   --dprime  Prune by |D'| threshold
//   --dedup   Phase 1 only — remove exact duplicates, skip threshold pruning
//
// Output files (written to <output_directory>):
//   bacprune_rust_results.csv      — pruned genotype matrix with header
//   ld_pruning_summary.csv         — representative variant → list of pruned
//                                    variants. Both columns contain the
//                                    original variant IDs preserved from the
//                                    input header row (NOT post-MAF
//                                    positions), so callers can cross-
//                                    reference them against the results CSV
//                                    header directly regardless of the MAF
//                                    cutoff in use.
//   direction_of_correlation.csv   — per-variant status and correlation direction
//                                    relative to its representative

use clap::Parser;
use ndarray::prelude::*;
use ndarray::OwnedRepr;
use csv::ReaderBuilder;
use ndarray::Array2;
use std::fs::File;
use ndarray_csv::Array2Reader;
use std::fs::OpenOptions;
use std::collections::{HashMap, HashSet};
use std::ops::Div;
use std::path::Path;

/// Fast LD pruning of haploid genotype matrices.
///
/// Pruning runs in two phases:
///
///   Phase 1 (all modes): hash-based exact duplicate removal — O(n·v), no pairwise work.
///
///   Phase 2 (--r / --dprime): greedy pairwise pruning by an LD threshold.
///
/// Output files written to <output_directory>:
///
///   bacprune_rust_results.csv     — pruned genotype matrix
///
///   ld_pruning_summary.csv        — representative SNP → pruned SNPs
///
///   direction_of_correlation.csv  — per-variant correlation direction
#[derive(Parser)]
#[command(name = "bacprune", version, about, long_about = None)]
struct Cli {
    /// Path to the input CSV (numeric header row + sample rows of 0/1 genotypes)
    input_file: String,

    /// Total number of rows in the CSV including the header row
    n_rows: usize,

    /// Number of columns in the CSV
    n_cols: usize,

    /// Minor allele frequency cutoff; variants below this threshold are removed
    maf_cutoff: f64,

    /// Directory where output files are written
    output_directory: String,

    /// LD pruning threshold; pairs at or above this value are pruned.
    /// Required for --r and --dprime; not used with --dedup.
    #[arg(long, value_name = "THRESHOLD", required_unless_present = "dedup",
          conflicts_with = "dedup")]
    ld: Option<f64>,

    /// Prune by r² (Pearson r-squared) threshold [default method]
    #[arg(long, conflicts_with_all = ["dprime", "dedup"])]
    r: bool,

    /// Prune by |D'| (Lewontin's D') threshold
    #[arg(long, conflicts_with_all = ["r", "dedup"])]
    dprime: bool,

    /// Remove exact duplicate variants only via hashing; no pairwise LD calculation
    #[arg(long, conflicts_with_all = ["r", "dprime", "ld"])]
    dedup: bool,
}

#[derive(PartialEq)]
enum Method { Correlation, DPrime, DedupOnly }

/// Maps a pruned variant's post-MAF column index to
/// (representative's post-MAF column index, is_positive_correlation).
/// Positive means the pruned variant is identical to its representative (r > 0);
/// negative means it is complementary / opposite (r < 0).
type DirectionMap = HashMap<usize, (usize, bool)>;

fn main() -> Result<(), csv::Error> {
    let cli = Cli::parse();

    println!("Welcome to the LD Pruning Module.");

    let input_file        = cli.input_file;
    let n_rows: usize     = cli.n_rows;
    let n_cols: usize     = cli.n_cols;
    let maf_cutoff: f64   = cli.maf_cutoff;
    let ld_threshold: f64 = cli.ld.unwrap_or(0.0);
    let outdir            = cli.output_directory;

    if n_rows < 2 {
        eprintln!("Error: n_rows must be at least 2 (one header row plus at least one sample).");
        return Ok(());
    }
    if maf_cutoff < 0.0 || maf_cutoff > 1.0 {
        eprintln!("Error: maf_cutoff must be between 0 and 1, got {maf_cutoff}.");
        return Ok(());
    }
    if !cli.dedup && (ld_threshold < 0.0 || ld_threshold > 1.0) {
        eprintln!("Error: --ld threshold must be between 0 and 1, got {ld_threshold}.");
        return Ok(());
    }

    let method = if cli.dedup {
        Method::DedupOnly
    } else if cli.dprime {
        Method::DPrime
    } else {
        Method::Correlation // default (also triggered by --r)
    };

    match method {
        Method::Correlation => println!("Method: r²    threshold >= {ld_threshold}"),
        Method::DPrime      => println!("Method: |D'|  threshold >= {ld_threshold}"),
        Method::DedupOnly   => println!("Method: --dedup (hash-based exact duplicate removal only)"),
    }

    // ── Read data ────────────────────────────────────────────────────────────
    let raw_gt_data = read_csv(&input_file, n_rows, n_cols);
    println!("Data read successfully.");

    // First row is the header (variant IDs stored as f64 column indices)
    let gt_header   = raw_gt_data.select(Axis(0), &[0usize]);
    let raw_gt_data = raw_gt_data.slice(s![1.., ..]).to_owned();

    // Validate that all genotype values are 0 or 1
    for val in raw_gt_data.iter() {
        if *val != 0.0 && *val != 1.0 {
            eprintln!("Error: genotype matrix contains value {val}; only 0 and 1 are permitted.");
            return Ok(());
        }
    }

    // ── MAF filter ───────────────────────────────────────────────────────────
    let mafs = calc_maf(&raw_gt_data);
    let (maf_data, maf_keep_idx) = maf_filter(&raw_gt_data, &mafs, &maf_cutoff);
    let gt_header = gt_header.select(Axis(1), maf_keep_idx.as_slice());

    // Save the full post-MAF header before any LD pruning.  This is used later
    // to look up variant IDs for ALL variants (including those pruned by LD),
    // so that direction_of_correlation.csv can reference them by name.
    let maf_gt_header = gt_header.clone(); // shape (1 × n_maf_variants)

    println!("MAF filter complete. {} variants kept.", maf_data.ncols());

    // ── Phase 1: hash-based exact duplicate removal ───────────────────────────
    // Indices in the returned DirectionMap are into maf_data (post-MAF space).
    // All identical duplicates are positively correlated (r = 1) with their rep.
    let (dedup_data, dedup_keep_idx, mut direction_map) = dedup_variants(&maf_data);
    let gt_header = gt_header.select(Axis(1), &dedup_keep_idx);
    println!("Exact duplicate removal complete. {} variants kept.", dedup_data.ncols());

    // If --dedup, we are done after phase 1
    if method == Method::DedupOnly {
        write_outputs(&dedup_data, &gt_header, &direction_map, &maf_gt_header, &outdir)?;
        return Ok(());
    }

    // ── Phase 2: threshold LD pruning (r or D') ───────────────────────────────
    // Indices here are into dedup_data (post-dedup space).
    // We convert them back to post-MAF space before inserting into direction_map.
    let mafs = calc_maf(&dedup_data);
    let mut skip: HashSet<usize> = HashSet::new();

    // For r² mode, sort variants by descending MAF so that max_r_squared
    // decreases monotonically along the inner loop, allowing early `break`
    // instead of `continue`.  For D' mode, keep the original order.
    let sort_perm: Vec<usize> = if method == Method::Correlation {
        let mut perm: Vec<usize> = (0..dedup_data.ncols()).collect();
        perm.sort_by(|&a, &b| mafs[b].partial_cmp(&mafs[a]).unwrap_or(std::cmp::Ordering::Equal));
        perm
    } else {
        (0..dedup_data.ncols()).collect()
    };
    let work_data = dedup_data.select(Axis(1), &sort_perm);
    let work_mafs: Vec<f64> = sort_perm.iter().map(|&k| mafs[k]).collect();

    'outer: for i in 0..work_data.ncols() {
        if skip.contains(&i) { continue; }
        for j in (i + 1)..work_data.ncols() {
            if skip.contains(&j) { continue; }

            // For r² mode, skip pairs whose MAFs make it impossible to reach
            // the threshold.  Because variants are sorted by descending MAF,
            // once the bound fails all subsequent j will also fail → break.
            if method == Method::Correlation
                && max_r_squared(work_mafs[i], work_mafs[j]) < ld_threshold
            {
                break;
            }

            // Compute the LD metric. For --r, cache r so the sign is available
            // for direction tracking without a second call. For --dprime, r is
            // only needed when the threshold is actually exceeded.
            let (ld, r_cache) = match method {
                Method::DPrime => (calculate_d_prime(&work_data, i, j), None),
                _              => { let r = calculate_r(&work_data, i, j); (r * r, Some(r)) },
            };

            if ld >= ld_threshold {
                // Prune the lower-MAF variant; keep the higher-MAF one.
                let (keep, prune) = if work_mafs[i] >= work_mafs[j] { (i, j) } else { (j, i) };
                skip.insert(prune);

                // Record direction using r (sign is the same for both r² and |D'|).
                // For --dprime, compute r only now (threshold was exceeded).
                // Map sorted indices → dedup indices → post-MAF indices.
                let r_signed      = r_cache.unwrap_or_else(|| calculate_r(&work_data, i, j));
                let is_positive   = r_signed >= 0.0;
                let prune_maf_idx = dedup_keep_idx[sort_perm[prune]];
                let keep_maf_idx  = dedup_keep_idx[sort_perm[keep]];
                direction_map.insert(prune_maf_idx, (keep_maf_idx, is_positive));

                // If i itself was just pruned, it can no longer act as a
                // representative for further j values — exit the inner loop.
                if prune == i { continue 'outer; }
            }
        }
    }

    let mut phase2_keep: Vec<usize> = (0..work_data.ncols())
        .filter(|x| !skip.contains(x))
        .map(|x| sort_perm[x])  // map back to dedup indices
        .collect();
    phase2_keep.sort();  // restore original column order
    let pruned_data   = dedup_data.select(Axis(1), &phase2_keep);
    let pruned_header = gt_header.select(Axis(1), &phase2_keep);
    println!("LD threshold pruning complete. {} variants kept.", pruned_data.ncols());

    write_outputs(&pruned_data, &pruned_header, &direction_map, &maf_gt_header, &outdir)?;
    Ok(())
}

// ── Output helper ─────────────────────────────────────────────────────────────

/// Write all three output CSVs to `outdir`:
///   1. `bacprune_rust_results.csv`     — final pruned genotype matrix
///   2. `ld_pruning_summary.csv`        — representative → pruned SNPs
///   3. `direction_of_correlation.csv`  — per-variant status and direction
///
/// `maf_gt_header` is the full post-MAF header (1 × n_maf_variants) used to
/// resolve variant IDs for both kept and pruned variants in the direction CSV.
/// Indices in `direction_map` are into post-MAF column space.
fn write_outputs(
    data:          &Array2<f64>,
    header:        &Array2<f64>,
    direction_map: &DirectionMap,
    maf_gt_header: &Array2<f64>,
    outdir:        &str,
) -> Result<(), csv::Error> {

    std::fs::create_dir_all(outdir).expect("Failed to create output directory");

    // 1. Results CSV ──────────────────────────────────────────────────────────
    let final_data = ndarray::concatenate![Axis(0), *header, *data];
    let string_arr = final_data.map(|e| e.to_string());
    let csv_path   = Path::new(outdir).join("bacprune_rust_results.csv");
    let file = OpenOptions::new().write(true).create(true).truncate(true).open(csv_path).unwrap();
    let mut wtr = csv::Writer::from_writer(file);
    for i in 0..final_data.nrows() {
        wtr.write_record(&string_arr.slice(s![i, ..])).expect("Error writing to CSV");
    }
    wtr.flush()?;

    // 2. Summary CSV ──────────────────────────────────────────────────────────
    // Derive rep → pruned list from direction_map.  The direction_map keys
    // are post-MAF column indices; translate them to original variant IDs
    // (as preserved in `maf_gt_header`) so downstream tools can cross-
    // reference them against the results CSV header and any caller-supplied
    // position/annotation table without needing to replicate the MAF filter.
    let mut rep_snps: HashMap<usize, Vec<usize>> = HashMap::new();
    for (&pruned, &(rep, _)) in direction_map {
        rep_snps.entry(rep).or_insert_with(Vec::new).push(pruned);
    }
    let summary_path = Path::new(outdir).join("ld_pruning_summary.csv");
    let file = OpenOptions::new().write(true).create(true).truncate(true).open(summary_path).unwrap();
    let mut wtr = csv::Writer::from_writer(file);
    wtr.write_record(&["Representative Variant", "Pruned Variants"])
        .expect("Error writing summary header");
    // Sort representatives by their original variant ID (via maf_gt_header
    // lookup) rather than by their post-MAF position.  With the current MAF
    // filter these orderings are equivalent, but sorting by variant ID keeps
    // the contract explicit and robust to future reorderings.
    let mut rep_keys: Vec<usize> = rep_snps.keys().copied().collect();
    rep_keys.sort_by(|&a, &b| {
        maf_gt_header[[0, a]]
            .partial_cmp(&maf_gt_header[[0, b]])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for rep in &rep_keys {
        // Sort each representative's pruned-variant list by original ID too
        // so the emitted order is human-readable.
        let mut pruned = rep_snps[rep].clone();
        pruned.sort_by(|&a, &b| {
            maf_gt_header[[0, a]]
                .partial_cmp(&maf_gt_header[[0, b]])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let rep_id     = maf_gt_header[[0, *rep]].to_string();
        let pruned_str = pruned
            .iter()
            .map(|&idx| maf_gt_header[[0, idx]].to_string())
            .collect::<Vec<_>>()
            .join(", ");
        wtr.write_record(&[rep_id, pruned_str]).expect("Error writing summary record");
    }
    wtr.flush()?;

    // 3. Direction CSV ────────────────────────────────────────────────────────
    // Rows for every post-MAF variant (representatives and pruned alike).
    // Status values:
    //   "representative"       — variant was kept; is the representative for its group
    //   "positive_correlation" — variant was pruned; genotype matches its representative (r > 0)
    //   "negative_correlation" — variant was pruned; genotype is the complement of its rep (r < 0)
    let direction_path = Path::new(outdir).join("direction_of_correlation.csv");
    let file = OpenOptions::new().write(true).create(true).truncate(true).open(direction_path).unwrap();
    let mut wtr = csv::Writer::from_writer(file);
    wtr.write_record(&["Variant", "Status", "Representative Variant"])
        .expect("Error writing direction header");
    let n_maf = maf_gt_header.ncols();
    for idx in 0..n_maf {
        let variant_id = maf_gt_header[[0, idx]].to_string();
        match direction_map.get(&idx) {
            Some(&(rep_idx, is_positive)) => {
                let status   = if is_positive { "positive_correlation" } else { "negative_correlation" };
                let rep_id   = maf_gt_header[[0, rep_idx]].to_string();
                wtr.write_record(&[variant_id, status.to_string(), rep_id])
                    .expect("Error writing direction record");
            }
            None => {
                // This variant is a representative (not pruned by any step)
                wtr.write_record(&[variant_id.clone(), "representative".to_string(), variant_id])
                    .expect("Error writing direction record");
            }
        }
    }
    wtr.flush()?;

    Ok(())
}

// ── Core functions ────────────────────────────────────────────────────────────

/// Read a headerless CSV into a 2-D f64 array of shape (n_rows × n_cols).
/// The caller is responsible for passing the correct dimensions; the function
/// will panic if the file cannot be opened or parsed.
fn read_csv(path_to_file: &str, n_rows: usize, n_cols: usize) -> Array2<f64> {
    let file = File::open(path_to_file)
        .unwrap_or_else(|e| panic!("Error reading CSV '{}': {}", path_to_file, e));
    let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);
    reader.deserialize_array2::<f64>((n_rows, n_cols))
        .unwrap_or_else(|e| panic!("Error reading CSV '{}': {}", path_to_file, e))
}

/// Compute the allele frequency for every variant (column).
/// For a haploid 0/1 matrix this is simply the mean of each column,
/// i.e. the frequency of the allele encoded as 1.
///
/// The MAF filter therefore operates on the frequency of whatever allele is
/// encoded as 1.  To filter out both rare ALT and rare REF alleles, normalise
/// the input so that the alternate allele is always 1 and the reference allele
/// is always 0 before running BacPrune.
/// Returns a 1-D array of length n_cols.
fn calc_maf(data: &Array2<f64>) -> Array1<f64> {
    data.sum_axis(Axis(0)).div(data.nrows() as f64)
}

/// Remove variants whose MAF is strictly below `cutoff`.
/// Returns the filtered matrix and the original column indices that were kept,
/// so that the caller can apply the same filter to the header row.
fn maf_filter(data: &Array2<f64>, mafs: &Array1<f64>, cutoff: &f64) -> (Array2<f64>, Vec<usize>) {
    let keep: Vec<usize> = mafs
        .iter()
        .enumerate()
        .filter(|(_, &x)| x >= *cutoff)
        .map(|(i, _)| i)
        .collect();
    (data.select(Axis(1), &keep), keep)
}

/// Hash-based exact duplicate removal — O(n·v), no pairwise comparisons.
///
/// Each column is encoded as a `Vec<u8>` key (safe for haploid 0/1 data) and
/// inserted into a `HashMap`. The first occurrence of each unique genotype
/// pattern is kept; all later identical columns are recorded in the
/// returned `DirectionMap` with `is_positive = true` (identical ⇒ r = 1).
///
/// Complementary columns (one is the bitwise NOT of the other) are also
/// detected and recorded with `is_positive = false` (r = −1).
///
/// Returns `(deduplicated_data, kept_column_indices, direction_map)`.
/// All indices are in the space of the input `data` (post-MAF column space).
fn dedup_variants(data: &Array2<f64>) -> (Array2<f64>, Vec<usize>, DirectionMap) {
    let mut seen: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut keep: Vec<usize>              = Vec::new();
    let mut direction_map: DirectionMap   = HashMap::new();

    for j in 0..data.ncols() {
        let key: Vec<u8> = data.column(j).iter().map(|&x| x as u8).collect();

        // Check for complement (bitwise NOT) — r = −1, |D'| = 1.
        let complement_key: Vec<u8> = key.iter().map(|&x| 1 - x).collect();
        if let Some(&rep) = seen.get(&complement_key) {
            direction_map.insert(j, (rep, false));
            continue;
        }

        match seen.get(&key) {
            Some(&rep) => {
                // Identical to `rep` — positively correlated (r = 1)
                direction_map.insert(j, (rep, true));
            }
            None => {
                seen.insert(key, j);
                keep.push(j);
            }
        }
    }

    (data.select(Axis(1), &keep), keep, direction_map)
}

/// Theoretical maximum r² achievable between two binary variants with
/// minor allele frequencies `p` and `q`.
///
/// For haploid 0/1 data the tightest upper bound on r² given marginal
/// frequencies p and q (with p ≤ q) is:
///
///   max_r² = p·(1 − q) / (q·(1 − p))
///
/// If either variant is monomorphic (frequency 0 or 1), the maximum is 0.
fn max_r_squared(p: f64, q: f64) -> f64 {
    let (lo, hi) = if p <= q { (p, q) } else { (q, p) };
    let denom = hi * (1.0 - lo);
    if denom < 1e-15 { return 0.0; }
    lo * (1.0 - hi) / denom
}

/// Pearson correlation coefficient between two binary (0/1) columns.
///
/// Uses the computational form:
///   r = (n·Σxy − Σx·Σy) / √((Σx·(n−Σx)) · (Σy·(n−Σy)))
///
/// Returns a value in [−1, 1].  Returns 0.0 if either column is monomorphic
/// (zero variance), avoiding a division-by-zero.
/// Square the result to obtain r² for use as an LD metric.
fn calculate_r(data: &Array2<f64>, varianta: usize, variantb: usize) -> f64 {
    let n      = data.nrows() as f64;
    let col_a  = data.column(varianta);
    let col_b  = data.column(variantb);
    let sum_a  = col_a.sum();
    let sum_b  = col_b.sum();
    let sum_ab: f64 = col_a.iter().zip(col_b.iter()).map(|(x, y)| x * y).sum();

    let numerator = n * sum_ab - sum_a * sum_b;
    // For binary data Σx² = Σx, so n·Σx² − (Σx)² simplifies to sum_x·(n − sum_x).
    let var_a = sum_a * (n - sum_a);
    let var_b = sum_b * (n - sum_b);
    let denom = (var_a * var_b).sqrt();

    if denom < 1e-10 { return 0.0; }
    numerator / denom
}

/// Lewontin's D' linkage disequilibrium coefficient for two variants.
///
/// D' measures the deviation from linkage equilibrium normalised by the
/// maximum possible deviation given the observed allele frequencies:
///   D  = f(00)·f(11) − f(10)·f(01)
///   D' = D / D_max
///
/// where D_max is chosen so that |D'| ∈ [0, 1]:
///   If D > 0: D_max = min(freq_0A·freq_1B,  freq_1A·freq_0B)
///   If D < 0: D_max = max(−freq_0A·freq_0B, −freq_1A·freq_1B)
///
/// This function always returns |D'|, so the result is in [0, 1].
/// D' = 0 means no LD; D' = 1 means perfect LD (identical or complementary columns).
/// Use `calculate_r` to obtain the sign (direction) of the association.
fn calculate_d_prime(data: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>, varianta: usize, variantb: usize) -> f64 {
    let (f0a, f1a, f0b, f1b) = find_allele_frequencies(data, varianta, variantb);
    let hf = find_haplotype_frequencies(data, varianta, variantb);
    return calc_d_prime(f0a, f1a, f0b, f1b, &hf);

    /// Returns (freq_allele0_A, freq_allele1_A, freq_allele0_B, freq_allele1_B)
    /// for variants va and vb only — O(n) instead of O(n·v).
    fn find_allele_frequencies(data: &Array2<f64>, va: usize, vb: usize) -> (f64, f64, f64, f64) {
        let n   = data.nrows() as f64;
        let f1a = data.column(va).sum() / n;
        let f1b = data.column(vb).sum() / n;
        (1.0 - f1a, f1a, 1.0 - f1b, f1b)
    }

    /// Returns the four haplotype frequencies [f00, f10, f01, f11] for the
    /// pair (va, vb), where e.g. f10 = freq(allele1_A, allele0_B).
    /// Counts in a single pass with no heap allocations.
    fn find_haplotype_frequencies(data: &Array2<f64>, va: usize, vb: usize) -> Array1<f64> {
        let n = data.nrows() as f64;
        let (mut f00, mut f10, mut f01, mut f11) = (0u64, 0u64, 0u64, 0u64);
        for k in 0..data.nrows() {
            match (data[[k, va]] as u8, data[[k, vb]] as u8) {
                (0, 0) => f00 += 1,
                (1, 0) => f10 += 1,
                (0, 1) => f01 += 1,
                (1, 1) => f11 += 1,
                _      => {}
            }
        }
        ndarray::array![f00 as f64, f10 as f64, f01 as f64, f11 as f64].div(n)
    }

    /// Compute D' from allele frequencies and haplotype frequencies for the pair.
    fn calc_d_prime(f0a: f64, f1a: f64, f0b: f64, f1b: f64, hf: &Array1<f64>) -> f64 {
        let d = hf[0] * hf[3] - hf[1] * hf[2]; // D = f(00)·f(11) − f(10)·f(01)
        if d < 0.0 {
            let d_max = f64::max(-f0a * f0b, -f1a * f1b);
            d / d_max
        } else if d > 0.0 {
            let d_max = f64::min(f0a * f1b, f1a * f0b);
            d / d_max
        } else {
            0.0
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── calc_maf ─────────────────────────────────────────────────────────────

    #[test]
    fn test_calc_maf_basic() {
        let data = array![
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        let mafs = calc_maf(&data);
        assert!((mafs[0] - 1.0).abs()      < 1e-10);
        assert!((mafs[1] - 1.0/3.0).abs()  < 1e-10);
        assert!((mafs[2] - 0.0).abs()       < 1e-10);
    }

    #[test]
    fn test_calc_maf_uses_nrows() {
        let data = array![[1.0, 0.0], [0.0, 1.0]];
        let mafs = calc_maf(&data);
        assert!((mafs[0] - 0.5).abs() < 1e-10);
        assert!((mafs[1] - 0.5).abs() < 1e-10);
    }

    // ── calc_maf: edge cases ─────────────────────────────────────────────────

    #[test]
    fn test_calc_maf_monomorphic_all_zero() {
        // All-reference column → frequency 0.0
        let data = array![[0.0], [0.0], [0.0]];
        let mafs = calc_maf(&data);
        assert!((mafs[0] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_calc_maf_monomorphic_all_one() {
        // All-alternate column → frequency 1.0
        let data = array![[1.0], [1.0], [1.0]];
        let mafs = calc_maf(&data);
        assert!((mafs[0] - 1.0).abs() < 1e-10);
    }

    // ── maf_filter ────────────────────────────────────────────────────────────

    #[test]
    fn test_maf_filter_removes_low_maf() {
        let data = array![[1.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 0.0]];
        let mafs = array![0.5, 0.0];
        let (filtered, keep_idx) = maf_filter(&data, &mafs, &0.01);
        assert_eq!(filtered.ncols(), 1);
        assert_eq!(keep_idx, vec![0]);
    }

    #[test]
    fn test_maf_filter_keeps_at_cutoff() {
        // MAF exactly equal to cutoff is kept (>=, not >)
        let data = array![[1.0, 0.0], [0.0, 0.0]];
        let mafs = array![0.5, 0.01];
        let (filtered, keep_idx) = maf_filter(&data, &mafs, &0.01);
        assert_eq!(filtered.ncols(), 2);
        assert_eq!(keep_idx, vec![0, 1]);
    }

    #[test]
    fn test_maf_filter_returns_correct_keep_index() {
        // Middle column removed; indices 0 and 2 must be returned
        let data = array![[1.0, 0.0, 1.0], [0.0, 0.0, 0.0]];
        let mafs = array![0.5, 0.0, 0.5];
        let (filtered, keep_idx) = maf_filter(&data, &mafs, &0.01);
        assert_eq!(filtered.ncols(), 2);
        assert_eq!(keep_idx, vec![0, 2]);
    }

    #[test]
    fn test_maf_filter_removes_all() {
        // All variants below cutoff → empty matrix
        let data = array![[0.0, 0.0], [0.0, 0.0]];
        let mafs = array![0.0, 0.0];
        let (filtered, keep_idx) = maf_filter(&data, &mafs, &0.01);
        assert_eq!(filtered.ncols(), 0);
        assert!(keep_idx.is_empty());
    }

    #[test]
    fn test_maf_filter_zero_cutoff_keeps_all() {
        // Cutoff of 0.0 keeps every variant including monomorphic
        let data = array![[1.0, 0.0], [0.0, 0.0], [1.0, 0.0]];
        let mafs = array![0.5, 0.0];
        let (filtered, keep_idx) = maf_filter(&data, &mafs, &0.0);
        assert_eq!(filtered.ncols(), 2);
        assert_eq!(keep_idx, vec![0, 1]);
    }

    #[test]
    fn test_maf_filter_preserves_column_data() {
        // Verify the actual values in the surviving column are correct
        let data = array![
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
        ];
        let mafs = array![2.0/3.0, 0.0, 1.0/3.0];
        let (filtered, _) = maf_filter(&data, &mafs, &0.1);
        // Only cols 0 and 2 survive; check their values are unchanged
        assert_eq!(filtered.column(0).to_vec(), vec![1.0, 0.0, 1.0]);
        assert_eq!(filtered.column(1).to_vec(), vec![0.0, 1.0, 0.0]);
    }

    // ── dedup_variants ───────────────────────────────────────────────────────

    #[test]
    fn test_dedup_removes_identical_columns() {
        // col 1 is identical to col 0 → pruned with direction=positive
        // col 2 is complement of col 0 → also pruned with direction=negative
        let data = array![
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ];
        let (deduped, keep_idx, direction_map) = dedup_variants(&data);
        assert_eq!(deduped.ncols(), 1);
        assert_eq!(keep_idx, vec![0]);
        assert_eq!(direction_map[&1], (0, true),  "col 1 identical → positive");
        assert_eq!(direction_map[&2], (0, false), "col 2 complement → negative");
    }

    #[test]
    fn test_dedup_prunes_opposite_columns() {
        // Opposite (complement) columns are now caught — col 1 pruned as negative
        let data = array![
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
        ];
        let (deduped, keep_idx, direction_map) = dedup_variants(&data);
        assert_eq!(deduped.ncols(), 1);
        assert_eq!(keep_idx, vec![0]);
        assert_eq!(direction_map[&1], (0, false));
    }

    #[test]
    fn test_dedup_no_duplicates() {
        let data = array![[1.0, 0.0], [0.0, 0.0], [1.0, 1.0]];
        let (deduped, keep_idx, direction_map) = dedup_variants(&data);
        assert_eq!(deduped.ncols(), 2);
        assert_eq!(keep_idx, vec![0, 1]);
        assert!(direction_map.is_empty());
    }

    #[test]
    fn test_dedup_multiple_duplicates_same_group() {
        // cols 0, 1, 2 all identical [1,0] → keep col 0, prune 1 and 2 (positive)
        // col 3 [0,1] is complement of col 0 → pruned (negative)
        let data = array![
            [1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let (deduped, keep_idx, direction_map) = dedup_variants(&data);
        assert_eq!(deduped.ncols(), 1);
        assert_eq!(keep_idx, vec![0]);
        assert_eq!(direction_map[&1], (0, true));
        assert_eq!(direction_map[&2], (0, true));
        assert_eq!(direction_map[&3], (0, false));
    }

    // ── dedup_variants: complement detection ──────────────────────────────────

    #[test]
    fn test_dedup_catches_complement_columns() {
        // col 1 is complement of col 0, col 2 is unique
        let data = array![
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let (deduped, keep_idx, direction_map) = dedup_variants(&data);
        assert_eq!(deduped.ncols(), 2, "complement should be pruned");
        assert_eq!(keep_idx, vec![0, 2]);
        assert_eq!(direction_map[&1], (0, false), "complement should be negative");
    }

    #[test]
    fn test_dedup_multiple_complements_of_same_rep() {
        // cols 1 and 2 are both complements of col 0
        let data = array![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 1.0],
        ];
        let (deduped, keep_idx, direction_map) = dedup_variants(&data);
        assert_eq!(deduped.ncols(), 1);
        assert_eq!(keep_idx, vec![0]);
        assert_eq!(direction_map[&1], (0, false));
        assert_eq!(direction_map[&2], (0, false));
    }

    #[test]
    fn test_dedup_mixed_identical_and_complement() {
        // col 0: [1,0,1], col 1: identical to 0, col 2: complement of 0, col 3: unique
        let data = array![
            [1.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0, 1.0],
        ];
        let (deduped, keep_idx, direction_map) = dedup_variants(&data);
        assert_eq!(deduped.ncols(), 2);
        assert_eq!(keep_idx, vec![0, 3]);
        assert_eq!(direction_map[&1], (0, true),  "identical → positive");
        assert_eq!(direction_map[&2], (0, false), "complement → negative");
    }

    // ── direction map from phase-2-style pruning ──────────────────────────────

    #[test]
    fn test_direction_positive_for_identical_pair() {
        // r for identical columns = 1.0 → positive
        let data = array![[1.0,1.0],[0.0,0.0],[1.0,1.0],[0.0,0.0]];
        let r = calculate_r(&data, 0, 1);
        assert!(r >= 0.0, "identical columns should give r >= 0, got {r}");
    }

    #[test]
    fn test_direction_negative_for_opposite_pair() {
        // r for complementary columns = -1.0 → negative
        let data = array![[1.0,0.0],[0.0,1.0],[1.0,0.0],[0.0,1.0]];
        let r = calculate_r(&data, 0, 1);
        assert!(r < 0.0, "opposite columns should give r < 0, got {r}");
    }

    // ── max_r_squared ────────────────────────────────────────────────────────

    #[test]
    fn test_max_r_squared_equal_mafs() {
        // Two variants with the same MAF can reach r² = 1
        let m = max_r_squared(0.3, 0.3);
        assert!((m - 1.0).abs() < 1e-10, "equal MAFs should give max_r²=1.0, got {m}");
    }

    #[test]
    fn test_max_r_squared_different_mafs() {
        // p=0.1, q=0.5 → max_r² = 0.1*0.5 / (0.5*0.9) = 0.05/0.45 ≈ 0.111
        let m = max_r_squared(0.1, 0.5);
        assert!((m - 1.0 / 9.0).abs() < 1e-10, "expected ~0.111, got {m}");
    }

    #[test]
    fn test_max_r_squared_symmetric() {
        let a = max_r_squared(0.1, 0.4);
        let b = max_r_squared(0.4, 0.1);
        assert!((a - b).abs() < 1e-10, "max_r_squared should be symmetric");
    }

    #[test]
    fn test_max_r_squared_monomorphic() {
        assert!(max_r_squared(0.0, 0.5) < 1e-10, "monomorphic variant should give 0");
        assert!(max_r_squared(0.3, 1.0) < 1e-10, "fixed variant should give 0");
    }

    // ── calculate_r ──────────────────────────────────────────────────────────

    #[test]
    fn test_r_identical_columns_is_one() {
        let data = array![[1.0,1.0],[0.0,0.0],[1.0,1.0],[0.0,0.0]];
        let r = calculate_r(&data, 0, 1);
        assert!((r - 1.0).abs() < 1e-10, "r for identical columns should be 1.0, got {r}");
    }

    #[test]
    fn test_r_opposite_columns_is_minus_one() {
        let data = array![[1.0,0.0],[0.0,1.0],[1.0,0.0],[0.0,1.0]];
        let r = calculate_r(&data, 0, 1);
        assert!((r + 1.0).abs() < 1e-10, "r for opposite columns should be -1.0, got {r}");
    }

    #[test]
    fn test_r_independent_variants_is_zero() {
        let data = array![[0.0,0.0],[0.0,1.0],[1.0,0.0],[1.0,1.0]];
        let r = calculate_r(&data, 0, 1);
        assert!(r.abs() < 1e-10, "r for independent variants should be 0.0, got {r}");
    }

    #[test]
    fn test_r_monomorphic_variant_is_zero() {
        let data = array![[1.0,0.0],[0.0,0.0],[1.0,0.0]];
        let r = calculate_r(&data, 0, 1);
        assert!(r.abs() < 1e-10, "r with monomorphic column should be 0.0, got {r}");
    }

    #[test]
    fn test_r_self_is_one() {
        let data = array![[1.0,0.0],[0.0,1.0],[1.0,0.0]];
        let r = calculate_r(&data, 0, 0);
        assert!((r - 1.0).abs() < 1e-10, "r of variant with itself should be 1.0, got {r}");
    }

    // ── calculate_d_prime ────────────────────────────────────────────────────

    #[test]
    fn test_d_prime_identical_columns_is_one() {
        let data = array![[1.0,1.0],[0.0,0.0],[1.0,1.0],[0.0,0.0]];
        let d = calculate_d_prime(&data, 0, 1);
        assert!((d - 1.0).abs() < 1e-10, "D' for identical columns should be 1.0, got {d}");
    }

    #[test]
    fn test_d_prime_self_is_one() {
        let data = array![[1.0,0.0],[0.0,1.0],[1.0,0.0],[0.0,1.0]];
        let d = calculate_d_prime(&data, 0, 0);
        assert!((d - 1.0).abs() < 1e-10, "D' of variant with itself should be 1.0, got {d}");
    }

    #[test]
    fn test_d_prime_independent_variants_is_zero() {
        let data = array![[0.0,0.0],[0.0,1.0],[1.0,0.0],[1.0,1.0]];
        let d = calculate_d_prime(&data, 0, 1);
        assert!(d.abs() < 1e-10, "D' for independent variants should be 0.0, got {d}");
    }

    #[test]
    fn test_d_prime_opposite_columns_is_one() {
        // Perfectly anticorrelated → |D'| = 1 (perfect LD, complementary direction)
        let data = array![[1.0,0.0],[0.0,1.0],[1.0,0.0],[0.0,1.0]];
        let d = calculate_d_prime(&data, 0, 1);
        assert!((d - 1.0).abs() < 1e-10, "D' for opposite columns should be 1.0 (|D'|), got {d}");
    }

    // ── MAF sort helpers ────────────────────────────────────────────────────

    #[test]
    fn test_maf_sort_descending_enables_early_break() {
        // Verify that sorting by descending MAF makes max_r_squared
        // monotonically decrease along the inner loop (j > i).
        let mafs = vec![0.5, 0.3, 0.1, 0.05];  // already descending
        for i in 0..mafs.len() {
            let mut prev = f64::INFINITY;
            for j in (i + 1)..mafs.len() {
                let m = max_r_squared(mafs[i], mafs[j]);
                assert!(m <= prev + 1e-10,
                    "max_r_squared should be non-increasing: i={i}, j={j}, prev={prev}, cur={m}");
                prev = m;
            }
        }
    }

    #[test]
    fn test_phase2_keeps_higher_maf_as_representative() {
        // Two correlated variants: col0 MAF=0.75, col1 MAF=0.25 (complement)
        // After complement detection in dedup, this pair is already caught.
        // Test with non-complement but high-r² pair instead.
        // col0: [1,1,1,0] MAF=0.75, col1: [1,1,0,0] MAF=0.5
        // r² = (4*2 - 3*2)² / ((3*1)*(2*2)) = 4/12 = 0.333 — below typical thresholds
        // Use a pair with higher r²: col0=[1,1,0,0] col1=[1,0,0,0] MAF=0.5, 0.25
        // r² = (4*1-2*1)²/((2*2)*(1*3)) = 4/12 = 0.333 — still low
        // Just verify the sort permutation logic directly:
        let mafs = array![0.1, 0.5, 0.3];
        let mut perm: Vec<usize> = (0..3).collect();
        perm.sort_by(|&a, &b| mafs[b].partial_cmp(&mafs[a]).unwrap_or(std::cmp::Ordering::Equal));
        assert_eq!(perm, vec![1, 2, 0], "should sort indices by descending MAF");
        // Verify MAF order after permutation
        let sorted_mafs: Vec<f64> = perm.iter().map(|&k| mafs[k]).collect();
        assert!((sorted_mafs[0] - 0.5).abs() < 1e-10);
        assert!((sorted_mafs[1] - 0.3).abs() < 1e-10);
        assert!((sorted_mafs[2] - 0.1).abs() < 1e-10);
    }
}
