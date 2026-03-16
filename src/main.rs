// LD Pruning in Rust
// This code currently only prunes for LD=1, as lower-D' pruning isn't required for our GWAS

//libraries:
use ndarray::prelude::*;
use ndarray::OwnedRepr;
use csv::ReaderBuilder;
use ndarray::Array2;
use std::fs::File;
use ndarray_csv::Array2Reader;
use std::cmp::Ordering;
use std::fs::OpenOptions;
use std::collections::{HashMap, HashSet};
use std::{env, ops::Div};
use std::path::Path;

fn main() -> Result<(), csv::Error> {
    println!("Welcome to the LD Pruning Module.");

    // Get the input file path, nrows (including header), and ncols from command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 6 {
        println!("Usage: {} <input_file> <n_rows> <n_cols> <maf_cutoff> <output_directory>", args[0]);
        return Ok(());
    }
    let input_file = &args[1];
    let n_rows: usize = args[2].parse().expect("Please provide a valid number for n_rows");
    let n_cols: usize = args[3].parse().expect("Please provide a valid number for n_cols");
    let cutoff: f64 = args[4].parse().expect("Please provide a MAF cutoff (variants with a minor allele frequency below this cutoff will be pruned out)");
    let outdir = &args[5];

    // Read in data
    let raw_gt_data = read_csv(input_file, n_rows, n_cols);
    println!("Your data has been successfully read in. Sit tight while we run your analysis.");

    // Extract header row (first row), then drop it from the data
    let gt_header = raw_gt_data.select(Axis(0), &[0usize]);
    let raw_gt_data = raw_gt_data.slice(s![1..raw_gt_data.nrows(), ..]).to_owned();

    //Calculate MAF
    let mafs = calc_maf(&raw_gt_data);
    println!("Minor allele frequencies were successfully calculated.");

    //Discard variants with MAF below cutoff
    let (filtered_gt_data, keep_index) = maf_prune(&raw_gt_data, &mafs, &cutoff);
    let gt_header = gt_header.select(Axis(1), keep_index.as_slice());
    println!("Data were successfully filtered by MAF.");

    //Update MAF list after discarding variants
    let mafs = calc_maf(&filtered_gt_data);

    // create skip index (using HashSet package for a set)
    let mut skip_index: HashSet<usize> = HashSet::new();
    // create index to track SNPs that are pruned out and their representative SNP
    let mut rep_snps: HashMap<usize, Vec<usize>> = HashMap::new();

    //When a variant is pruned, its index number is added to this set
    //Before trying to compare two variants, the loop first checks that neither is in the index
    //This means that variants that have already been pruned will be skipped over during future iterations,
    //without causing any indexing issues
    //The skip index also acts as the record of which variants should be pruned out of the dataset

    for i in 0..filtered_gt_data.ncols() {
        for j in i + 1..filtered_gt_data.ncols() {
            if (mafs[i] - mafs[j]).abs() < 1e-6 && !skip_index.contains(&i) && !skip_index.contains(&j) {
                let rowsums_ij = filtered_gt_data.select(Axis(1), &[i, j]).sum_axis(Axis(1));
                if rowsums_ij.iter().all(|&x| x == 0.0 || x == 2.0) {
                    // SNPs i and j are in perfect LD; keep i as the representative, prune j
                    skip_index.insert(j);
                    rep_snps.entry(i).or_insert_with(Vec::new).push(j);
                }
            }
        }
    }

    //turn skip index into keep index (so can use in pruning .select() function)
    let keep_index: Vec<usize> = (0..filtered_gt_data.ncols()).filter(|x| !skip_index.contains(x)).collect();

    //LD PRUNE PHASE 1
    // Prune the LD=1 variants out!
    let ldbelow1_gt_data = filtered_gt_data.select(Axis(1), keep_index.as_slice());
    println!("LD pruning phase 1 has been completed.");

    //ADD HEADER BACK IN
    let gt_header = gt_header.select(Axis(1), keep_index.as_slice());
    let ldbelow1_gt_data = ndarray::concatenate![Axis(0), gt_header, ldbelow1_gt_data];

    let string_arr = ldbelow1_gt_data.map(|e| e.to_string());

    // construct the full path to the results CSV file
    let csv_path = Path::new(outdir).join("bacprune_rust_results.csv");

    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open(csv_path)
        .unwrap();
    let mut wtr = csv::Writer::from_writer(file);

    for i in 0..ldbelow1_gt_data.nrows() {
        wtr.write_record(&string_arr.slice(s![i, ..])).expect("Error in writing to .csv");
    }

    wtr.flush()?;

    // Write representative SNPs and pruned SNPs to a new .csv
    let rep_snp_path = Path::new(outdir).join("ld_pruning_summary.csv");
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open(rep_snp_path)
        .unwrap();
    let mut wtr = csv::Writer::from_writer(file);
    wtr.write_record(&["Representative SNP (base 0 indexing)", "Pruned SNPs (base 0 indexing)"]).expect("Error writing header to CSV");
    for (rep_snp, pruned_snps) in rep_snps {
        let rep_snp_str = rep_snp.to_string();
        let pruned_snps_str = pruned_snps.iter().map(|snp| snp.to_string()).collect::<Vec<String>>().join(", ");
        wtr.write_record(&[rep_snp_str, pruned_snps_str]).expect("Error writing record to CSV");
    }
    wtr.flush()?;

    Ok(())
}

// MAIN FUNCTION DEFINITIONS

fn read_csv(path_to_file: &str, n_rows: usize, n_cols: usize) -> Array2<f64> {
    let file = File::open(path_to_file).expect("File not found :(");
    let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);
    reader.deserialize_array2::<f64>((n_rows, n_cols)).expect("Failed to unwrap .csv file.")
}

fn calc_maf(data: &Array2<f64>) -> Array1<f64> {
    let num_individuals = data.nrows() as f64;
    let calcmafs = data.sum_axis(Axis(0));
    calcmafs.div(num_individuals)
}

fn maf_prune(data: &Array2<f64>, mafs: &Array1<f64>, cutoff: &f64) -> (Array2<f64>, Vec<usize>) {
    // Find index of each column with a MAF at or above the cutoff
    let keep_these_index = mafs
        .into_iter()
        .enumerate()
        .filter(|(_, x)| x >= &cutoff)
        .map(|(index, _)| index)
        .collect::<Vec<_>>();

    let filtered_data = data.select(Axis(1), &keep_these_index);
    (filtered_data, keep_these_index)
}

fn sort_by_maf(mafs: &Array1<f64>, data: &Array2<f64>) -> Array2<f64> {
    let sorted_mafs_index = sort(mafs);
    return data.select(Axis(1), sorted_mafs_index.as_slice());

    // Define function used to sort MAF vector and find its index
    fn sort(arr: &Array1<f64>) -> Vec<usize> {
        let mut out = (0..arr.len()).collect::<Vec<usize>>();
        out.sort_by(|&a_idx, &b_idx| {
            let a = arr[a_idx];
            let b = arr[b_idx];
            match (a.is_nan(), b.is_nan()) {
                (true, true) => Ordering::Equal,
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                (false, false) => a.partial_cmp(&b).unwrap(),
            }
        });
        out
    }
}

fn calculate_d_prime(data: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>, varianta: usize, variantb: usize) -> f64 {
    let allele_frequencies = find_allele_frequencies(data);
    let haplotype_frequencies = find_haplotype_frequencies(data, varianta, variantb);
    return calc_d_prime(&allele_frequencies, &haplotype_frequencies, varianta, variantb);

    fn find_allele_frequencies(data: &Array2<f64>) -> Array2<f64> {
        let n = data.nrows() as f64;
        let allele_frequencies1 = data.t().sum_axis(Axis(1)); // col sums = count of allele 1
        let allele_frequencies0 = n - &allele_frequencies1;   // count of allele 0
        let allele_frequencies = ndarray::stack![Axis(0), allele_frequencies0, allele_frequencies1];
        allele_frequencies.div(n) // allele_freqs0 in first row, allele_freqs1 in second row
    }

    fn find_haplotype_frequencies(data: &Array2<f64>, varianta: usize, variantb: usize) -> Array1<f64> {
        let n = data.nrows() as f64;
        let var_a = data.slice(s![.., varianta]);
        let var_b = data.slice(s![.., variantb]);
        let ab = ndarray::stack![Axis(0), var_a, var_b];

        let arr_onezero = Array::from_shape_vec((2).f(), vec![1.0, 0.0]).unwrap();
        let arr_zeroone = Array::from_shape_vec((2).f(), vec![0.0, 1.0]).unwrap();

        // Count each of the four possible haplotypes (00, 10, 01, 11)
        let bothzero = ab.axis_iter(Axis(1)).filter(|&x| x == Array::zeros(2)).count() as f64;
        let onezero  = ab.axis_iter(Axis(1)).filter(|&x| x == arr_onezero).count() as f64;
        let zeroone  = ab.axis_iter(Axis(1)).filter(|&x| x == arr_zeroone).count() as f64;
        let bothone  = ab.axis_iter(Axis(1)).filter(|&x| x == Array::ones(2)).count() as f64;

        let haplotype_frequencies = ndarray::array![bothzero, onezero, zeroone, bothone];
        haplotype_frequencies.div(n)
    }

    fn calc_d_prime(allele_freqs: &Array2<f64>, haplotype_freqs: &Array1<f64>, varianta: usize, variantb: usize) -> f64 {
        // D = f(00)*f(11) - f(10)*f(01)
        let d_score = (haplotype_freqs[[0]] * haplotype_freqs[[3]]) - (haplotype_freqs[[1]] * haplotype_freqs[[2]]);

        if d_score < 0.0 {
            // D_max = max(-f(A)*f(B), -f(a)*f(b))
            let d_max_set = vec![
                -1.0 * allele_freqs[[0, varianta]] * allele_freqs[[0, variantb]],
                -1.0 * allele_freqs[[1, varianta]] * allele_freqs[[1, variantb]],
            ];
            let d_max = d_max_set.iter().max_by(|a, b| a.total_cmp(b)).expect("Oops");
            d_score / d_max

        } else if d_score > 0.0 {
            // D_max = min(f(A)*f(b), f(a)*f(B))
            let d_max_set = vec![
                allele_freqs[[0, varianta]] * allele_freqs[[1, variantb]],
                allele_freqs[[1, varianta]] * allele_freqs[[0, variantb]],
            ];
            let d_max = d_max_set.iter().min_by(|a, b| a.total_cmp(b)).expect("Oops");
            d_score / d_max

        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── calc_maf ─────────────────────────────────────────────────────────────

    #[test]
    fn test_calc_maf_basic() {
        // 3 samples × 3 variants
        // col 0: all 1s  → MAF = 1.0
        // col 1: one 1   → MAF = 1/3
        // col 2: all 0s  → MAF = 0.0
        let data = array![
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        let mafs = calc_maf(&data);
        assert!((mafs[0] - 1.0).abs() < 1e-10);
        assert!((mafs[1] - 1.0 / 3.0).abs() < 1e-10);
        assert!((mafs[2] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_calc_maf_uses_nrows() {
        // MAF must be relative to the actual number of rows, not a hardcoded constant
        let data = array![[1.0, 0.0], [0.0, 1.0]];
        let mafs = calc_maf(&data);
        assert!((mafs[0] - 0.5).abs() < 1e-10);
        assert!((mafs[1] - 0.5).abs() < 1e-10);
    }

    // ── maf_prune ────────────────────────────────────────────────────────────

    #[test]
    fn test_maf_prune_removes_low_maf() {
        // col 0 MAF=0.5, col 1 MAF=0.0 → col 1 pruned
        let data = array![[1.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 0.0]];
        let mafs = array![0.5, 0.0];
        let (filtered, keep_idx) = maf_prune(&data, &mafs, &0.01);
        assert_eq!(filtered.ncols(), 1);
        assert_eq!(keep_idx, vec![0]);
    }

    #[test]
    fn test_maf_prune_keeps_at_cutoff() {
        // A variant with MAF exactly equal to the cutoff should be kept
        let data = array![[1.0, 0.0], [0.0, 0.0]];
        let mafs = array![0.5, 0.01];
        let (filtered, keep_idx) = maf_prune(&data, &mafs, &0.01);
        assert_eq!(filtered.ncols(), 2);
        assert_eq!(keep_idx, vec![0, 1]);
    }

    #[test]
    fn test_maf_prune_returns_correct_keep_index() {
        let data = array![[1.0, 0.0, 1.0], [0.0, 0.0, 0.0]];
        let mafs = array![0.5, 0.0, 0.5];
        let (filtered, keep_idx) = maf_prune(&data, &mafs, &0.01);
        assert_eq!(filtered.ncols(), 2);
        assert_eq!(keep_idx, vec![0, 2]);
    }

    // ── sort_by_maf ──────────────────────────────────────────────────────────

    #[test]
    fn test_sort_by_maf_ascending() {
        let data = array![
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ];
        // col 0 → maf 0.9, col 1 → maf 0.1, col 2 → maf 0.5
        // sorted ascending: col1, col2, col0
        let mafs = array![0.9, 0.1, 0.5];
        let sorted = sort_by_maf(&mafs, &data);
        // After sorting the first column should be what was col 1
        assert_eq!(sorted.column(0).to_vec(), data.column(1).to_vec());
        assert_eq!(sorted.column(1).to_vec(), data.column(2).to_vec());
        assert_eq!(sorted.column(2).to_vec(), data.column(0).to_vec());
    }

    #[test]
    fn test_sort_by_maf_already_sorted() {
        let data = array![[0.0, 1.0, 1.0], [0.0, 0.0, 1.0]];
        let mafs = array![0.0, 0.5, 1.0];
        let sorted = sort_by_maf(&mafs, &data);
        assert_eq!(sorted, data);
    }

    // ── calculate_d_prime ────────────────────────────────────────────────────

    #[test]
    fn test_d_prime_identical_columns_is_one() {
        // Two identical columns are in perfect LD → D' = 1
        let data = array![
            [1.0, 1.0],
            [0.0, 0.0],
            [1.0, 1.0],
            [0.0, 0.0],
        ];
        let d = calculate_d_prime(&data, 0, 1);
        assert!((d - 1.0).abs() < 1e-10, "D' should be 1.0, got {d}");
    }

    #[test]
    fn test_d_prime_self_is_one() {
        // A variant compared with itself is always D' = 1
        let data = array![
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ];
        let d = calculate_d_prime(&data, 0, 0);
        assert!((d - 1.0).abs() < 1e-10, "D' of variant with itself should be 1.0, got {d}");
    }

    #[test]
    fn test_d_prime_independent_variants_is_zero() {
        // All four haplotypes equally represented → D = 0 → D' = 0
        let data = array![
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [1.0, 1.0],
        ];
        let d = calculate_d_prime(&data, 0, 1);
        assert!(d.abs() < 1e-10, "D' for independent variants should be 0.0, got {d}");
    }

    #[test]
    fn test_d_prime_opposite_columns_is_one() {
        // col 1 = NOT col 0 → perfectly anticorrelated, but still perfect LD.
        // This implementation returns |D'| (always ≥ 0), so the result is 1.0.
        // Alleles are fully predictable from each other, just complementary.
        let data = array![
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ];
        let d = calculate_d_prime(&data, 0, 1);
        assert!((d - 1.0).abs() < 1e-10, "D' for opposite columns should be 1.0 (|D'|), got {d}");
    }
}
