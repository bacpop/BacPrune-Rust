// Rust prune 4

//NOTE TO SLEF:
//Remember that Rust indexing starts at 0 - make sure that all frequencies etc. are correctly indexed!!

//libraries:
use std::{error::Error, io, process, ops::Div};
use ndarray_rand::rand_distr::num_traits::zero;
use polars::{prelude::*, frame::row::Row};
//extern crate csv;
use ndarray::prelude::*;
use ndarray::OwnedRepr;

use csv::ReaderBuilder;
use ndarray::Array2;
use std::fs::File;
use ndarray_csv::Array2Reader;

use std::cmp::Ordering;


fn main() {
    println!("Welcome to the LD Pruning Module.");

    // Read in data
    let raw_gt_data = read_csv("3000_gts.csv");
    println!("{:?}", raw_gt_data);
    println!("Your data has been successfully read in. Sit tight while we run your analysis.");

    //To-do
    // Automatically find the number of rows and columns (samples and variants) in the df
    // Use those numbers for further calculations
    let n_rows = raw_gt_data.nrows();
    let n_cols = raw_gt_data.ncols();

    //Calculate MAF
    let mafs = calc_maf(&raw_gt_data);
    println!("Minor allele frequencies were successfully calculated.");
    
    //Discard variants with MAF below cutoff
    let cutoff = 0.01f64;
    let filtered_gt_data = maf_prune(&raw_gt_data, &mafs, &cutoff);
    println!("Data were successfully filtered by MAF.");
    
    //Update V after discarding variants
    let n_cols = filtered_gt_data.ncols();

    //Update MAF list after discarding variants
    let mafs = calc_maf(&filtered_gt_data);

    // If you wanted to ONLY prune out those with LD=1, you could do it with just the following steps:
    //1. Sort by MAF
    let sorted_gt_data = sort_by_maf(&mafs, &filtered_gt_data);
    //2. For each pair of SNPs with same MAF:
    //2. Pairwise comparison of all individuals checking that there are at least one pair of indivs that is not perfectly 0/0 or 1/1; the second you find that pair, break
    //3. Else (meaning if that pair doesn't exist), then prune out the lower MAF sample (see SNPPrune paper pg.2)

    for i in 0..n_cols {
        for j in 0..n_cols {
            if mafs[i] == mafs[j] {
                let rowsums_ij = sorted_gt_data.select(Axis(1), &[i,j]).sum_axis(Axis(0));
                println!("{:?}", rowsums_ij);
                break;
                //sum row of the two individuals
                //if sum =0 or 2, homo
                //if sum =1, hetero, break
            }
        }
    }
    //println!("D prime matrix: {:?}", pruned_data)

    //If you wanted to be able to set a threshold other than just LD=1, you could run the above and then (with the reduced dataset) run everything after this:

    //Next: make it so you do not compare samples with very different MAFs
    //Maybe by sorting by MAF and only iterating over varBs that are within a certain range of varA

    //Generate pairwise matrix of D' values
    //let mut d_prime_matrix = Array::zeros((n_cols, n_cols));
    //for i in 0..n_cols {
    //    for j in 0..n_cols {
    //        let d_prime_score = calculate_d_prime(&filtered_gt_data, i, j);
    //        let got = std::mem::replace(&mut d_prime_matrix[[i,j]], d_prime_score);
    //    }
    //}
    //println!("D prime matrix: {:?}", d_prime_matrix)
    
    //LD Prune (return dataset minus the pruned variants)
    

}

// MAIN FUNCTION DEFINITIONS

fn read_csv(path_to_file: &str) -> Array2<f64> {
    //1000 cols, 604 rows (incl. headers)
    let file = File::open(path_to_file).expect("File not found :(");
    let mut reader = ReaderBuilder::new().has_headers(true).from_reader(file);
    reader.deserialize_array2::<f64>((603, 1000)).expect("Failed to unwrap .csv file.")
}

fn calc_maf(data:&Array2<f64>) -> Array1<f64> {
    let calcmafs: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 1]>> = data.sum_axis(Axis(0));
    let calcmafs: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 1]>> = calcmafs.div(603f64);
    return calcmafs;
}

fn maf_prune(data:&Array2<f64>, mafs:&Array1<f64>, cutoff:&f64) -> Array2<f64> {
    // Find index of each column with a MAF at or above the cutoff
    let keep_these_index = mafs
        .into_iter()
        .enumerate()
        .filter(|(_, x)| x >= &cutoff)
        .map(|(index, _)| index)
        .collect::<Vec<_>>();

    let filtered_data = data.select(Axis(1), &keep_these_index);
    return filtered_data;
}

fn sort_by_maf(mafs:&Array1<f64>, data:&Array2<f64>) -> Array2<f64> {
    //Find index of sorted mafs
    let sorted_mafs_index = sort(mafs);
    //Sort data columns by index
    let sorted_data = data.select(Axis(1), sorted_mafs_index.as_slice());
    return sorted_data;

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

fn calculate_d_prime(data: &ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>, varianta:usize, variantb:usize) -> f64 {
    //Calculate the allele frequencies
    let allele_frequencies: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 2]>> = find_allele_frequencies(&data);

    //Calculate haplotype frequencies for  given pair of variants
    let haplotype_frequencies = find_haplotype_frequencies(&data, varianta, variantb);

    //Calculate D prime score for a given pair of variants
    let d_prime = calc_d_prime(&allele_frequencies, &haplotype_frequencies, varianta, variantb);
    return d_prime;

    // FUNCTION DEFINITIONS
    
    fn find_allele_frequencies(data:&Array2<f64>) -> Array2<f64> {
        let allele_frequencies1 = data.t().sum_axis(Axis(1)); //col sums
        let allele_frequencies0 = 603f64 - allele_frequencies1.clone();
        let allele_frequencies: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 2]>> = ndarray::stack![Axis(0), allele_frequencies0, allele_frequencies1];
        let allele_frequencies: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 2]>> = allele_frequencies.div(603f64);
        return allele_frequencies; //allele_freqs0 in first row, allele_freqs1 in second row
    }
    
    fn find_haplotype_frequencies(data:&Array2<f64>, varianta:usize, variantb:usize) -> Array1<f64> { // variantA/variantB is the index number of the variant
        // Create array (ab) containing just variant A and variant B from the dataset
        // To make this more efficient, could take values directly from the dataset for next steps
        // instead of creating new array
        let var_a = data.slice(s![.., varianta]);
        let var_b = data.slice(s![.., variantb]);
        let ab: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 2]>> = ndarray::stack![Axis(0), var_a, var_b];
    
        // Create arrays for comparisons in counting step (next step)
        let arr_onezero = Array::from_shape_vec((2).f(), vec![1.0,0.0])
            .unwrap();
        let arr_zeroone = Array::from_shape_vec((2).f(), vec![0.0,1.0])
            .unwrap();
    
        // Count each of the four possible haplotypes (00, 10, 01, 11)
        let bothzero = ab.axis_iter(Axis(1)).filter(|&x| x == Array::zeros(2)).count() as f64;
        let onezero = ab.axis_iter(Axis(1)).filter(|&x| x == arr_onezero).count() as f64;
        let zeroone = ab.axis_iter(Axis(1)).filter(|&x| x == arr_zeroone).count() as f64;
        let bothone = ab.axis_iter(Axis(1)).filter(|&x| x == Array::ones(2)).count() as f64;
    
        // Calculate haplotype frequencies using counts
        let haplotype_frequencies = ndarray::array![bothzero, onezero, zeroone, bothone];
        let haplotype_frequencies = haplotype_frequencies.div(603f64);
    
        return haplotype_frequencies;
    }

    //Calculate D prime score for a given pair of variants
    fn calc_d_prime(allele_freqs:&Array2<f64>, haplotype_freqs:&Array1<f64>, varianta:usize, variantb:usize) -> f64 {
        
        //Caclulate D score for given pair of variants
        // D = f(00)*f(11)-f(10)*f(01)
        let d_score = (haplotype_freqs[[0]]*haplotype_freqs[[3]]) - (haplotype_freqs[[1]]*haplotype_freqs[[2]]);

        if d_score < 0.0 {
            //Get D max set: -f(A)f(B) and -f(a)f(b)
            let d_max_set:Vec<f64> = vec![-1f64*(allele_freqs[[0,varianta]])*(allele_freqs[[0, variantb]]), -1f64*(allele_freqs[[1, varianta]])*(allele_freqs[[1, variantb]])];
            //Get Dmax value (max of D max set)
            let d_max = d_max_set.iter().max_by(|a,b| a.total_cmp(b)).expect("Oops");
            //Calculate D prime score
            let d_prime_score = d_score/d_max;
            return d_prime_score;

        } else if d_score > 0.0 {
            //Get D max set: f(A)f(b) and f(a)f(B)
            let d_max_set:Vec<f64> = vec![(allele_freqs[[0,varianta]])*(allele_freqs[[1, variantb]]), (allele_freqs[[1, varianta]])*(allele_freqs[[0, variantb]])];
            //Get Dmax value (min of D max set)
            let d_max = d_max_set.iter().min_by(|a,b| a.total_cmp(b)).expect("Oops");
            //Calculate D prime score
            let d_prime_score = d_score/d_max;
            return d_prime_score;

        } else {
            let d_prime_score = d_score;
            println!("D prime score should equal zero: {:?}", d_prime_score);
            return d_prime_score;
        }
    }

}
