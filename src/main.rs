// Rust prune 4

//NOTE TO SLEF:
//Remember that Rust indexing starts at 0 - make sure that all frequencies etc. are correctly indexed!!

//libraries:
use std::{error::Error, io, process, ops::Div};
use polars::{prelude::*, frame::row::Row};
//extern crate csv;
use ndarray::prelude::*;
use ndarray::OwnedRepr;

use csv::ReaderBuilder;
use ndarray::Array2;
use std::fs::File;
use ndarray_csv::Array2Reader;
use ndarray::Zip;

use ndarray_csv::Array2Writer;

use std::cmp::Ordering;
use std::fs::OpenOptions;

//use serde::Serialize;

use std::collections::HashSet;
use csv::Writer;
use csv::WriterBuilder;

 
fn main() -> Result<(), csv::Error> {
    println!("Welcome to the LD Pruning Module.");

    //to read in data, need to know the number of rows and columns in your data
    let n_rows = 603;
    let n_cols = 1000; //including header

    // Read in data
    let raw_gt_data = read_csv("3000_gts.csv", n_rows, n_cols);
    println!("{:?}", raw_gt_data);
    println!("Your data has been successfully read in. Sit tight while we run your analysis.");

    //To-do
    // Automatically find the number of rows and columns (samples and variants) in the df
    // Use those numbers for further calculations

    //Calculate MAF
    let mafs = calc_maf(&raw_gt_data);
    println!("Minor allele frequencies were successfully calculated.");
    
    //Discard variants with MAF below cutoff
    let cutoff = 0.01f64;
    let filtered_gt_data = maf_prune(&raw_gt_data, &mafs, &cutoff);
    println!("Data were successfully filtered by MAF.");

    //Update MAF list after discarding variants
    let mafs = calc_maf(&filtered_gt_data);

    // If you wanted to ONLY prune out those with LD=1, you could do it with just the following steps:
    //1. Sort by MAF
    let sorted_gt_data = sort_by_maf(&mafs, &filtered_gt_data);
    //2. For each pair of SNPs with same MAF:
    //2. Pairwise comparison of all individuals checking that there are at least one pair of indivs that is not perfectly 0/0 or 1/1; the second you find that pair, break
    //3. Else (meaning if that pair doesn't exist), then prune out the lower MAF sample (see SNPPrune paper pg.2)

    //Create skip index (using HashSet package for a set)
    let mut skip_index: HashSet<usize> = HashSet::new();
    //When a variant is pruned, its index number is added to this set
    //Before trying to compare two variants, the loop first checks that neither is in the index
    //This means that variants that have already been pruned will be skipped over during furture iterations,
    //without causing any indexing issues
    //The skip index also acts as the record of which variants should be pruned out of the dataset

    for i in 0..filtered_gt_data.ncols() {
        for j in i..filtered_gt_data.ncols() {
            //when i = j (are the same variant) in the for loops, the snp will prune itself out
            // hence why needed to add the && i!=j condition
            if mafs[i] == mafs[j] && i != j && !skip_index.contains(&i) && !skip_index.contains(&j) { //inside this loop is one snp vs one snp
                //sum row of the two SNPs
                let rowsums_ij = sorted_gt_data.select(Axis(1), &[i,j]).sum_axis(Axis(1));
                //if any rowsums are equal to 1, then the SNP pair is not in perfect LD (so break)
                //if no rowsums are 1 (meaning all rowsums are 0 or 2, aka all pairs of individuals are 00 or 11),
                //then prune one of the snps
                if rowsums_ij.iter().any(|&i| i==1.0) == true {
                    break;
                } else {
                    //add to skip index (to be pruned)
                    skip_index.insert(i);
                }
            }
        }
    }
    
    //turn skip index into keep index (so can use in pruning .select() function)
    //there is probably a better way to do this but this is quick and dirty
    let skip_index: Vec<_> = skip_index.into_iter().collect();
    let keep_index: Vec<usize> = (0..(sorted_gt_data.ncols()-1)).collect::<Vec<_>>().into_iter().filter(|x| skip_index.contains(x) == false).collect::<Vec<usize>>();

    //LD PRUNE PHASE 1
    // Prune the LD=1 variants out!
    let ldbelow1_gt_data = sorted_gt_data.select(Axis(1), keep_index.as_slice());
    println!("LD pruning phase 1 has been completed.");

   //If user wants to prune for LD=1 cases, program ends here.
   //(Note that the program doesn't prune ALL LD=1 cases, just some)
   //If user wants to prune for LD<1 cases, continue:
   //(with the reduced dataset) run everything after this:

    //Next: make it so you do not compare samples with very different MAFs
    //Maybe by sorting by MAF and only iterating over varBs that are within a certain range of varA

    //Make another skip/prune index
    let mut prune_index: HashSet<usize> = HashSet::new();
    //Set LD threshold
    let ld_threshold = 0.99f64;

    for i in 0..ldbelow1_gt_data.ncols() {
        for j in i..ldbelow1_gt_data.ncols() {
            if i != j && !prune_index.contains(&i) && !prune_index.contains(&j) {
                let d_prime_score = calculate_d_prime(&ldbelow1_gt_data, i, j);
                if d_prime_score >= ld_threshold {
                    prune_index.insert(i);
                }
            }
        }
    }

    //turn skip index into keep index (so can use in pruning .select() function)
    //there is probably a better way to do this but this is quick and dirty
    let prune_index: Vec<_> = prune_index.into_iter().collect();
    let keep_prune_index: Vec<usize> = (0..(ldbelow1_gt_data.ncols()-1)).collect::<Vec<_>>().into_iter().filter(|x| prune_index.contains(x) == false).collect::<Vec<usize>>();

    //LD PRUNE PHASE 2
    // Prune the LD>=threshold variants out!
    let full_prune_gt_data = ldbelow1_gt_data.select(Axis(1), &keep_prune_index.as_slice());
    println!("LD pruning phase 2 has been completed.");

    //I think it's pruning out every single variant? (instead of leaving one)

    //Write results to .csv
    // Make sure to include column names (headers)
    // could maybe use .skip(1) to skip the first row throughout the program?
    //full_prune_gt_data.records();
    
    //let headers = full_prune_gt_data.headers()?;

    //let output_path = "rust_results.csv";

    let string_arr = full_prune_gt_data.map(|e| e.to_string());
    println!("String array: {:?}", string_arr);

    //let mut file = Writer::from_path(output_path)?;
    
    //for i in 0..full_prune_gt_data.nrows() {
    //    file.write_record(&string_arr[[i,..]]);
    //}

    //let nrows = &full_prune_gt_data.nrows();
    //let ncols = &full_prune_gt_data.ncols();


    //let file = File::create("rust_results.csv")?;
    //let mut writer = WriterBuilder::new().has_headers(false).from_writer(file);
    //writer.serialize_array2(&full_prune_gt_data)?;


    let mut file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open("test.csv")
        .unwrap();
    let mut wtr = csv::Writer::from_writer(file);

    for i in 0..full_prune_gt_data.nrows() {
        wtr.write_record(&string_arr.slice(s![i, ..]));
    }

    wtr.flush()?;
    
    Ok(())

}

// MAIN FUNCTION DEFINITIONS

fn read_csv(path_to_file: &str, n_rows:usize, n_cols:usize) -> Array2<f64> {
    let file = File::open(path_to_file).expect("File not found :(");
    let mut reader = ReaderBuilder::new().has_headers(true).from_reader(file);
    reader.deserialize_array2::<f64>((n_rows,n_cols)).expect("Failed to unwrap .csv file.")
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
        let allele_frequencies0 = 603f64 - &allele_frequencies1;
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

