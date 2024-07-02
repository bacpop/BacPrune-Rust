// LD Pruning in Rust
// This code currently only prunes for LD=1, as lower-D' pruning isn't required for our GWAS

//libraries:
use ndarray::prelude::*;
use ndarray::OwnedRepr;
use csv::ReaderBuilder;
use ndarray::Array2;
use ndarray_rand::rand_distr::num_traits::ToPrimitive;
use std::fs::File;
use ndarray_csv::Array2Reader;
use std::cmp::Ordering;
use std::fs::OpenOptions;
use std::collections::{HashMap, HashSet};
use std::{env, ops::Div};
use std::path::Path;

//use rgsl::{
//    randist::t_distribution::{tdist_P, tdist_Q},
//    statistics::correlation,
//};

fn main() -> Result<(), csv::Error> {
    println!("Welcome to the LD Pruning Module.");
/* 
    use std::env;
    let key = "RUSTFLAGS";
    env::set_var(key, "/Users/lilyjacqueline/mambaforge/pkgs/gsl-2.7.1-hdbe807d_1/bin/gsl-config");
*/
    // Get the input file path, nrows (including header), and ncols from command line arguments
    let args: Vec<String> = env::args().collect(); // Changed: Reading command-line arguments.
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
    println!("{:?}", raw_gt_data);
    println!("Your data has been successfully read in. Sit tight while we run your analysis.");

    //turn "has headers" to false in read.csv function, copy header into seperate vector, then remove it from gt (just remove first row)
    let gt_header = raw_gt_data.select(Axis(0), &[0usize]);
    let raw_gt_data = raw_gt_data.slice(s![1..raw_gt_data.nrows(),..]).to_owned();
    
    //remove variants from header using skip index
    //concat header with final pruned gt data and export as usual
    //check that the headers match the variants by comparing header+SNP to original dataset's header+SNP (maybe by using match function to find identical columns?)

    //Calculate MAF
    let mafs = calc_maf(&raw_gt_data);
    println!("Minor allele frequencies: {:?}", mafs);
    println!("Minor allele frequencies were successfully calculated.");
    
    //Discard variants with MAF below cutoff
    //let filtered_gt_data = maf_prune(&raw_gt_data, &mafs, &cutoff);
    let (filtered_gt_data, keep_index) = maf_prune(&raw_gt_data, &mafs, &cutoff);
    let gt_header = gt_header.select(Axis(1), keep_index.as_slice());
    println!("Kept indices after MAF filtering: {:?}", keep_index);
    println!("Data were successfully filtered by MAF.");

    //Update MAF list after discarding variants
    let mafs = calc_maf(&filtered_gt_data);


    //Spearman's coefficient
    //Works for binary (SNP, AMR), continuous (E-test), and ordinal (MICs) data
    //fn correlationscore(gtdata:&Array2<f64>, phendata:&Array1<f64>, variant:usize) {
        //to use continuous data, must first rank
        //otherwise, it's fine to use binary AMR data or categorical ordinal MIC data directly (skip to next section)

        //calculate correlation coefficient
        //let disqrd:Array1<f64> = (&gtdata.slice(s![..,variant]) - phendata).iter().map(|x| x.powf(2.0)).collect();
        //let corrcoeff = 1-(((6usize*disqrd.sum().to_usize().unwrap()))/(disqrd.len().pow(3)-disqrd.len()));
        
        //calculate p value
        //use z_table::{lookup_with, reverse_lookup_with};
        //use rand_distr::{StudentT, Distribution};

        //let ttest = (corrcoeff*(disqrd.len()-2).sqrt())/((1-corrcoeff.pow(2)).sqrt());
        //let tcrit = 

        //TRYING rgsl stuff
        //let gtdat = gtdata.slice(s![..,variant]).to_vec();
        //let phendata = phendata.to_vec();

        //let r = correlation(&gtdat, 1, &phendata, 1, gtdata.len_of(Axis(0)));

        //let df = (gtdata.len_of(Axis(0)) - 2) as f64;
        //let statistic = df.sqrt() * r / (1.0 - r.powi(2)).sqrt();
        //let p_value:f64 = 2.0 * tdist_P(statistic, df).min(tdist_Q(statistic, df));
        //return p_value;

        //TRYING linear regression by hand



    //}

/* 
    //read in phenotype data
    let phenotype_data = read_csv("resistances.csv", 603, 2);
    let phenotype_data = phenotype_data.slice(s![..,1]);
*/

    //use correlation function and phenotype data data to calculate p values for each variant
    //for i in 0..=filtered_gt_data.nrows() {
    //    let pval = correlationscore(&filtered_gt_data.to_owned(), &phenotype_data.to_owned(), i);
    //    println!("P value for variant {:?} is: {:?}", i, pval)
    //}


    // If you wanted to ONLY prune out those with LD=1, you could do it with just the following steps:
    //1. Sort by MAF
    //let sorted_gt_data = sort_by_maf(&mafs, &filtered_gt_data);
    //let mafs = calc_maf(&sorted_gt_data);

    //2. For each pair of SNPs with same MAF:
    //2. Pairwise comparison of all individuals checking that there are at least one pair of indivs that is not perfectly 0/0 or 1/1; the second you find that pair, break
    //3. Else (meaning if that pair doesn't exist), then prune out the lower MAF sample (see SNPPrune paper pg.2)

    // create skip index (using HashSet package for a set)
    let mut skip_index: HashSet<usize> = HashSet::new();
    // create index to track SNPs that are pruned out and their representative SNP
    let mut rep_snps: HashMap<usize, Vec<usize>> = HashMap::new();
        
    //When a variant is pruned, its index number is added to this set
    //Before trying to compare two variants, the loop first checks that neither is in the index
    //This means that variants that have already been pruned will be skipped over during furture iterations,
    //without causing any indexing issues
    //The skip index also acts as the record of which variants should be pruned out of the dataset

    for i in 0..filtered_gt_data.ncols() {
        for j in i + 1..filtered_gt_data.ncols() {
            if (mafs[i] - mafs[j]).abs() < 1e-6 && (i != j) && !skip_index.contains(&i) && !skip_index.contains(&j) {
                let rowsums_ij = filtered_gt_data.select(Axis(1), &[i, j]).sum_axis(Axis(1));
                if rowsums_ij.iter().all(|&x| x == 0.0 || x == 2.0) {
                    println!("Pruning column {} due to perfect LD with column {}", j, i);
                    // SNPs i and j are in perfect LD, prune the one with the second one
                    skip_index.insert(j); // Prune the second SNP
                    rep_snps.entry(i).or_insert_with(Vec::new).push(j); // add pruned SNP to representative index
                }
            }
        }
    }
    
    //turn skip index into keep index (so can use in pruning .select() function)
    let keep_index: Vec<usize> = (0..filtered_gt_data.ncols()).filter(|x| !skip_index.contains(x)).collect();
    println!("Kept indices after LD pruning phase 1: {:?}", keep_index);

    //LD PRUNE PHASE 1
    // Prune the LD=1 variants out!
    let ldbelow1_gt_data = filtered_gt_data.select(Axis(1), keep_index.as_slice());
    println!("LD pruning phase 1 has been completed.");

/*
   //If user wants to prune for LD=1 cases, program ends here.
   //(Note that the program doesn't prune ALL LD=1 cases, just some)
   //If user wants to prune for LD<1 cases, continue:
   //(with the reduced dataset) run everything after this:

    //Next: make it so you do not compare samples with very different MAFs
    //Maybe by sorting by MAF and only iterating over varBs that are within a certain range of varA

    //Update MAF list after discarding variants
    let mafs = calc_maf(&ldbelow1_gt_data);

    //Make another skip/prune index
    let mut prune_index: HashSet<usize> = HashSet::new();
    //Set LD threshold
    let ld_threshold = 0.95f64;

    for i in 0..ldbelow1_gt_data.ncols() {
        for j in i..ldbelow1_gt_data.ncols() {
            //add if loop here with SNPrune calculation of how close mafs have to be to run rest
            //quick and dirty method here: if (mafs[i] - mafs[j]).abs() <= 0.05 
            if i != j && !prune_index.contains(&i) && !prune_index.contains(&j) && (mafs[i] - mafs[j]).abs() <= 0.1 {
                //here you have to calculate the AFs, haplotypes freqs, etc. for the d_prime_score function
                //would be more efficient to calculate these once, and reference that matrix
                let d_prime_score = calculate_d_prime(&ldbelow1_gt_data, i, j);
                if d_prime_score >= ld_threshold {
                    prune_index.insert(i);
                }
            }
        }
        println!("{:?} of 198248 variants have been completed.", i);
        println!("That's {:?}% complete.", i.div(198248));
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

*/
    //Write results to .csv
    // Make sure to include column names (headers)
    // could maybe use .skip(1) to skip the first row throughout the program?
    //full_prune_gt_data.records();


    
    //ADD HEADER BACK IN
    let gt_header = gt_header.select(Axis(1), keep_index.as_slice());
    let ldbelow1_gt_data = ndarray::concatenate![Axis(0), gt_header, ldbelow1_gt_data];


    let string_arr = ldbelow1_gt_data.map(|e| e.to_string());
    //let string_arr = full_prune_gt_data.map(|e| e.to_string());
    println!("String array: {:?}", string_arr);

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

/*
    for i in 0..full_prune_gt_data.nrows() {
        wtr.write_record(&string_arr.slice(s![i, ..])).expect("Error in writing to .csv");
    }
*/

    wtr.flush()?;

    // Write representative SNPs and pruned SNPs to a new .csv
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open("ld_pruning_summary.csv")
        .unwrap();
    let mut wtr = csv::Writer::from_writer(file);
    wtr.write_record(&["Representative SNP", "Pruned SNPs"]).expect("Error writing header to CSV");
    for (rep_snp, pruned_snps) in rep_snps {
        let rep_snp_str = rep_snp.to_string();
        let pruned_snps_str = pruned_snps.iter().map(|snp| snp.to_string()).collect::<Vec<String>>().join(", ");
        wtr.write_record(&[rep_snp_str, pruned_snps_str]).expect("Error writing record to CSV");
    }
    wtr.flush()?;

    Ok(())

}

// MAIN FUNCTION DEFINITIONS

fn read_csv(path_to_file: &str, n_rows:usize, n_cols:usize) -> Array2<f64> {
    let file = File::open(path_to_file).expect("File not found :(");
    let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);
    reader.deserialize_array2::<f64>((n_rows,n_cols)).expect("Failed to unwrap .csv file.")
}

fn calc_maf(data:&Array2<f64>) -> Array1<f64> {
    let num_individuals = data.nrows() as f64;
    println!("Number of rows in calc_maf (should equal 10): {:?}", num_individuals);
    let calcmafs: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 1]>> = data.sum_axis(Axis(0));
    let calcmafs: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 1]>> = calcmafs.div(num_individuals);
    return calcmafs;
}

fn maf_prune(data:&Array2<f64>, mafs:&Array1<f64>, cutoff:&f64) -> (Array2<f64>, Vec<usize>) {
    // Find index of each column with a MAF at or above the cutoff
    let keep_these_index = mafs
        .into_iter()
        .enumerate()
        .filter(|(_, x)| x >= &cutoff)
        .map(|(index, _)| index)
        .collect::<Vec<_>>();

    let filtered_data = data.select(Axis(1), &keep_these_index);
    return (filtered_data, keep_these_index);
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
