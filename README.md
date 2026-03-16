# BacPrune-Rust

BacPrune-Rust performs linkage disequilibrium (LD) pruning for haploid genotype matrices using D' and r scores with optional additional filtering by minor allele frequency (MAF). Written in Rust and designed to scale to matrices with millions of variants and thousands of isolates, such as the complete CRyPTIC Consortium _M. tuberculosis_ dataset.

## Installation

### Dependencies

BacPrune-Rust is a self-contained compiled binary with no runtime dependencies. Building from source requires [Rust](https://www.rust-lang.org/tools/install) ≥ 1.78 (stable toolchain, edition 2021).

### Installing with conda (recommended)

BacPrune-Rust is available on the bioconda channel. Install into an existing environment:
```
conda install -c bioconda bacprune
```
or create a new dedicated environment (recommended):
```
conda create -n bacprune -c bioconda bacprune
```

### Building from source

To build from source (requires Rust ≥ 1.78, stable toolchain, edition 2021), run the following:

```
git clone https://github.com/bacpop/BacPrune-Rust
cd BacPrune-Rust
cargo build --release
```
The compiled binary will be located at `target/release/bacprune`.

You can also install the Python package from source using [pip](https://pip.pypa.io/en/stable/getting-started/) (requires Rust ≥ 1.78 and [maturin](https://github.com/PyO3/maturin)):

```
git clone https://github.com/bacpop/BacPrune-Rust
cd BacPrune-Rust
pip install maturin
pip install .
```

## Usage

BacPrune offers three modes: pruning by r² (Pearson r-squared), pruning by |D'| (Lewontin's D'), and pruning only variants with identical presence and absence (perfect positive correlation) which is determined by variant hashes. r² is the default and is commonly used for LD pruning in GWAS. D' is an alternative that captures both strong positive and strong negative correlations and may be preferred when the direction of correlation matters.

```
Usage: bacprune [OPTIONS] <INPUT_FILE> <N_ROWS> <N_COLS> <MAF_CUTOFF> <OUTPUT_DIRECTORY>

Arguments:
  <INPUT_FILE>        Path to the input CSV (numeric header row + sample rows of 0/1 genotypes)
  <N_ROWS>            Total number of rows in the CSV including the header row
  <N_COLS>            Number of columns in the CSV
  <MAF_CUTOFF>        Minor allele frequency cutoff; variants below this threshold are removed
  <OUTPUT_DIRECTORY>  Directory where output files are written

Options:
      --ld <THRESHOLD>  LD pruning threshold; pairs at or above this value are pruned. Required for --r and --dprime; not used with --dedup
      --r               Prune by r² (Pearson r-squared) threshold [default method]
      --dprime          Prune by |D'| (Lewontin's D') threshold
      --dedup           Remove exact duplicate variants only via hashing; no pairwise LD calculation
  -h, --help            Print help
  -V, --version         Print version
```

> [!TIP]
> `n_rows` includes the header row

### Formatting the genotype matrix

The genotype matrix should be in CSV format, with the first row being a header of (integer) variant IDs, and subsequent rows being samples with values `0` (reference) or `1` (alternate allele).

> [!NOTE]
> MAF filtering operates on the frequency of the allele encoded as `1`; thus, a variant with an allele encoded as `0` that is below the MAF threshold will never be filtered. If you wish to filter low-frequency reference alleles, normalize your genotype matrix so that the alternate allele is always `1` and the reference allele is always `0`.

### Example Commands

```
# r² (default) pruning for r²≥0.8, filter variants with MAF<5%
bacprune genotypes.csv 613 87092 0.05 ./results --ld 0.8

# D' pruning for |D'|=1, filter variants with MAF<5%
bacprune genotypes.csv 613 87092 0.05 ./results --ld 1 --dprime

# exact duplicates only, filter variants with MAF<1%
bacprune genotypes.csv 613 87092 0.01 ./results --dedup 
```

## Output files

All files are written to `<output_directory>`.

- `bacprune_rust_results.csv`: pruned genotype matrix in the same format as the input
- `ld_pruning_summary.csv`: one row per representative SNP, listing the post-MAF-filter column indices (base 0) of all variants pruned from its group
- `direction_of_correlation.csv`: one row per post-MAF variant; `Status` is `representative`, `positive_correlation` (identical genotype pattern, r > 0), or `negative_correlation` (complement pattern, r < 0)

## Algorithm

**Phase 1 (all modes):** each variant column is hashed as a `Vec<u8>` and inserted into a hash map. Exact duplicates are removed in a single O(n·v) pass. Complementary columns hash differently and are both retained.

**Phase 2 (`--r` and `--dprime` only):** greedy pairwise scan over remaining variants. For each pair exceeding the `--ld` threshold, the lower-MAF variant is pruned. The sign of r determines the correlation direction recorded in `direction_of_correlation.csv` for both modes.

**Pearson r²:** `r = (n·Σxy − Σx·Σy) / √(Σx·(n−Σx) · Σy·(n−Σy))` — threshold applied to r². The sign of r determines the correlation direction recorded in `direction_of_correlation.csv`.

**D':** `D = f(00)·f(11) − f(10)·f(01)`, normalised by D_max. Always returns |D'| ∈ [0, 1]. Correlation direction is determined separately via r since D' discards the sign.
