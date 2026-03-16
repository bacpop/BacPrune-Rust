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

### Installing with pip

```
pip install bacprune
bacprune --version
```

### Building from source

```
git clone https://github.com/bacpop/BacPrune-Rust
cd BacPrune-Rust
cargo build --release
```

The compiled binary will be located at `target/release/bacprune`.

### Installing Python package from source

To install the Python package from source (requires Rust ≥ 1.78 and [maturin](https://github.com/PyO3/maturin)):

```
git clone https://github.com/bacpop/BacPrune-Rust
cd BacPrune-Rust
pip install maturin
pip install .
```

## Usage

BacPrune offers three modes: pruning by D' score, pruning by Pearson r, and pruning only variants with identical presence and absence (perfect positive correlation). D' score is recommended for use in GWAS applications where LD pruning is used to reduce nonidentifiability arising from correlations between variants, as this will prune both strong positive and strong negative correlations.

```
# remove by D' or r threshold
bacprune <input_file> <n_rows> <n_cols> <maf_cutoff> <output_directory> --ld <threshold> [--r|--dprime]

# only remove variants with identical presence/absence
bacprune <input_file> <n_rows> <n_cols> <maf_cutoff> <output_directory> --dedup
```

> [!TIP]
> `n_rows` includes the header row

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
      --r               Prune by |Pearson r| threshold
      --dprime          Prune by |D'| (Lewontin's D') threshold [default method]
      --dedup           Remove exact duplicate variants only via hashing; no pairwise LD calculation
  -h, --help            Print help
  -V, --version         Print version
```

### Input format

A CSV of numeric values. The first row is a header of (numeric) variant identifiers; subsequent rows are samples, with values `0` (reference) or `1` (alternate allele).

```
1,2,3,4
0,1,0,1
1,1,0,0
0,0,1,1
```

> [!NOTE]
> The MAF filter operates on the frequency of the allele encoded as `1`. This means it removes variants where the alternate allele is rare, but will not remove variants where the *reference* allele is rare. To filter out both rare ALT and rare REF alleles, normalise your input so that the alternate allele is always `1` and the reference allele is always `0` before running BacPrune.

### Example Commands

```
bacprune genotypes.csv 613 87092 0.05 ./results --ld 1           # D' (default) pruning for D'=1, filter variants with MAF<5%
bacprune genotypes.csv 613 87092 0 ./results --ld 0.8 --r        # Pearson r pruning for |r|≥0.8, no MAF filtering
bacprune genotypes.csv 613 87092 0.01 ./results --dedup          # exact duplicates only, filter variants with MAF<1%
```

## Output files

All files are written to `<output_directory>`.

**`bacprune_rust_results.csv`** — pruned genotype matrix in the same format as the input.

**`ld_pruning_summary.csv`** — one row per representative SNP, listing the post-MAF-filter column indices (base 0) of all variants pruned from its group.

**`direction_of_correlation.csv`** — one row per post-MAF variant. `Status` is `representative`, `positive_correlation` (identical genotype pattern, r > 0), or `negative_correlation` (complement pattern, r < 0).

## Algorithm

**Phase 1 (all modes):** each variant column is hashed as a `Vec<u8>` and inserted into a hash map. Exact duplicates are removed in a single O(n·v) pass. Complementary columns hash differently and are both retained.

**Phase 2 (`--r` and `--dprime` only):** greedy pairwise scan over remaining variants. For each pair exceeding the `--ld` threshold, the lower-MAF variant is pruned. The sign of r determines the correlation direction recorded in `direction_of_correlation.csv`.

**Pearson r:** `r = (n·Σxy − Σx·Σy) / √(Σx·(n−Σx) · Σy·(n−Σy))` — threshold applied to |r|.

**D':** `D = f(00)·f(11) − f(10)·f(01)`, normalised by D_max. Always returns |D'| ∈ [0, 1]. Correlation direction is determined separately via r since D' discards the sign.
