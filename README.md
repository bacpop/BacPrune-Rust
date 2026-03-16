# BacPrune-Rust

BacPrune-Rust is a fast linkage disequilibrium (LD) pruning tool for haploid genotype matrices, written in Rust.
It is the primary, up-to-date implementation of BacPrune and is designed to handle large genomic datasets efficiently.

---

## Background

Earlier versions of BacPrune exist in other languages but are no longer actively developed:

- **R version** (`r_and_stan_versions/bac_prune.R`) — functionally correct but too slow for large datasets.
- **Stan version** (`r_and_stan_versions/bac_prune.Stan`) — an early experimental attempt; not a pruning tool.

---

## Installation

Requires [Rust](https://www.rust-lang.org/tools/install) (edition 2021, stable toolchain).

```bash
git clone <repo-url>
cd BacPrune-Rust
cargo build --release
```

The compiled binary will be at `target/release/bacprune`.

---

## Usage

```
bacprune <input_file> <n_rows> <n_cols> <maf_cutoff> <ld_threshold> <output_directory> [--r|--dprime]
bacprune <input_file> <n_rows> <n_cols> <maf_cutoff> <output_directory> --dedup
```

### Positional arguments

| Argument           | Description |
|--------------------|-------------|
| `input_file`       | Path to the input CSV file (see Input format below) |
| `n_rows`           | Total number of rows in the CSV **including the header row** |
| `n_cols`           | Number of columns in the CSV |
| `maf_cutoff`       | Minor allele frequency cutoff; variants with MAF below this are removed before LD pruning |
| `ld_threshold`     | LD threshold; variant pairs at or above this value are pruned (not used with `--dedup`) |
| `output_directory` | Directory to write output files into |

### Method flags (choose one)

| Flag        | Description |
|-------------|-------------|
| *(default)* or `--r`    | Prune by **Pearson r** threshold. Keeps the higher-MAF variant from each pair with \|r\| ≥ threshold. |
| `--dprime`  | Prune by **D'** threshold. Keeps the higher-MAF variant from each pair with \|D'\| ≥ threshold. |
| `--dedup`   | **Hash-based deduplication only.** Removes exact duplicate columns (identical genotype vectors) without any pairwise LD calculation. Much faster for a first-pass clean-up; `ld_threshold` is not required. |

### Examples

```bash
# Prune with Pearson r >= 0.8 (default)
bacprune genotypes.csv 604 1000 0.01 0.8 ./results

# Prune with D' >= 0.95
bacprune genotypes.csv 604 1000 0.01 0.95 ./results --dprime

# Remove exact duplicates only (no LD threshold needed)
bacprune genotypes.csv 604 1000 0.01 ./results --dedup
```

---

## Input format

The input file must be a **CSV with no row names** and numeric (float) values throughout, including the header row.
The header row should contain variant identifiers encoded as numbers (e.g. base-1 column indices).

- **Rows**: samples (individuals)
- **Columns**: variants (SNPs)
- **Values**: haploid genotype calls encoded as `0` (reference allele) or `1` (alternate allele)

Example (3 samples, 4 variants, with a numeric header):
```
1,2,3,4
0,1,0,1
1,1,0,0
0,0,1,1
```

---

## Output files

Both output files are written to `<output_directory>`.

| File | Description |
|------|-------------|
| `bacprune_rust_results.csv` | Pruned genotype matrix in the same format as the input (numeric header row included) |
| `ld_pruning_summary.csv`    | Two-column CSV recording which representative SNP was kept for each group of pruned SNPs (base-0 indexing) |

---

## Algorithm

Pruning runs in two phases:

1. **Exact duplicate removal (all modes):** Each variant column is hashed to a `Vec<u8>` key. The first occurrence of each unique genotype pattern is kept; all later identical columns are pruned. This is O(n·v) and performs no pairwise comparisons.

2. **LD threshold pruning (`--r` / `--dprime` only):** A greedy pairwise pass over the deduplicated variants. For each pair (i, j):
   - Compute \|r\| (Pearson correlation) or \|D'\| (Lewontin's normalised disequilibrium coefficient).
   - If the score ≥ `ld_threshold`, the variant with the **lower MAF** is pruned and the higher-MAF variant is recorded as the representative.

In both phases, pruned variants and their representatives are recorded in the summary CSV.

### LD metrics

**Pearson r (`--r`):**
The standard correlation coefficient between two binary vectors. r = 0 means the variants are independent; |r| = 1 means they are in perfect LD (identical or complementary). The pruning threshold is applied to |r|.

**D' (`--dprime`):**
Lewontin's normalised disequilibrium coefficient. D' = 0 means no LD; |D'| = 1 means perfect LD (no recombinant haplotypes observed). This implementation returns |D'|, so the value is always in [0, 1].

---

## Running tests

```bash
cargo test
```
