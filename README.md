# BacPrune-Rust

**BacPrune-Rust** is a fast linkage disequilibrium (LD) pruning tool for haploid genotype matrices, written in Rust.
It is the primary, up-to-date implementation of BacPrune and is designed to scale to large bacterial/haploid genomic datasets.

---

## Contents

- [Installation](#installation)
  - [pip](#pip)
  - [bioconda](#bioconda)
  - [Build from source](#build-from-source)
- [Usage](#usage)
- [Input format](#input-format)
- [Output files](#output-files)
- [Algorithm](#algorithm)
- [Running tests](#running-tests)
- [Other versions](#other-versions)

---

## Installation

### pip

Pre-built wheels are published to PyPI for Linux and macOS.
The `bacprune` executable is placed on `PATH` automatically.

```bash
pip install bacprune
```

To install from source (requires Rust ≥ 1.78 and [maturin](https://github.com/PyO3/maturin)):

```bash
pip install maturin
pip install .
```

### bioconda

```bash
conda install -c bioconda bacprune
```

### Build from source

Requires [Rust](https://www.rust-lang.org/tools/install) (stable toolchain, edition 2021).

```bash
git clone https://github.com/bacpop/BacPrune-Rust
cd BacPrune-Rust
cargo build --release
# binary is at target/release/bacprune
```

---

## Usage

```
bacprune <input_file> <n_rows> <n_cols> <maf_cutoff> <ld_threshold> <output_directory> [--r|--dprime]
bacprune <input_file> <n_rows> <n_cols> <maf_cutoff> <output_directory> --dedup
```

### Positional arguments

| Argument             | Description |
|----------------------|-------------|
| `input_file`         | Path to the input CSV (see [Input format](#input-format)) |
| `n_rows`             | Total number of rows in the CSV **including the header row** |
| `n_cols`             | Number of columns |
| `maf_cutoff`         | Minor allele frequency cutoff; variants below this threshold are removed before LD pruning |
| `ld_threshold`       | LD pruning threshold; pairs at or above this value are pruned (not required for `--dedup`) |
| `output_directory`   | Directory where output files are written |

### Method flags

Choose one; the default when no flag is given is `--r`.

| Flag       | LD metric | Pairwise? | `ld_threshold` required? |
|------------|-----------|-----------|--------------------------|
| `--r`      | \|Pearson r\| | Yes | Yes |
| `--dprime` | \|D'\| | Yes | Yes |
| `--dedup`  | Exact hash match | No (O(n·v)) | No |

### Examples

```bash
# Prune with |r| >= 0.8  (default method)
bacprune genotypes.csv 604 1000 0.01 0.8 ./results

# Prune with |r| >= 0.8  (explicit flag)
bacprune genotypes.csv 604 1000 0.01 0.8 ./results --r

# Prune with |D'| >= 0.95
bacprune genotypes.csv 604 1000 0.01 0.95 ./results --dprime

# Remove exact duplicate variants only (no LD threshold)
bacprune genotypes.csv 604 1000 0.01 ./results --dedup
```

---

## Input format

The input must be a **headerless-style CSV** where every value — including the first row — is numeric (f64).
The first row is treated as a header of variant identifiers (e.g. base-1 column indices); subsequent rows are samples.

- **Rows**: samples (individuals)
- **Columns**: variants (SNPs / loci)
- **Values**: haploid genotype calls encoded as `0` (reference allele) or `1` (alternate allele)

Example (3 samples, 4 variants, numeric header):
```
1,2,3,4
0,1,0,1
1,1,0,0
0,0,1,1
```

`n_rows` must include the header row (so a file with 3 samples has `n_rows = 4`).

---

## Output files

All files are written to `<output_directory>`.

### `bacprune_rust_results.csv`

The pruned genotype matrix in the same format as the input (numeric header row + sample rows, one variant per column).

### `ld_pruning_summary.csv`

Two-column CSV recording the representative kept for each group of redundant variants.
Indices are in post-MAF-filter column space (base 0).

| Column | Description |
|--------|-------------|
| `Representative SNP (base 0 indexing)` | Column index of the variant that was kept |
| `Pruned SNPs (base 0 indexing)` | Comma-separated column indices of all variants pruned from this group |

### `direction_of_correlation.csv`

One row per variant that passed the MAF filter (both kept and pruned).
Records whether each pruned variant was positively or negatively correlated with its representative.

| Column | Description |
|--------|-------------|
| `Variant` | Variant identifier (from the input header row) |
| `Status` | `representative`, `positive_correlation`, or `negative_correlation` |
| `Representative Variant` | Identifier of the representative; equals `Variant` for representatives |

`positive_correlation` means the pruned variant's genotype pattern matches its representative (r > 0, i.e. identical allele calls).
`negative_correlation` means the pruned variant is the complement of its representative (r < 0, i.e. allele calls are flipped at every sample).

---

## Algorithm

Pruning runs in two phases:

### Phase 1 — exact duplicate removal (all modes)

Each variant column is encoded as a `Vec<u8>` key and inserted into a hash map.
The first occurrence of each unique genotype pattern is kept; all later identical columns are recorded as pruned with `positive_correlation`.
This is **O(n·v)** where n = samples and v = variants — no pairwise comparisons are performed.

Note: complementary columns (one is the bitwise NOT of the other) encode to different keys and are both retained by phase 1.

### Phase 2 — LD threshold pruning (`--r` and `--dprime` only)

A greedy pairwise scan over the variants remaining after phase 1.
For each pair (i, j):

1. Compute |r| (Pearson r) or |D'| (Lewontin's D').
2. If the score ≥ `ld_threshold`, the variant with the **lower MAF** is pruned; the higher-MAF variant is kept as the representative.
3. The sign of r (regardless of whether `--r` or `--dprime` was used) determines the correlation direction written to `direction_of_correlation.csv`.

### LD metrics

**Pearson r (`--r`, default)**

Standard correlation coefficient between two binary vectors:

```
r = (n·Σxy − Σx·Σy) / √(Σx·(n−Σx) · Σy·(n−Σy))
```

r = 0 → independent; |r| = 1 → perfect LD (identical or complementary).
The pruning threshold is applied to |r|.

**D' (`--dprime`)**

Lewontin's normalised disequilibrium coefficient:

```
D  = f(00)·f(11) − f(10)·f(01)
D' = D / D_max
```

where D_max is the maximum possible |D| given the observed allele frequencies.
This implementation returns |D'|, so results are always in [0, 1].
D' = 0 → no LD; D' = 1 → perfect LD (identical or complementary columns).
The correlation direction (positive / negative) is determined separately using r, since D' discards the sign.

---

## Running tests

```bash
cargo test
```

---

## Other versions

Earlier implementations are kept in `r_and_stan_versions/` for reference but are no longer actively developed:

| Version | Notes |
|---------|-------|
| **R** (`bac_prune.R`) | Functionally correct but too slow for large datasets |
| **Stan** (`bac_prune.Stan`) | Early experimental attempt; not a pruning tool |
