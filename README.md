# BacPrune-Rust

BacPrune-Rust performs linkage disequilibrium (LD) pruning for haploid genotype matrices using D' and r scores. It was written in Rust and designed to scale to matrices with millions of variants, such as the CRyPTIC Consortium _M. tuberculosis_ dataset.

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

```bash
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

Three modes are available. D' is recommended for GWAS applications as it prunes both positive and negative correlations.

```
bacprune <input_file> <n_rows> <n_cols> <maf_cutoff> <output_directory> --ld <threshold> [--r|--dprime]
bacprune <input_file> <n_rows> <n_cols> <maf_cutoff> <output_directory> --dedup
```

| Argument | Description |
|----------|-------------|
| `input_file` | Path to input CSV |
| `n_rows` | Number of rows **including the header row** |
| `n_cols` | Number of columns |
| `maf_cutoff` | Variants with allele frequency below this are removed |
| `output_directory` | Directory for output files |
| `--ld <threshold>` | Prune pairs at or above this score (required for `--r`/`--dprime`) |

| Flag | Method | `--ld` required? |
|------|--------|-----------------|
| `--dprime` | \|D'\| — default | Yes |
| `--r` | \|Pearson r\| | Yes |
| `--dedup` | Exact hash match only | No |

Run `bacprune --help` for full usage.

### Input format

A CSV where every value is numeric (f64). The first row is a header of variant identifiers; subsequent rows are samples, with values `0` (reference) or `1` (alternate allele).

```
1,2,3,4
0,1,0,1
1,1,0,0
0,0,1,1
```

`n_rows` includes the header (3 samples → `n_rows = 4`).

> [!NOTE]
> The MAF filter operates on the frequency of the allele encoded as `1`. This means it removes variants where the alternate allele is rare, but will not remove variants where the *reference* allele is rare (i.e. high-frequency alternates). To filter out both rare ALT and rare REF alleles, normalise your input so that the alternate allele is always `1` and the reference allele is always `0` before running BacPrune.

### Examples

```bash
bacprune genotypes.csv 604 1000 0.01 ./results --ld 0.95           # D' (default)
bacprune genotypes.csv 604 1000 0.01 ./results --ld 0.8 --r        # Pearson r
bacprune genotypes.csv 604 1000 0.01 ./results --dedup             # exact duplicates only
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
