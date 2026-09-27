# DspikeIn 1.3.2

## Bug fixes
- `validate_spikein_clade()`: fixed the `could not find function "DNAStringSet"`
  error in R CMD check (BioC 3.24). Alignment now uses `DECIPHER::AlignSeqs()`
  instead of `msa::msa()`, whose S4 constructor calls `DNAStringSet()`
  unqualified and fails when Biostrings is not attached.
- `validate_spikein_clade()`: terminal branch lengths are now matched to the
  correct tips (previously taken in edge order, not tip order).
- `validate_spikein_clade()`: the reported clade bootstrap is now the support
  of the sample clade's MRCA node rather than the maximum support in the tree;
  support is computed on bipartitions because NJ trees are unrooted.
- `validate_spikein_clade()`: edge-length annotations are now drawn with
  `ape::edgelabels()` (they were previously drawn on nodes).

## Improvements
- `validate_spikein_clade()`: input validation, unique tip labels, tree rooted
  on the first reference, plots recorded off-screen (no stray devices or
  `Rplots.pdf`), new `node_support` and `branch_table` outputs.
- Removed the `msa` dependency (one fewer package in Imports).

# DspikeIn 1.3.1

# DspikeIn 1.3.0

# DspikeIn 1.2.0

# DspikeIn 1.1.0

# DspikeIn 0.99.29
 
# DspikeIn 0.99.28
- Bump version 

# DspikeIn 0.99.27
- Bump version 


# DspikeIn 0.99.26
- Bump version 

# DspikeIn 0.99.25
- Bump version 


# DspikeIn 0.99.24
- Bump version 


# DspikeIn 0.99.23
- merging

# DspikeIn 0.99.22
- merging

# DspikeIn 0.99.21

# DspikeIn 0.99.20
- function added
- reviewer's comments implemented 

# Version: 0.99.11
- Adding TSE Examples
- Including TSE vignettes

# DspikeIn 0.99.10

_No major changes recorded._

# DspikeIn 0.99.9

- Cleaned up Rcheck.

# DspikeIn 0.99.9

- Documented both `physeq` and `tse` in `R/data.R` with proper `@format`, `@source`, and usage examples.
- Fixed system hostname configuration to remove “unknown” from shell prompt.
- Verified that vignettes build cleanly and only required `doc/` appears in tarball.
- Cleaned up redundant or empty directories during build to pass Bioconductor checks.

# DspikeIn 0.99.8

## Changes

- Streamlined the package by **removing non-essential helper and pipeline functions**.
- Retained **core spike-in quantification methods** only: validation, scaling, bias correction, and absolute abundance estimation.
- Removed dependencies on `DESeq2`, `edgeR`, and downstream differential analysis tools.
- Updated `DESCRIPTION`, `NAMESPACE`, and internal documentation to reflect core-only functionality.
- Added clarifying notes in the vignette and package description about the full pipeline availability on the `MGhotbi` branch.
- Reorganized internal structure for improved reproducibility and modular design.
- Updated `R` version dependency to `R (>= 4.5.0)` to comply with Bioconductor guidelines.
- Refined `biocViews` to match core functionality focus and removed unused categories.

# DspikeIn 0.99.7

_No major changes recorded._

# DspikeIn 0.99.6

## Changes

- Fixed error caused by missing `geom_star()` in `ggstar`.
- Updated example sections to conditionally load required packages.
- Corrected object reconstruction logic in `random_subsample_WithReductionFactor()`.
- Replaced deprecated functions and improved Bioconductor compliance.
- Added proper LICENSE formatting to eliminate NOTE.
- Documented dataset provenance and processing steps in `inst/scripts/data_preparation.md`, including SRA BioProject accessions for embedded data.

# DspikeIn 0.99.0

## Initial submission

- Submitted package to Bioconductor.
- Added spike-in analysis workflows.
- Included synthetic and real microbiome datasets.
