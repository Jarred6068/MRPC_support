# README: HiC Chromatin Interaction Analysis Tools

This collection of R functions enables the analysis of chromatin interactions for eQTL trios (SNP-cis gene-trans gene) using Hi-C data.

## Overview

These functions are designed to analyze potential chromatin interactions between genetic variants and their associated cis and trans genes. The primary goal is to determine if physical chromatin interactions might explain trans-eQTL relationships discovered in GTEx data.

## Key Functions

### Helper Functions

- `loadRData()`: Loads an RData file and returns its contents
- `find_trans()`: Retrieves attributes for trans genes
- `find_cis()`: Retrieves attributes for cis genes

### Core Analysis Functions

- `get_trio_attr()`: Obtains all attributes for a SNP-cis gene-trans gene trio in a specified tissue
- `extract_hic()`: Extracts chromatin interaction data from .hic files for specified chromosomal regions
- `Resample_interactions()`: Performs Monte Carlo resampling to establish a null distribution of chromatin interactions
- `interaction_check()`: Primary function that combines all steps to evaluate if significant chromatin interactions exist between SNP locations and trans genes

## Workflow

1. Identify trios of interest (variant-cis gene-trans gene relationships)
2. Use `get_trio_attr()` to obtain genomic locations and other attributes
3. Extract Hi-C data for specific regions using `extract_hic()`
4. Compare observed interaction counts to null distribution using `Resample_interactions()`
5. Calculate statistical significance with multiple testing correction

## Dependencies

- strawr
- ggpubr
- plyr
- qvalue

## Example Usage

```r
# Check interactions for specific trios in lymphoblastoid cells
results <- interaction_check(
  hic.filename = "/path/to/hic/file.hic",
  trios = c(3324, 4551),
  resolution = 10000,
  search.size = 100000,
  tiss = "CellsEBVtransformedlymphocytes",
  FDR = "ADDIS",
  plot.h = TRUE
)
```

## Hi-C Data Files

The following Hi-C data files from ENCODE are used in the analysis:

```
------------------------------->>Meta_DATA_HiC_Files<<----------------------------------
File_Name            Tissue
ENCFF366ERB.hic        lung
ENCFF355OWW.hic        lymphoblastoid cells
ENCFF894RRQ.hic        fibroblast cells
ENCFF569RJM.hic        skin
```

These files correspond to the tissues analyzed in the script:
- "CellsCulturedfibroblasts" → ENCFF894RRQ.hic
- "SkinNotSunExposed" → ENCFF569RJM.hic
- "Lung" → ENCFF366ERB.hic
- "CellsEBVtransformedlymphocytes" → ENCFF355OWW.hic

## Output

The `interaction_check()` function returns:
- `summary.table1`: Results excluding NAs from p-value calculations
- `summary.table2`: Results including NAs in p-value calculations
- `data`: Raw resampled data for each trio analyzed

## Notes

- This pipeline was designed for analyzing GTEx Version 8 data
- Default search size around gene locations is 100kb
- Multiple testing correction methods include Holm-Bonferroni, Benjamini-Hochberg, and Benjamini-Yekutieli

## File Paths

The script contains hardcoded paths to data files. Update these paths before running in your environment:
- Gene metadata: `/mnt/ceph/jarredk/...`
- GTEx data: `/mnt/lfs2/mdbadsha/peer_example/...`
- Output directories: `/mnt/ceph/jarredk/HiC_Analyses/...`
