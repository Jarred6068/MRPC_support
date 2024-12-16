# MRPC_support
Code repository for Kvamme, Jarred, et al. "Types of Cis-and Trans-Gene Regulation of Expression Quantitative Trait Loci Across Human Tissues." bioRxiv (2022): 2022-01.

The following analyses were performed:

- In the folder `GSEA`, `performGSEA.R` performed gene set enrichment analysis for cis-genes and trans-genes in trios of each regulatory type.  Input data is also provided in this folder.

- In the folder `TFDatabases`, `compareWithTFDBs_TFLink.R` looked at pairs of transcription factors and target genes in mediation trios to see whether they appear in the TFLink database.  Additional R functions and intermediate output files are provided in this folder.
