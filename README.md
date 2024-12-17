# MRPC_support
Code repository for Kvamme, Jarred, et al. "Types of Cis-and Trans-Gene Regulation of Expression Quantitative Trait Loci Across Human Tissues." bioRxiv (2022): 2022-01.

The following analyses were performed:

- In the folder `./ADDIS_CisTrans`, the MRPC-ADDIS and LOND methods were applied to trios from 48 tissues in the GTEx consortium. 
    - The subdirectory `./Data` contains intermediate data files and results see `README_Data.txt`. 
    - The subdirectory `./Tables` contains intermediate results
    - The subdirectory `./Documents` contains two .Rmd files which both contain code and description of the gene type enrichment analysis of of MRPC-ADDIS and MRPC-LOND results for trios identified as mediators.  
    - The script `ADDIS_rerun.R` performs the trio analysis with MRPC-ADDIS for all 48 GTEx tissues
    - The script `gen_all_mediation_files.R` creates results tables for trios inferred cis or trans mediation by MRPC-ADDIS
    - the scripts `ADDIS_Post_Analysis_processing.R` and `AL_genetabV2.R` contain supporting functions used in analyzing and producing results.
    - Additional R functions and intermediate output files are provided in this folder.

- In the folder `./Analysis with GMAC` the GMAC function from Yang et al, 2017 "Identifying cis-mediators for trans-eQTLs across many
human tissues using genomic mediation analysis" was applied to trios from the top 5 five tissues by sample size from GTEx and compared with results from MRPC. 
    - The subdirectory `run_GMAC_scripts` contains scripts used to run GMAC on the trios in each tissue. See `README_gmac.txt` 
    - The sudirectory `/sim_scripts` contains the scripts and functions for comparing GMAC and MRPC and simulated mediation trios. See `GMAC_Simulations_README.txt` 
    - The subdirectory `./Lab Notes` contains intermediate progress reports for the analses of trios with GMAC
    - The script `GMACpostproc.R` contains a set of helper functions used across post inference analyses and generating results. 
    - The script `GMACanalysis.R` contains wrapper functions used for running GMAC. 
    - Additional R scripts and intermediate output files are provided in this folder

- In the folder `./GSEA`, `performGSEA.R` performed gene set enrichment analysis for cis-genes and trans-genes in trios of each regulatory type.  Input data is also provided in this folder.

- In the folder `./TFDatabases`, `compareWithTFDBs_TFLink.R` looked at pairs of transcription factors and target genes in mediation trios to see whether they appear in the TFLink database.  Additional R functions and intermediate output files are provided in this folder.

- In the folder `./HiC`, the script `Run_HiC_DA.R` investigated enrichment of HiC reads among cis and trans mediation trios. The data for this analysis were obtained from Encode - See `DataFilesMeta`.txt
    - The script `HiC_search.R` contains R functions which support the HiC enrichment analysis.
    - The script `Run_HiC_DA.R` performs the analysis checking for interaction enrichment among trios inferred as mediation by MRPC-ADDIS 
    - The subdirectory `./Results` contains intermediate results from this analysis

- The folder `./Manuscript`, contains R scripts to create tables and figures of results, supplementary figures, and other intermediate outputs. This folder also contains the manuscript document, intermediate reports, other materials needed to produce manuscript results

- The folder `./SuppTables` contains the tables included in the supplementary material of the manuscript.  
