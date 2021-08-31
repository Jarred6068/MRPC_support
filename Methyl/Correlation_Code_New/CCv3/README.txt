
###################
Table of Contents
###################

- R scripts
- Input Data
- Output Data


#################################################
# R scripts
#################################################

### Overview ######
This folder contains three R scripts:
- Example_run.R: runs a small example with just 100 SNPs, which finishes in ~1 min.
- Run_script.R: runs on all the SNPs for all the chromosomes.
- calcCOR_v3.R: contains the functions called by the above two R scripts.

### Example_run.R ###
### Running the example and checking the output ######
Copy and paste the lines in Example_run.R into an R session.  
It will load six data files. The one for methylation will take ~20 sec to load.
It then runs the correlation calculation function, which takes 10-15 sec.
The folder /Input/ contains the input data (see Input Data section).
The output files will be stored in the folder /Example_output/ and can be compared to the sample output files M_chr1_correlations_Example.txt and E_chr1_correlations_Example.txt (See Output Data section).

### Run_script.R ###
### Running on all the genotypes #####
The only difference between running the example and running on all the genotypes is the genotype data files.
The user will need to specify the filenames for the genotype data and genotype metadata in Run_script.R.
This script runs on each chromosome individually.
The folder /Input/ contains the input data (see Input Data section).
The output files will be stored in the empty folder /Output/.
Can be parallelized for efficiency (see description in Run_script.R).

### calcCOR_v3.R ###
Example_script.R and Run_script.R call functions in this script.


###############################################
# Input data
###############################################
- Data files in /Input/:
* gene expression data (Expression_data_BM_aligned.Rdata): a matrix of 606 individuals in rows and 20,578 genes/probes in columns.
* gene expression metadata (meta_E.Rdata): a meta data matrix of 20,578 probes in rows and 5 columns: gene_ID, chr, gene.start, gene.end,   ILMN_ID].
* methylation data (Mdata.Rdata): a matrix of 606 individuals in rows and 363,516 probes in columns.
* methylation metadata (meta_M.Rdata): a meta data matrix of 363,516 probes and 3 columns: Probe_ID, chr, coordinate.
* genotype data: a data matrix with individuals in the rows, and SNPs in the columns.
* genotype metadata: a data matrix with SNPs in the rows, and 3 columns: SNP_ID, Chromosome #, and Chromosome Coordinate.

- Format of the example genotype data file (/Input/Example_Genotype_data.txt):

           SNP_ID     bsgs_1 bsgs_2 bsgs_3 bsgs_4     ...
  SNP_1942635_chr_6      2      0      1      1       ... 
  SNP_2610867_chr_15     2      0      2      1       ...
         .               .      .      .      .       ...
         .               .      .      .      .       ...
         .               .      .      .      .       ...
  SNP_2237944_chr_12     2      2      0      1       ...     



- Format of the example genotype metadata file (/Input/Example_Genotype_Metadata.txt):

SNP_ID                  Chr    Coordinate
SNP_1942635_chr_6        6     24423331
SNP_2610867_chr_15      15     82155698
       .                 .        .
       .                 .        .
       .                 .        .
SNP_2237944_chr_12      12     62931549


###################################
# Output data
###################################
- Running Run_script.R will produce output files and store them in the empty folder /Output/.
- Running Example_run.R will produce two output files and store them in /Example_Output/.
- The sample output files in /Example_Output/ are M_chr1_correlations_Example.txt and E_chr1_correlations_Example.txt.
- The format of E_chr1_correlations_Example.txt is shown below:

Express_Gene_ID	            SNP_ID	             Cor
ILMN_1672389:HLA-DQB1	SNP_428936_chr_1	-0.037259
ILMN_1683980:LOC644482	SNP_1287122_chr_1	-0.017493
ILMN_1683990:LOC645914	SNP_1287122_chr_1	-0.054687

- The format of M_chr1_correlations_Example.txt is shown below:

Methyl_Probe_ID	      SNP_ID	           Cor
cg03374965	SNP_428936_chr_1	-0.02847
cg03374976	SNP_428936_chr_1	 0.019684
cg03375002	SNP_428936_chr_1	-0.031407

                 













