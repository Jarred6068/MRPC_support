#################################################
# R scripts
#################################################
### Overview ######
This folder contains three R scripts:
- Example_run.R: to run a small example with just 100 SNPs, which finishes in ~1 min
- Run_script.R: runs on all the SNPs for all the chromosomes
- calcCOR_v2.R: contains the functions called by the above two R scripts

### Running the example and checking the output ######
Copy and paste the lines in Example_run.R into an R session.  
It will load six data files.  The one for methylation will take ~20 sec to load.
It then runs the correlation calculation function, which takes 10-15 sec.
It generates two output text files based on the chromosome number i.e M_chr1_correlations.txt and E_chr1_correlations.txt.
The expected output files from this test run are located under Example_output/.

### Running on all the genotypes #####
The only difference between the example and running on all the genotypes is the genotype data file.
A fake genotype data file and its fake metadata file are included to illustrate the file format and usage (please see explanation on the data files below).

### Functions in calcCOR_v2.R #####

** A short example script Example_run.R is included as well as example input files for the genotype data and associated genotype meta data. 
The example runs quickly and therefore will allow the user to familiarize themselves with the code before use. There are also two  
.txt files showing the output generated from Example_run.R : "M_chr1_correlations.txt" which is the output file for the methylation data  
and "E_chr1_correlations.txt" for the expression data

** The accompanying file Run_script.R is a simple script used to run the algorithm and is set to run the search on each  
chromosome individually. As some of the generated tables may be quite large this script should help the user if they wish to  
parallelize the operations. Note that there are two subfolders within "/Correlation_Code/": a subfolder labeled "Input" contains 
the Methylation, Expression, Metadata, and Example input files. The user supplied Genotype data, once formatted, should be saved 
to this folder. The other included folder, "Output", is an empty location for the output tables to saved. Note that the only input(s) the user must make 
is to enter the correct path for the variable(s) "input.location = "/path/Correlation_Code/" and "output.location = "/path/Correlation_Code/"
which identifies the designated "Input" and "Output" subfolders. Additionally, the variable "bp.set.range=500000" is pre-set to allow a 
500,000 bp tolerance on either side of the probe/gene start location, but it may be adjusted down as needed for memory purposes.  




using calcCOR_v2.R:

** This set of functions works as a search algorithm which extracts the SNPs from the Genotype data which are in close 
proximity to each probe in the Expression and Methylation data sets and calculates the correlation between them. The 
Output are two tables for each chromosome (one for methylation and one for expression) which are procedurally generated 
within the function. This table contains three columns: Probe/Gene ID, SNP ID, and Correlation Coefficient Cor(SNP, Probe/Gene)


** Additional inputs for the function include the parameter "bp.range"  which is a user specified tolerance for the 
maximum alotted distance between probes and SNPs for them to be considered "close" to one another. Lastly, the parameter
"fn" gives the path directory for the output tables. 





###############################################
# Input data
###############################################
- Data files needed as input:
* genotype data: a data matrix with individuals in the rows, and the SNPs in the columns.
* genotype metadata: a data matrix with SNPs in the rows, and three columns: xx and xxx.
* gene expression data (filename): a data matrix ...
* gene expression metadata (filename): a data matrix ...
* methylation data (filename): ...
* methylation metadata (filename): ...

- Format of the example genotype data file (/Correlation_Code/Input/Example_Genotype_data.txt):




- Format of the example genotype metadata file (/Correlation_Code/Input/Example_Genotype_Metadata.txt):





###################################
# Output data
###################################
- Two output files are generated: xxx1 and xxx2.
- The format of xxx1 is shown below:
Example Expression output for chromosome 3

Express_Gene_ID	          SNP_ID	     	   Cor
ILMN_1651329:LPP	SNP_402_chr_3	 	 0.039774
ILMN_1651329:LPP	SNP_9178_chr_3		-0.017883
ILMN_1651329:LPP	SNP_13067_chr_3	 	 0.053581

- The format of xxx2 is shown below:
Example Methylation output for chromosome 3

Methyl_Probe_ID   SNP_ID                    Cor
cg00013409      SNP_479_chr_3          -0.0849196
cg00013409      SNP_1783_chr_3  	0.0150559
cg00013409      SNP_2187_chr_3  	0.0159960


** The primary inputs are named expression, genotype, and methylation matrices with probes/SNPs as columns 
and subjects/patients as rows. Next are the 3 "meta data" matrices containing the chromosomal coordinate information
for probes/Genes in all matrices. For the Expression and Methylation data, the corresponding metadata matrices should
have probes/genes in rows and 3 columns: Probe_ID, Chromosome #, and coordinate (start). The metadata for the genotype 
matrix should contain two columns: chromosome # and coordinate (in that order). **It is important that each row of the genotype metadata
matrix align with each column of the genotype data matrix** such as:

(a) Ex genotype metadata

              sim.chr sim.chrpos
SNP_1_chr_16       16   76868343
SNP_2_chr_8         8   97060256
     .              .       .
     .              .       .
     .              .       .
SNP_10000_chr_7     7   23543299


(b) Ex genotype data

        SNP_1_chr_16 SNP_2_chr_8 SNP_3_chr_14 SNP_4_chr_19 SNP_5_chr_6	   ...    SNP_10000_chr_7
bsgs_1             2           2            2            1           2	   ...                  2
bsgs_2             2           0            1            0           2	   ...                  1
bsgs_3             2           2            0            1           0	   ...                  0
bsgs_4             0           2            2            0           0     ...                  0
bsgs_5             1           2            2            2           1     ...                  1
bsgs_6             0           0            2            1           1     ...                  1
bsgs_7             1           2            0            0           2     ...                  2
bsgs_8             2           2            1            1           1     ...                  2
bsgs_9             0           2            1            0           0     ...                  0
bsgs_10            1           0            0            0           1     ...                  0
                 

Note that the rows in (a) match the columns in (b). 













