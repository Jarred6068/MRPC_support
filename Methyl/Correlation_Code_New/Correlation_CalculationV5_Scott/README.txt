
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
- Example_run_v5a.R: runs a small example with just 1000 SNPs, which finishes in ~1 min.
- Run_script_v5a.R: runs on all the SNPs for all the chromosomes.
- calcCOR_v3_v5a.R: contains the functions called by the above two R scripts.

### Example_run_v5a.R ###
### Running the example and checking the output ######
Copy and paste the lines in Example_run.R into an R session.  
It will load six data files. The one for methylation will take ~20 sec to load.
It then runs the correlation calculation function, which takes 10-15 sec.
The folder /Input/ contains the input data (see Input Data section).
The output files will be stored in the folder /Example_output/ and can be compared to the sample output files M_chr1_correlations_Example.txt and E_chr1_correlations_Example.txt (See Output Data section).

### Run_script_v5a.R ###
### Running on all the genotypes #####
The only difference between running the example and running on all the genotypes is the genotype data files.
The user will need to specify the filenames for the genotype data and genotype metadata in Run_script.R.
This script runs on each chromosome individually.
The folder /Input/ contains the input data (see Input Data section).
The output files will be stored in the empty folder /Output/.
Can be parallelized for efficiency (see description in Run_script.R).

### calcCOR_v5a.R ###
Example_script.R and Run_script.R call functions in this script.


###############################################
# Input data
###############################################
- Data files in /Input/:
* gene expression data (Express.BM.aligned.Ilmn.aligned.Rdata): a matrix of 606 individuals in rows and 20,578 genes/probes in columns.
* gene expression metadata (Express.BM.aligned.Ilmn.aligned.MetaData.Rdata): a meta data matrix of 20,578 probes in rows and 5 columns: ID, chr, coordinate, Gene.name].
* methylation data (Mdata.Rdata): a matrix of 606 individuals in rows and 470,158 probes in columns.
* methylation metadata (meta_M.Rdata): a meta data matrix of 470,158 probes and 3 columns: Probe_ID, chr, coordinate.
* genotype data: a data matrix with individuals in the rows, and SNPs in the columns.
* genotype metadata: a data matrix with SNPs in the rows, and 5 columns: SNP_ID, Chromosome Number, and Chromosome Coordinate, Other Allele, Counted Allele.

- Format of the example genotype data file (/Input/Example_Genotype_data.txt):

           SNP_ID     bsgs_1 bsgs_2 bsgs_3 bsgs_4     ...
  SNP_130647807_chr4     2      0      1      1       ... 
  SNP_50536821_chrY      2      0      2      1       ...
         .               .      .      .      .       ...
         .               .      .      .      .       ...
         .               .      .      .      .       ...
  SNP_193452519_chr2     2      2      0      1       ...     



- Format of the example genotype metadata file (/Input/Example_Genotype_Metadata.txt):

SNP_ID	chr	coordinate	OtherAllele	CountedAllele
SNP_130647807_chr4	4	130647807	C	A
SNP_50536821_chrY	Y	50536821	A	G
	.		.	    .		.	.
	.		.	    .		.       .
	.		.	    .		.	.
SNP_193452519_chr2	2	193452519	C	T


###################################
# Output data
###################################
- Running Run_script.R will produce output files and store them in the empty folder /Output/.
- Running Example_run.R will produce two output files and store them in /Example_Output/.
- The sample output files in /Example_Output/ are M_chr1_correlations_Example.txt and E_chr1_correlations_Example.txt.
- The format of E_chr1_correlations_Example.txt is shown below:

Express_Gene_ID	SNP_ID	chr	coordinate	Cor	CountedAllele	OtherAllele
ILMN_1651229	NA	NA	NA	NA	NA	NA
ILMN_1651278	SNP_37600171_chr1	0.077811	1	37600171	A	C
ILMN_1651385	NA	NA	NA	NA	NA	NA
ILMN_1651492	NA	NA	NA	NA	NA	NA
ILMN_1651554	NA	NA	NA	NA	NA	NA
ILMN_1651800	NA	NA	NA	NA	NA	NA
ILMN_1651828	NA	NA	NA	NA	NA	NA
ILMN_1651872	NA	NA	NA	NA	NA	NA
ILMN_1651949	NA	NA	NA	NA	NA	NA
ILMN_1652205	NA	NA	NA	NA	NA	NA
ILMN_1652459	NA	NA	NA	NA	NA	NA
ILMN_1652677	NA	NA	NA	NA	NA	NA
ILMN_1652758	NA	NA	NA	NA	NA	NA
ILMN_1652929	SNP_151704009_chr1	0.028359	1	151704009	T	A
ILMN_1653146	NA	NA	NA	NA	NA	NA
ILMN_1653367	NA	NA	NA	NA	NA	NA
ILMN_1653429	NA	NA	NA	NA	NA	NA
ILMN_1653466	NA	NA	NA	NA	NA	NA
ILMN_1653494	NA	NA	NA	NA	NA	NA
ILMN_1653496	NA	NA	NA	NA	NA	NA
ILMN_1653504	NA	NA	NA	NA	NA	NA
ILMN_1653514	NA	NA	NA	NA	NA	NA
ILMN_1653524	NA	NA	NA	NA	NA	NA
ILMN_1653618	SNP_77906363_chr1	-0.004054	1	77906363	A	G

- The format of M_chr1_correlations_Example.txt is shown below:

Methyl_Probe_ID	SNP_ID	Cor	chr	coordinate	CountedAllele	OtherAllele
cg00000165	NA	NA	NA	NA	NA	NA
cg00000363	NA	NA	NA	NA	NA	NA
cg00000957	SNP_6244975_chr1	-0.015595	1	6244975	G	C
cg00001349	NA	NA	NA	NA	NA	NA
cg00001364	NA	NA	NA	NA	NA	NA
cg00001446	NA	NA	NA	NA	NA	NA
cg00001534	NA	NA	NA	NA	NA	NA
cg00001583	NA	NA	NA	NA	NA	NA
cg00001593	NA	NA	NA	NA	NA	NA
cg00002028	NA	NA	NA	NA	NA	NA
cg00002593	NA	NA	NA	NA	NA	NA
cg00002646	NA	NA	NA	NA	NA	NA
cg00002719	NA	NA	NA	NA	NA	NA
cg00002808	NA	NA	NA	NA	NA	NA
cg00002810	NA	NA	NA	NA	NA	NA
cg00002837	NA	NA	NA	NA	NA	NA
cg00003187	NA	NA	NA	NA	NA	NA
cg00003202	SNP_151704009_chr1	0.03503	1	151704009	T	A
cg00003287	NA	NA	NA	NA	NA	NA

                 













