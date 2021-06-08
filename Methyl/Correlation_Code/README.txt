using calcCOR_v2_final:

1.) This set of functions works as a search algorithm which extracts the SNPs from the Genotype data which are in close 
proximity to each probe in the Expression and Methylation data sets and calculates the correlation between them. The 
Output are two tables for each chromosome (one for methylation and one for expression) which are procedurally generated 
within the function. This table contains three columns: Probe/Gene ID, SNP ID, and Correlation Coefficient Cor(SNP, Probe/Gene)

Example Expression output for chromosome 3

Express_Gene_ID	          SNP_ID	     	   Cor
ILMN_1651329:LPP	SNP_402_chr_3	 	 0.039774
ILMN_1651329:LPP	SNP_9178_chr_3		-0.017883
ILMN_1651329:LPP	SNP_13067_chr_3	 	 0.053581
ILMN_1651329:LPP	SNP_13729_chr_3	 	-0.015380
ILMN_1651329:LPP	SNP_13902_chr_3		-0.023204
ILMN_1651329:LPP	SNP_15250_chr_3	 	 0.033120
ILMN_1651329:LPP	SNP_15864_chr_3	 	 0.024005

Example Methylation output for chromosome 3

Methyl_Probe_ID   SNP_ID                    Cor
cg00013409      SNP_479_chr_3          -0.0849196
cg00013409      SNP_1783_chr_3  	0.0150559
cg00013409      SNP_2187_chr_3  	0.0159960
cg00013409      SNP_5083_chr_3  	0.0173898
cg00013409      SNP_9178_chr_3  	0.0373560
cg00013409      SNP_12153_chr_3        -0.0417238
cg00013409      SNP_13902_chr_3 	0.0089846
cg00013409      SNP_15250_chr_3 	0.0575382
cg00013409      SNP_17563_chr_3 	0.0408202


2.) The primary inputs are named expression, genotype, and methylation matrices with probes/SNPs as columns 
and subjects/patients as rows. Next are the 3 "meta data" matrices containing the chromosomal coordinate information
for probes/Genes in all matrices. For the Expression and Methylation data, the corresponding metadata matrices should
have probes/genes in rows and 3 columns: Probe_ID, Chromosome #, and coordinate (start). The metadata for the genotype 
matrix should contain two columns: chromosome # and coordinate (in that order). **It is important that each row of the genotype metadata
matrix align with each column of the genotype data matrix** such as:

(a) Ex genotype metadata

              sim.chr sim.chrpos
SNP_1_chr_16       16   76868343
SNP_2_chr_8         8   97060256
SNP_3_chr_14       14   77303865
SNP_4_chr_19       19   22429168
SNP_5_chr_6         6  102339318
SNP_6_chr_14       14   88219023
SNP_7_chr_X         X  100268374
SNP_8_chr_13       13   12170231
SNP_9_chr_4         4  148975113
SNP_10_chr_18      18   67025553
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

3.) Additional inputs for the function include the parameter "bp.range"  which is a user specified tolerance for the 
maximum alotted distance between probes and SNPs for them to be considered "close" to one another. Lastly, the parameter
"fn" gives the path directory for the output tables. 

4.) The accompanying file Run_script.R is a simple script used to run the algorithm and is set to run the search on each 
chromosome individually. As some of the generated tables may be quite large this script should help the user if they wish to 
parallelize the operations. Note that the only input(s) the user must make is to enter the correct path for the variable(s)
"load.location = /path/ " and "save.location = /path/ " which idenitfies the designated folder(s) for the input and output files.
Additionally, the variable "bp.set.range=500000" is pre-set to allow a 500,000 bp tolerane on either side of the probe/gene start
location, but it may be adjusted down as needed for memory purposes. 

5.) A short example script EXrun.R is included as well as example .txt files for the genotype data and associated genotype meta data.
The example runs quickly and therefore will allow the user to familiarize themselves with the code before use. There are also two 
.txt files showing the output generated from EXrun.R : "M_chr1_correlations.txt" which is the output file for the methylation data 
and "E_chr1_correlations.txt" for the expression data in the format described in (1.) above











