using calcCOR():

1.) This function is a search algorithm which extracts the SNPs from the Genotype data which are in close proximity to each 
probe in the Expression, and Methylation data sets and calculates the correlation between them. The Output are two tables 
for each chromosome (one for methylation and one for expression) which are procedurally generated within the function. This 
table contains three columns: Probe/Gene ID, SNP ID, and Correlation Coefficient Cor(SNP, Probe/Gene)

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
matrix should contain at two columns: chromosome # and coordinate

3.) Additional inputs for the function include the parameter "bp.range"  which is a user specified tolerance for the 
maximum alotted distance between probes and SNPs for them to be considered close to one another. Lastly, the parameter
"fn" gives the path directory for the output tables. 

4.) A working example of the code can be found with the accompanying script EXcalcCOR.R 