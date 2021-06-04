using calcCOR():

1.) This function is a search algorithm which extracts from the Genotype, Expression, and Methylation data matrics
the SNPs which are in close proximity to each probe and calculates the correlation between them. The Output are two tables 
for each chromosome (one for methylation and one for expression) which is procedurally generated within the function. This 
table contains three columns: Probe/Gene ID, SNP ID, and Correlation Coefficient Cor(SNP, Probe/Gene)

2.) The primary inputs are named expression, genotype, and methylation matrices with probes/SNPs as columns 
and subjects/patients as rows. Next are the 3 "meta data" matrices containing the chromosomal coordinate information
for probes/Genes in all matrices. For the Expression and Methylation data, the corresponding metadata matrices should
have probes/genes in rows and 3 columns: Probe_ID, Chromosome #, and coordinate (start). The metadata for the genotype 
matrix should contain at least two columns: chromosome # and coordinate

3.) Additional inputs for the function include the parameter "bp.range"  which is a user specified tolerance for the 
maximum alotted distance between probes and SNPs for them to be considered close to one another. Lastly, the parameter
"fn" gives the path directory for the output tables. 

4.) A working example of the code can be found with the accompanying script EXcalcCOR.R 