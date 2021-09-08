GMAC analysis procedure:

Background: Both MRPC-LOND and MRPC-ADDIS techniques inferred a large number of trans mediated trios. The trans-mediation model has been 
previously identified but not typically thought to be the common mediation model for the relationship between an eQTL, a cis transcript and and trans 
transcript. Since this result is surprising relative to the existing literature, we sought to apply another method for inferring mediation on the 
trios we tested in GTEx. The Genomic Medaition analysis with Adaptive Confounding (GMAC) algorithm allows for adaptive selection of potential 
confounders from a large pool of known and unknown potentially confounding variables for each mediation test and, like MRPC, GMAC also assumes the 
principle of mendelian randomization

Goals:
1. Run GMAC on a selected number of tissues and compare the number of inferred mediation trios to that of MRPC-LOND and MRPC-ADDIS
	- each tissue was ran twice in GMAC: first with the cis probes as mediators and second with trans probes as mediators
2. determine the number of specific Cis and Trans mediated trios that all three methods had in common (if any)

#####################Methods##############################

File Usage and Data Assembly:

1. The files used to obtain the necessary information for each tissue applied to GMAC were:

	./data.snp.cis.trans.final."tissue.name".V8.unique.snps.RData - The trio matrix for each tissue
	./PCs.matrix."tissue.name".RData - the PC matrix for each tissue
	./"tissue.name".v8.covariates.txt - known covariates for each tissue

   - the tisses selected for analysis were the first 5 tissues with the largest sample sizes (see table 1)

2. The GMAC algorithm requires input of both known and unknown potential confounders. We used the known confounders of sex, age, and pcr test
   as the known set of confounders and the full set of principle components obtained from the PCA on the expression matrix for each tissue as the 
   pool of unknown potential confounders

3. A function was used to assemble the above files into a list. 
	*The expression and SNP matrices were extracted from the trio matrix. 
	*the SNP matrix consisted of the columns pertaining to the unique SNPs in the trio matrix from each tissue. 
	*The expression matrix was retained as the cis and trans expression probes from the trio matrix after removing all SNP columns. 

   	*The trio index map was constructed by converting the column names of the trio matrix into an (n X 3) matrix where each 
	  consists of a trio with the SNPs in the first column and cis and trans expression probes in the remaining two columns respectively. The index
	  map was created by matching the name in each row of this (n X 3) array to the colnames of the unique SNP and expression matrices

4. Many of the SNP's contained at least one missing value: 
	*Imputation using Multiple Correspondance Analysis (MCA - missMDA package in r) was applied to the unique SNP matrix prior to use in GMAC
	*Becasue the Expression data was PEER normalized it contained no missing values

GMAC Settings:

*GMAC was run on all 5 tissues using first 50 permutations (nperm = 50) and later 500 (as the algorithm is slow for such a large data set). 
*The return p-values were the nomial p-values (nominal.p = TRUE)
*fdr_filter which is the fdr threshold at which common child and intermediate variables are filtered (left default at fdr_filter=0.1)
*fdr which is the fdr rate used to select confounders was left at default (fdr=0.05)
*all other settings were left to their default values

Post-Processing Results:

1. GMAC outputs 2 sets of p-values for the mediation test (1) the p-value considering only known confounders and (2) the p-value adjusted for known and 
   selected covariates (adjusted for any/all selected PC's) we used only the latter

2. The rate of type 1 errors for the p-values retained in (1.) was controlled using an FDR of 10% and the resulting q-values were returned

3. based on the q-values, the number of remaining significant tests were counted and the significant trios were kept and non-significant removed.

4. the significant trios retained were matched against the mediation trios in LOND and ADDIS to see how many were in common with GMAC


TABLE 1. - 

Tissue	                        Total.Num.Samples	GMAC.Total.M1T1	         GMAC.Total.M1T2	LOND.total.M1	 ADDIS.total.M1
AdiposeSubcutaneous		   11850		     453	             351	             123	      130
ArteryTibial			   11471		     372	             293	              90	      113
MuscleSkeletal			   10257		     493	             354	             104	      103
SkinSunExposed			   13045		     520	             400	             122	      123
WholeBlood			   8823		             552	             484	             131	      113






















