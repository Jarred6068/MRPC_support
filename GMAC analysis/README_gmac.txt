GMAC analysis procedure:

Background: Both MRPC-LOND and MRPC-ADDIS techniques inferred a large number of trans mediated trios. The trans-mediation model has been 
previously identified but not typically thought to be the common mediation model for the relationship between an eQTL, a cis transcript and and trans 
transcript. Since this result is surprising relative to the existing literature, we sought to apply another method for inferring mediation on the 
trios we tested in GTEx. The Genomic Medaition analysis with Adaptive Confounding (GMAC) algorithm allows for adaptive selection of potential 
confounders from a large pool of known and unknown potentially confounding variables for each mediation test and. Like MRPC, GMAC also assumes the 
principle of mendelian randomization for the SNP

Goals:
1. Run GMAC on a selected number of tissues and compare the number of inferred mediation trios to that of MRPC-LOND and MRPC-ADDIS
	- each tissue was ran twice in GMAC: first with the cis probes as mediators and second with trans probes as mediators
2. determine the number of specific Cis and Trans mediated trios that all three methods had in common (if any)
3. Identify and understand differences in inferred model types between the two methods

#####################Methods##############################

File Usage and Data Assembly:

1. The files used to obtain the necessary information for each tissue applied to GMAC were:

	./data.snp.cis.trans.final."tissue.name".V8.unique.snps.RData - The trio matrix for each tissue
	./PCs.matrix."tissue.name".RData - the PC matrix for each tissue
	./"tissue.name".v8.covariates.txt - known covariates for each tissue

   - the tisses selected for analysis were the first 5 tissues with the largest sample sizes (see table 1)

2. The GMAC algorithm requires input of both known and unknown potential confounders. We used the known confounders of sex, platform, and pcr test
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
*fdr_filter which is the fdr threshold at which common child and intermediate variables are filtered left default value (fdr_filter=0.1)
*fdr which is the fdr rate used to select confounders was left at default (fdr=0.05)
*all other settings were left to their default values

5. GMAC was run twice on each set of trios:
	*Firstly run with the cis gene as the mediator
	*Second run with the trans gene as the mediator


Post-Processing Results:

1. GMAC outputs 2 sets of p-values for the mediation test (1) the p-value considering only known confounders and (2) the p-value adjusted for known and 
   selected covariates (adjusted for any/all selected PC's) we used only the latter

2. The rate of type 1 errors for the p-values retained in (1.) was controlled using an FDR of 10% and the resulting q-values were returned
   (i.e the function run.postproc() in GMACpostproc.R)

3. based on the q-values, the number of remaining significant tests were counted and the significant trios were kept and non-significant removed.

4. the significant trios retained were matched against the mediation (M1T1 and M2T2) trios in LOND and ADDIS to see how many were in common with GMAC


TABLE 1. - 

Tissue	                        Total.Num.Samples	GMAC.Total.M1T1	         GMAC.Total.M1T2	LOND.total.M1	 ADDIS.total.M1
AdiposeSubcutaneous		   11850		     453	             351	             123	      130
ArteryTibial			   11471		     372	             293	              90	      113
MuscleSkeletal			   10257		     493	             354	             104	      103
SkinSunExposed			   13045		     520	             400	             122	      123
WholeBlood			   8823		             552	             484	             131	      113



Checking Numerical Stability of GMAC Results by Simulation of the Mediation Test:

* Mediation Test refers to the regression of the trans gene with the full set of confounders by GMAC:

		Tj = b0 + b1Ci + b2Li + AXij + e

* The numerical stability of the coefficients for the SNP and Cis gene were checked by simulation for mediation models inferred by GMAC that 
  were not inferred to be mediation by ADDIS

* Each mediation test (regression) was done using all the adaptively selected and known confounders. 

* From each regression, the coefficients, b, and residual standard error, sigma, were retained and the expression of the trans gene was simulated via

		Tnew = b0 + b1Ci + b2Li + AXij + e

where 

		e ~ norm( 0, sigma)

		

* a new regression was conducted using the simulated trans gene as the response.
* this process was repeated 1000 times and the distribution of the estimated coefficients were retained.
* plots of the trios tested are in the folder /mnt/ceph/jarredk/GMACanalysis/perm_regplots/
* RESULT: instability was typically only visible in the SNP coefficient (sometimes the cis gene and sometimes both) 
	  however, the wide distribution of coefficients was not reserved to the mediation test regression but was also
	  visible in the addis regression using a much smaller set of covariates (inconclusive) 


Further Post Processing Analysis and Notes:

***We realized that GMAC reports significant mediation based only on the edge between the cis and trans gene: Ci ---> Tj
   regardless of the possible effect from SNP to cis gene: Li ---> Ci or SNP to trans gene: Li ---> Tj. As a result, many of the 
   the classifications for model types agree such as those reported M4 by ADDIS. This is because the M4 (fully connected) model
   detects an undirected edge between the cis and trans gene Ci --- Tj. 

***GMAC is unable to detect a true model 3 because model three lacks a Ci --- Tj edge and is therefore going to be classified as not significant
   according to the mediation test 

***The p-value obtained from the general regression (Tj = b0 + b1Ci + b2Li + AXij + e) can be larger (sometimes by an order of magnitude or two)
   than the p-value observed from the permuted regression: This is the typical scenario causing differences between model classes based on the 
   Ci ---> Tj edge

***The regression under ADDIS selectes confounders from top 10 PCs of the expression matrix (those having a strong association with both the cis and trans) 
   genes. Therefore, let this set of confounders be Zij such that it is a subset of Xij. The regression using adaptively selected confounders can
   occasionally reverse the decision on a weakly positive edge between two nodes. An example are trios classified M3 under the ADDIS regression 
							Tj = b0 + b1Ci + b2Li + BZij + e
   
   but classified as M0 under when using the full set of GMAC selected confounders. 



                     MO M1 M2  M3  M4 Other Total GMAC Inferred
AdiposeSubcutaneous 109 43  7 126 132     3                 420
ArteryTibial         89 25  5 139 110     0                 368
MuscleSkeletal      126 31  7 183 118     3                 468
SkinSunExposed      106 33 12 185 138     0                 474
WholeBlood          122 55 22 152 162     0                 513


***The write up for this folder is GMACwriteup.pdf and can be found in the Manuscript Directory








