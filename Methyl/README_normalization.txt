

Outline of the normalization procedures for the BSGS Methylation and Gene Expression
data sets. 


==================================Expression===========================================

The BSGS expression data in file: 
/mnt/ceph/jarredk/Methyl/ExpressData/new_data_GE_R2_bsgs_delete.csv

The related covariates obtained from web scrape of GEO in:
/mnt/ceph/jarredk/Methyl/ExpressData/BSGS_GEO_Accession_data2.csv

original BSGS data reads in with probes in rows and subjects in columns (expression and p-values)

data was broken into 2 matrices: one of pvalues and another with expression data.

after reading in the data was transposed so probes were in columns and subjects rows

Data Dimension: 47323  1727

--------------------------------filtering-Expression-----------------------------------

the expression matrix was filtered such that probes in the expression matrix that corresponded 
to p-values with less than 10% total of p-values being significant (at alpha = 0.05) were removed.

after p-value filtering, subjects with missing covariate information were removed (4 rows total)
from both the covariate and filtered expression matrix. This aligned the covariates and expression 
matrix such that they represented the same set of subjects for analysis. 

Next, the expression data was scaled and quantile normalized for PEER. Lastly the scaled, quantile normalized 
expression data and aligned covariates were passed to PEER for further normalization.

The PEER residuals were saved as:
"/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs.RData"
or 
"/mnt/ceph/jarredk/Methyl/ExpressData/Peerdata.bsgs.CSV"

Final Data Dim:  858 24317

==================================Methylation==========================================

The BSGS expression data in file: 
"/mnt/ceph/jarredk/Methyl/new_data_M_average_delete.csv"

related covariates obtained from web scrape of GEO in:
/mnt/ceph/jarredk/Methyl/ExpressData/BSGS_GEO_Accession_data2.csv

The subject ID's key obtained from Australian Team in:
"/mnt/ceph/jarredk/Methyl/GSE53195_GSE56105_ID_equivalence.txt"

original BSGS data reads in with probes in rows and subjects in columns (methylation and p-values)

the first 3 columns of the feshly read-in data correspond to: Probe Reference ID, UCSC Gene Name, and Row Numbering
which were removed. 

data was broken into 2 matrices: one of pvalues and another with methylation data.

after reading in the data was transposed so probes were in columns and subjects rows

Data Dimension: 485577 1231

--------------------------------filtering-Methylation----------------------------------

the methylation matrix was filtered such that any probe with more than 11 subjects with missing data
or more than 5 subjects with detection p-values > 0.001 were removed. 

after filtering, the data was aligned with the covariates using the Subject ID key provided by the Australian team
and subjects with missing covariate information were removed from both the covariate and filtered methylation
matrix. 

The data was normalized by regressing out the variation due to the available covariates on each column of the 
filtered methylation matrix. First each column was transformed using the logit transformation. values on the boundary 
(beloning to exactly 0 or 1) were given adjusted by adding/subtracting a pseudocount: 0 --> 0.001, 1 --> 0.999. After 
transformation, each column was regressed on age, sex, age^2, sex*age, and sex*age^2 and the residuals were retained.

The trasnformed residuals/normalized methylation matrix is in:
/mnt/ceph/jarredk/Methyl/MethylData.RegressResids.Rdata

(note: the first row of this matrix contains the UCSC gene names for each probe) 

Final Data Dinemsion:  607 470158






