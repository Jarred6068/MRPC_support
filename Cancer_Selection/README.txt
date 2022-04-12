4 data methylation data sets
====================================================================================================================
Name         |     Count    |  Sample Size    |     GSE series #      |       Tissue Labels File:                  |
-------------|--------------|-----------------|-----------------------|--------------------------------------------|
STEM         |      2       |     38          |      116754           |       *./GSE116754_meta_info.csv           |
             |              |     195         |      59091            |       *./GSE59091_meta_info.csv            |
-------------|--------------|-----------------|-----------------------|--------------------------------------------|
Healthy      |      1       |     121         |      101961           |       *./GSE101961_meta_info.csv           |
-------------|--------------|-----------------|-----------------------|--------------------------------------------|
Breast Cancer|      1       |     188         |      75067            |       *./GSE75067_sample_annotations.txt   |
-------------|--------------|-----------------|-----------------------|--------------------------------------------|


=================================================
---------***Normalization--Process***------------
=================================================
-All datasets were normalized by logit transformation of the proportion of methylated sites and additional factors were 
adjusted for in the linear regression of the log proportion of methylation:

*Stem cell 59091:
log.p = Sex 

*Stem cell 116754:
log.p = sample.type

*Healthy Breast tissue:
log.p = Race + Age + BMI + Age^2 + BMI^2 + Race*Age + Race*BMI + Race*Age^2 + Race*BMI^2

*Cancerous Breast Tissue:
log.p = Sample.Section + Sample.Plate + Sample.Section*Sample.Plate

the residuals from each regression were retained and kept for the analysis

**sample sizes within each histology for breast cancer sample sizes:        
               Grade
  =============================
  ER_status |  1  |   2  |  3
  ----------|-----|------|-----
  er_neg    |  2  |   8  |  60
  ----------|-----|------|-----
  er_pos    |  18 |   42 |  26

***stem cell data 59091 had the fewest genes: all data sets aligned using stem cell data 59091 as target.
   after alignment only HESC and non-reprogrammed cells were kept for SC data.
***rows that did not contain a gene name were omitted from the final data set using the HM450k manifest as reference

final dimension: 365,860 X 333
SC59091: 365,860 X 47
SC116754: 365,860 X 9
Healthy: 365,860 X 121
BC: 365,860 X 155

=================================================
---------***Calculating--Centroids***------------
=================================================
-Illumina 450k human methylation assay contains 485K probes and the manifest that comes standard 
 with the assay contains gene names for 365K of theses probes (119,652 with missing UCSC Ref_names)
-We filter out all probes with missing UCSC Ref_names 

-we calculated the gene centroids for the remaining probes by averaging over all probes into the 
 same gene (accounting for aliases) gene centroids data dimension:
 21,244 genes X 333 samples

=================================================
-------***Detecting-Differential-Genes***--------
=================================================
-We used 1-way anova across each gene centroid where the explanatory variable is a factor with two levels
 "Normal" or "Cancerous" indicating the state of the tissue for each individual. We retain genes 
 with significant F from this test as differentially methylated between normal and cancereous tissues

-We use a rejection threshold of alpha = 0.05 and the qvalue method to control the FDR at 5%


