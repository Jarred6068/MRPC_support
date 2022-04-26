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
-All datasets were normalized by logit transformation of the proportion of methylated sites and additional factors were adjusted for in the linear regression of the log proportion of methylation:

*Stem cell 59091:
log.p = Sex

*Stem cell 116754:
log.p = sample.type

*Healthy Breast tissue:
log.p = poly(Age,2)*poly(BMI,2)*Race

*Cancerous Breast Tissue:
log.p = Sample.Section*Sample.Plate

the residuals from each regression were retained and kept for downstream analysis

**sample sizes within each histology group for
               Grade
  =============================
  ER_status |  1* |  2*  |  3*
  ----------|-----|------|-----
  er_neg    |  2  |   8  |  60
  ----------|-----|------|-----
  er_pos    |  18 |   42 |  26

***stem cell data 59091 had the fewest genes: all data sets aligned using stem cell data 59091 as the target.
   after alignment only HESC and non-reprogrammed cells were kept for Stem Cell data.
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
-We used 1-way anova across each gene centroid where the explanatory variable is a factor with 7 levels
 "Normal" and the six combinations of grade and ER status (we excluded STEM cells from the differential analysis)
 We identify all genes with significant p-value form the omnibus F test as differentially methylated between normal and cancereous tissues

-We use a rejection threshold of alpha = 0.05 and the qvalue method to control the FDR at 5%


=================================================
---------***Tree-Construction-Methods***---------
=================================================

-starting with the matrix of identified differentially expressed genes dimension 333 obs X 1164 differentially expressed gene centroids

-we start by obtaining the group centers for each representative tissue type (STEM cell groups + healthy + an       cancer grade and ER status groups)

    dermal.fibroblast  | endothelial.precursor |      er_neg1
  |-----------------------------------------------------------------|
  |         23         |          5            |         2          |
  |-----------------------------------------------------------------|
  |       er_neg2      |       er_neg3         |      er_pos1       |
  |-----------------------------------------------------------------|
  |          8         |         60            |        18          |
  |-----------------------------------------------------------------|
  |       er_pos2      |       er_pos3         | foreskin.fibroblast|
  |-----------------------------------------------------------------|
  |          42        |         26            |         6          |
  |-----------------------------------------------------------------|
  |         HESC       |       Normal          |                    |
  |-----------------------------------------------------------------|
  |          22        |        121            |                    |
   -----------------------------------------------------------------


-Once the group centers are obtained the data is of dimension 11 tissue types X 1164 genes

-trees were constructed by calculating the distance between tissue groups across genes. We used
 10,000 bootstrappings and taking the concensus tree...Details below

-for each bootstrap:
      (1) sample 70% of the differentially methylated genes (without replacement)
      (2) calculate the distance between the tissue groups among the sampled genes
      (3) construct the histology tree
-repeate 10,000 times

-obtain the concensus tree

**we used combinations of distance metrics and tree reconstruction methods, namely
  DISTANCE: 1-R^2, euclidean, Minkowski p=1, maximum
  TREE METHODS: Neighbor Joining, Minimum Evolution Weighted Least Squares

=================================================
---------***PCA-and-Genetic-Diversity***---------
=================================================



































