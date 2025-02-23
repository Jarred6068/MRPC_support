##########################################
# Perform gene set enrichment analysis
##########################################
library(gprofiler2)

setwd("~/Documents/CisTrans/Data")
all.trios <- read.table(gzfile("all_trios_inferred_models.txt.gz"), header = TRUE, sep = "\t")
dim (all.trios)
# 323423     17

all.trios.filtered <- all.trios[which (all.trios$cis.gene.type=="protein_coding" | all.trios$cis.gene.type=="lncRNA"),]
all.trios.filtered <- all.trios.filtered[which (all.trios.filtered$trans.gene.type=="protein_coding" | all.trios.filtered$trans.gene.type=="lncRNA"),]
dim (all.trios.filtered)
#[1] 108446     17

all.trios <- all.trios.filtered

# M0
id.m0 <- which (all.trios$lond.model=="M0")
genes.cis.m0 <- unique (all.trios$cis.gene.name[id.m0])
genes.trans.m0 <- unique (all.trios$trans.gene.name[id.m0])

genes.cis.m0.go.all <- gost(query = genes.cis.m0, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m0.go.all$result)
# 3796   15
genes.cis.m0.go.all$result$query <- rep ("m0.cis", nrow(genes.cis.m0.go.all$result))
gostplot(genes.cis.m0.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m0.go.all <- gost(query = genes.trans.m0, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.trans.m0.go.all$result)
# 486  15
genes.trans.m0.go.all$result$query <- rep ("m0.trans", nrow(genes.trans.m0.go.all$result))
gostplot(genes.trans.m0.go.all, capped = TRUE, interactive = TRUE)


# M1
id.m1 <- which (all.trios$lond.model=="M1")
genes.cis.m1 <- unique (all.trios$cis.gene.name[id.m1])
genes.trans.m1 <- unique (all.trios$trans.gene.name[id.m1])

genes.cis.m1.go.all <- gost(query = genes.cis.m1, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m1.go.all$result)
# 107  15
genes.cis.m1.go.all$result$query <- rep ("m1.cis", nrow(genes.cis.m1.go.all$result))
gostplot(genes.cis.m1.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m1.go.all <- gost(query = genes.trans.m1, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.trans.m1.go.all$result)
# 27 15
genes.trans.m1.go.all$result$query <- rep ("m1.trans", nrow(genes.trans.m1.go.all$result))
gostplot(genes.trans.m1.go.all, capped = TRUE, interactive = TRUE)


# M2
id.m2 <- which (all.trios$lond.model=="M2")
genes.cis.m2 <- unique (all.trios$cis.gene.name[id.m2])
genes.trans.m2 <- unique (all.trios$trans.gene.name[id.m2])

genes.cis.m2.go.all <- gost(query = genes.cis.m2, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m2.go.all$result)
# 8  15
genes.cis.m2.go.all$result$query <- rep ("m2.cis", nrow(genes.cis.m2.go.all$result))
gostplot(genes.cis.m2.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m2.go.all <- gost(query = genes.trans.m2, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.trans.m2.go.all$result)
# no result
#genes.trans.m2.go.all$result$query <- rep ("m2.trans", nrow(genes.trans.m2.go.all$result))
#gostplot(genes.trans.m2.go.all, capped = TRUE, interactive = TRUE)


# M3
id.m3 <- which (all.trios$lond.model=="M3")
genes.cis.m3 <- unique (all.trios$cis.gene.name[id.m3])
genes.trans.m3 <- unique (all.trios$trans.gene.name[id.m3])

genes.cis.m3.go.all <- gost(query = genes.cis.m3, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m3.go.all$result)
# 1852   15
genes.cis.m3.go.all$result$query <- rep ("m3.cis", nrow(genes.cis.m3.go.all$result))
gostplot(genes.cis.m3.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m3.go.all <- gost(query = genes.trans.m3, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.trans.m3.go.all$result)
# 663  15
genes.trans.m3.go.all$result$query <- rep ("m3.trans", nrow(genes.trans.m3.go.all$result))
gostplot(genes.trans.m3.go.all, capped = TRUE, interactive = TRUE)


# M4
id.m4 <- which (all.trios$lond.model=="M4")
genes.cis.m4 <- unique (all.trios$cis.gene.name[id.m4])
genes.trans.m4 <- unique (all.trios$trans.gene.name[id.m4])

genes.cis.m4.go.all <- gost(query = genes.cis.m4, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m4.go.all$result)
# 7
genes.cis.m4.go.all$result$query <- rep ("m4.cis", nrow(genes.cis.m4.go.all$result))
gostplot(genes.cis.m4.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m4.go.all <- gost(query = genes.trans.m4, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.trans.m4.go.all$result)
# 4
gostplot(genes.trans.m4.go.all, capped = TRUE, interactive = TRUE)


# Plot_genes_enrichment_all.pdf
gostplot(genes.cis.m0.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.trans.m0.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.cis.m1.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.trans.m1.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.cis.m2.go.all, capped = TRUE, interactive = FALSE)
#gostplot(genes.trans.m2.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.cis.m3.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.trans.m3.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.cis.m4.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.trans.m4.go.all, capped = TRUE, interactive = FALSE)

# count the total number of genes for each enrichment analysis
length (genes.cis.m0)
# 15703
length (genes.trans.m0)
# 12866
length (genes.cis.m1)
# 1224
length (genes.trans.m1)
# 1011
length (genes.cis.m2)
# 673
length (genes.trans.m2)
# 474
length (genes.cis.m3)
# 14639
length (genes.trans.m3)
# 15929
length (genes.cis.m4)
# 845
length (genes.trans.m4)
# 657

# Merge all the enrichments 
genes.go.all <- rbind (genes.cis.m0.go.all$result, genes.trans.m0.go.all$result, genes.cis.m1.go.all$result, 
                       genes.trans.m1.go.all$result, genes.cis.m2.go.all$result, genes.trans.m2.go.all$result,
                       genes.cis.m3.go.all$result, genes.trans.m3.go.all$result, genes.cis.m4.go.all$result)
dim (genes.go.all)
# 6946   15

# The 'parents' column is a list.  Need to convert it to a vector
genes.go.all$parents <- sapply(genes.go.all$parents, function(x) paste(x, collapse = ","))

# Write to output
library(openxlsx)
write.xlsx(genes.go.all, "mediation_genes_GSEA_all.xlsx")
