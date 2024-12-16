##########################################
# Perform gene set enrichment analysis
##########################################
library(gprofiler2)

all.trios <- read.table(gzfile("all_trios_inferred_models.txt.gz"), header = TRUE, sep = "\t")
dim (all.trios)
# 323423     17

# M0
id.m0 <- which (all.trios$lond.model=="M0")
genes.cis.m0 <- unique (all.trios$cis.gene.name[id.m0])
genes.trans.m0 <- unique (all.trios$trans.gene.name[id.m0])

genes.cis.m0.go.all <- gost(query = genes.cis.m0, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m0.go.all$result)
# 4715   15
genes.cis.m0.go.all$result$query <- rep ("m0.cis", nrow(genes.cis.m0.go.all$result))
gostplot(genes.cis.m0.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m0.go.all <- gost(query = genes.trans.m0, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.trans.m0.go.all$result)
# 236  15
genes.trans.m0.go.all$result$query <- rep ("m0.trans", nrow(genes.trans.m0.go.all$result))
gostplot(genes.trans.m0.go.all, capped = TRUE, interactive = TRUE)


# M1
id.m1 <- which (all.trios$lond.model=="M1")
genes.cis.m1 <- unique (all.trios$cis.gene.name[id.m1])
genes.trans.m1 <- unique (all.trios$trans.gene.name[id.m1])

genes.cis.m1.go.all <- gost(query = genes.cis.m1, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m1.go.all$result)
# 178  15
genes.cis.m1.go.all$result$query <- rep ("m1.cis", nrow(genes.cis.m1.go.all$result))
gostplot(genes.cis.m1.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m1.go.all <- gost(query = genes.trans.m1, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.trans.m1.go.all$result)
# 3 15
genes.trans.m1.go.all$result$query <- rep ("m1.trans", nrow(genes.trans.m1.go.all$result))
gostplot(genes.trans.m1.go.all, capped = TRUE, interactive = TRUE)


# M2
id.m2 <- which (all.trios$lond.model=="M2")
genes.cis.m2 <- unique (all.trios$cis.gene.name[id.m2])
genes.trans.m2 <- unique (all.trios$trans.gene.name[id.m2])

genes.cis.m2.go.all <- gost(query = genes.cis.m2, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m2.go.all$result)
# 108  15
genes.cis.m2.go.all$result$query <- rep ("m2.cis", nrow(genes.cis.m2.go.all$result))
gostplot(genes.cis.m2.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m2.go.all <- gost(query = genes.trans.m2, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.trans.m2.go.all$result)
# 1 15
genes.trans.m2.go.all$result$query <- rep ("m2.trans", nrow(genes.trans.m2.go.all$result))
gostplot(genes.trans.m2.go.all, capped = TRUE, interactive = TRUE)


# M3
id.m3 <- which (all.trios$lond.model=="M3")
genes.cis.m3 <- unique (all.trios$cis.gene.name[id.m3])
genes.trans.m3 <- unique (all.trios$trans.gene.name[id.m3])

genes.cis.m3.go.all <- gost(query = genes.cis.m3, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m3.go.all$result)
# 2696   15
genes.cis.m3.go.all$result$query <- rep ("m3.cis", nrow(genes.cis.m3.go.all$result))
gostplot(genes.cis.m3.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m3.go.all <- gost(query = genes.trans.m3, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.trans.m3.go.all$result)
# 419  15
genes.trans.m3.go.all$result$query <- rep ("m3.trans", nrow(genes.trans.m3.go.all$result))
gostplot(genes.trans.m3.go.all, capped = TRUE, interactive = TRUE)


# M4
id.m4 <- which (all.trios$lond.model=="M4")
genes.cis.m4 <- unique (all.trios$cis.gene.name[id.m4])
genes.trans.m4 <- unique (all.trios$trans.gene.name[id.m4])

genes.cis.m4.go.all <- gost(query = genes.cis.m4, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
dim (genes.cis.m4.go.all$result)
# 11
genes.cis.m4.go.all$result$query <- rep ("m4.cis", nrow(genes.cis.m4.go.all$result))
gostplot(genes.cis.m4.go.all, capped = TRUE, interactive = TRUE)

genes.trans.m4.go.all <- gost(query = genes.trans.m4, organism = "hsapiens", exclude_iea = FALSE, correction_method = "bonferroni", highlight = TRUE)
# no results
# dim (genes.trans.m4.go.all$result)
# gostplot(genes.trans.m4.go.all, capped = TRUE, interactive = TRUE)


# Plot_genes_enrichment_all.pdf
gostplot(genes.cis.m0.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.trans.m0.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.cis.m1.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.trans.m1.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.cis.m2.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.trans.m2.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.cis.m3.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.trans.m3.go.all, capped = TRUE, interactive = FALSE)
gostplot(genes.cis.m4.go.all, capped = TRUE, interactive = FALSE)

# count the total number of genes for each enrichment analysis
length (genes.cis.m0)
# 21569
length (genes.trans.m0)
# 24744
length (genes.cis.m1)
# 2586
length (genes.trans.m1)
# 2211
length (genes.cis.m2)
# 1797
length (genes.trans.m2)
# 1210
length (genes.cis.m3)
# 22042
length (genes.trans.m3)
# 29733
length (genes.cis.m4)
# 2068
length (genes.trans.m4)
# 1687

# Merge all the enrichments 
genes.go.all <- rbind (genes.cis.m0.go.all$result, genes.trans.m0.go.all$result, genes.cis.m1.go.all$result, 
                       genes.trans.m1.go.all$result, genes.cis.m2.go.all$result, genes.trans.m2.go.all$result,
                       genes.cis.m3.go.all$result, genes.trans.m3.go.all$result, genes.cis.m4.go.all$result)
dim (genes.go.all)
# 8367   15

# The 'parents' column is a list.  Need to convert it to a vector
genes.go.all$parents <- sapply(genes.go.all$parents, function(x) paste(x, collapse = ","))

# Write to output
library(openxlsx)
write.xlsx(genes.go.all, "mediation_genes_GSEA_all.xlsx")
