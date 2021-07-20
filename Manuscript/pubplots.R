

library(ggpubr)


lond.table=read.csv(file="C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/TableS1_Trios_analysis_GTEx_v8_allPCs_V4_LOND.csv")
addis.table=read.csv(file="C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/TableS2_Trios_analysis_GTEx_v8_allPCs_V3_ADDIS.csv")


longform=expand.grid(lond.table$Tissue.name, c("M0","M1","M2","M3","M4"))
longdata=as.vector(cbind(lond.table$Num.M0, lond.table$Num.M1, lond.table$Num.M2, lond.table$Num.M3, lond.table$Num.M4))
data.new.lond=cbind.data.frame(longform, longdata)
colnames(data.new.lond)=c("Tissue_Name", "Model_type","Inferred_Count")

data.new.lond$Model_type=as.factor(data.new.lond$Model_type)

ggbarplot(data.new.lond, x = "Tissue_Name", y = "Inferred_Count",
          fill = "Model_type", color = "Model_type", 
          palette = "uchicago",
          label = FALSE)+theme(axis.text.x = element_text(angle = 90))

























