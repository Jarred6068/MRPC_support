

#inferred mediation types across tissues (plots)
library(ggpubr)

#color palette choice 2:

#c("sienna1", "blue", "greenyellow", "grey52", "red3")

#read in tissue names and tables
tiss=read.csv("C:/Users/Bruin/Documents/GitHub/MRPC_support//ADDIS_ReRun/tissuenames.csv")

lond.table=read.csv(file="C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/TableS1_Trios_analysis_GTEx_v8_allPCs_V4_LOND.csv")
addis.table=read.csv(file="C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/TableS2_Trios_analysis_GTEx_v8_allPCs_V3_ADDIS.csv")

#adjust data format lond
longform=expand.grid(tiss[,2], c("M0","M1","M2","M3","M4","Others"))
longdata=as.vector(cbind(lond.table$Num.M0, lond.table$Num.M1, lond.table$Num.M2, lond.table$Num.M3, lond.table$Num.M4, lond.table$Num.Others))
data.new.lond=cbind.data.frame(longform, longdata)
colnames(data.new.lond)=c("Tissue_Name", "Model","Number of Trios")

data.new.lond$Model=as.factor(data.new.lond$Model)
data.new.lond$logtrios=log(data.new.lond$`Number of Trios`, base = 10)

#plot for lond
png("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F3_Breakdown of Inferred Trio Types (LOND).png", units = "in",
    width=12, height=8, res=1000)
#tiff("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F3_Breakdown of Inferred Trio Types (LOND).tiff", units = "in",
#    width=12, height=8, res=800)

ggbarplot(data.new.lond, x = "Tissue_Name", y = "Number of Trios",
          fill = "Model", color = "Model", 
          palette = c("cadetblue1","red4", "yellow","navyblue","salmon1", "chartreuse3"),
          label = FALSE, ylab = FALSE, order = sort(tiss[,2], decreasing = T))+
  theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7), 
        axis.text.y = element_text(size=9))+
  scale_y_continuous(breaks = seq(0, 14000, 2000), lim = c(0, 14000))+
  coord_flip()+
  ggtitle("Breakdown of Inferred Trio Types (LOND)")

dev.off()

ggsave("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F3_Breakdown of Inferred Trio Types (LOND).pdf", units = "in",
    width=12, height=8)

#adjust data format for addis:
longform=expand.grid(tiss[,2], c("M0","M1","M2","M3","M4","Others"))
longdata=as.vector(cbind(addis.table$Num.M0, addis.table$Num.M1, addis.table$Num.M2, addis.table$Num.M3, addis.table$Num.M4, addis.table$Num.Others))
data.new.addis=cbind.data.frame(longform, longdata)
colnames(data.new.addis)=c("Tissue_Name", "Model","Number of Trios")

data.new.addis$Model=as.factor(data.new.addis$Model)

#plot for addis

png("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F3_Breakdown of Inferred Trio Types (ADDIS).png", units = "in",
    width=12, height=8, res=1000)
#tiff("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F3_Breakdown of Inferred Trio Types (ADDIS).tiff", units = "in",
#    width=12, height=8, res=800)
ggbarplot(data.new.addis, x = "Tissue_Name", y = "Number of Trios",
          fill = "Model", color = "Model", 
          palette = c("cadetblue1", "red4", "yellow","navyblue","salmon1", "chartreuse3"),
          label = FALSE,ylab = FALSE,order = sort(tiss[,2], decreasing = T))+
  theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7), 
        axis.text.y = element_text(size=9))+
  scale_y_continuous(breaks = seq(0, 14000, 2000), lim = c(0, 14000))+
  coord_flip()+
  ggtitle("Breakdown of Inferred Trio Types (ADDIS)")

dev.off()

ggsave("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F3_Breakdown of Inferred Trio Types (ADDIS).pdf", units = "in",
       width=12, height=8)
#===========================================================================================
#--------------------Cis--Trans--inferred--M1--Mediator--breakdown--------------------------
#===========================================================================================
#source("C:/Users/Bruin/Documents/GitHub/MRPC_support/Reg_Net/AL_genetabV2.R")
#read in ADDIS_Post_Analysis_processing files
tiss=read.csv("C:/Users/Bruin/Documents/GitHub/MRPC_support//ADDIS_ReRun/tissuenames.csv")

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#load the master summary table for M1 trios (found in ./Reg_Net)
df=loadRData(fileName="C:/Users/Bruin/Documents/GitHub/MRPC_support/Reg_Net/loaded_AL_datatables.Rdata")

#partition the lond and addis data
addis.table=subset(df$mdf, df$mdf$Mediator=="AM1T1" | df$mdf$Mediator=="AM1T2")
lond.table=subset(df$mdf, df$mdf$Mediator=="LM1T1" | df$mdf$Mediator=="LM1T2")

#allocate space
addis.cis.sumtable=as.data.frame(matrix(0, nrow=dim(tiss)[1], ncol = 3))
addis.trans.sumtable=as.data.frame(matrix(0, nrow=dim(tiss)[1], ncol = 3))
lond.cis.sumtable=as.data.frame(matrix(0, nrow=dim(tiss)[1], ncol = 3))
lond.trans.sumtable=as.data.frame(matrix(0, nrow=dim(tiss)[1], ncol = 3))


#reformat data table so it matches the "longform" from above
for(i in 1:dim(tiss)[1]){
  
  addis.cis.sumtable[i,]=c(tiss[i,2], "via cis gene", dim(subset(addis.table, addis.table$Mediator=="AM1T1" & addis.table$tissue==tiss[i,2]))[1])
  addis.trans.sumtable[i,]=c(tiss[i,2], "via trans gene", dim(subset(addis.table, addis.table$Mediator=="AM1T2" & addis.table$tissue==tiss[i,2]))[1])
  lond.cis.sumtable[i,]=c(tiss[i,2], "via cis gene", dim(subset(lond.table, lond.table$Mediator=="LM1T1" & lond.table$tissue==tiss[i,2]))[1])
  lond.trans.sumtable[i,]=c(tiss[i,2], "via trans gene", dim(subset(lond.table, lond.table$Mediator=="LM1T2" & lond.table$tissue==tiss[i,2]))[1])
  
}

#combine, give colnames, and make type a factor:
addis.table.final=rbind.data.frame(addis.cis.sumtable, addis.trans.sumtable)
colnames(addis.table.final)=c("Tissue_Name", "Type of Mediation", "Number of Trios")
addis.table.final$`Type of Mediation`=as.factor(addis.table.final$`Type of Mediation`)
addis.table.final$`Number of Trios`=as.numeric(addis.table.final$`Number of Trios`)

lond.table.final=rbind.data.frame(lond.cis.sumtable, lond.trans.sumtable)
colnames(lond.table.final)=c("Tissue_Name", "Type of Mediation", "Number of Trios")
lond.table.final$`Type of Mediation`=as.factor(lond.table.final$`Type of Mediation`)
lond.table.final$`Number of Trios`=as.numeric(lond.table.final$`Number of Trios`)
#plot addis M1 inferred

#plot for addis M1 inferred

png("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F4_Breakdown of Inferred M1 Mediation Types (ADDIS).png", units = "in",
    width=12, height=8, res=1000)
#tiff("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F4_Breakdown of Inferred M1 Mediation Types (ADDIS).tiff", units = "in",
#    width=12, height=8, res=800)
ggbarplot(addis.table.final, x = "Tissue_Name", y = "Number of Trios",
          fill = "Type of Mediation", color = "Type of Mediation", 
          palette = c("navyblue","salmon1"),
          label = FALSE, ylab = FALSE, order = sort(tiss[,2], decreasing = T))+
  theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7), 
        axis.text.y = element_text(size=8))+
  scale_y_continuous(breaks = seq(0, 150, 10), lim = c(0, 150))+
  coord_flip()+
  ggtitle("Breakdown of Inferred M1 Mediation Types (ADDIS)")


dev.off()
ggsave("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F4_Breakdown of Inferred M1 Mediation Types (ADDIS).pdf", units = "in",
    width=12, height=8)



png("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F4_Breakdown of Inferred M1 Mediation Types (LOND).png", units = "in",
    width=12, height=8, res=1000)
#tiff("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F4_Breakdown of Inferred M1 Mediation Types (LOND).tiff", units = "in",
#    width=12, height=8, res=800)

#plot lond M1 inferred
ggbarplot(lond.table.final, x = "Tissue_Name", y = "Number of Trios",
          fill = "Type of Mediation", color = "Type of Mediation", 
          palette = c("navyblue","salmon1"),
          label = FALSE, ylab = FALSE, order = sort(tiss[,2], decreasing = T))+
  theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7), 
        axis.text.y = element_text(size=8))+
  scale_y_continuous(breaks = seq(0, 150, 10), lim = c(0, 150))+
  coord_flip()+
  ggtitle("Breakdown of Inferred M1 Mediation Types (LOND)")


dev.off()
ggsave("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F4_Breakdown of Inferred M1 Mediation Types (LOND).pdf", units = "in",
    width=12, height=8)



#===========================================================================================
#--------------------Cis--Trans--hists--of--shared--tissues---------------------------------
#===========================================================================================

source("C:/Users/Bruin/Documents/GitHub/MRPC_support/Reg_Net/AL_genetabV2.R")
tiss=read.csv("C:/Users/Bruin/Documents/GitHub/MRPC_support//ADDIS_ReRun/tissuenames.csv")

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#load the master summary table for M1 trios (found in ./Reg_Net)
df=loadRData(fileName="C:/Users/Bruin/Documents/GitHub/MRPC_support/Reg_Net/loaded_AL_datatables.Rdata")


#counts and bins for 4 subsets (Addis T1, Addis T2, Lond T1, Lond T2)
#uses functions in AL_genetabV2.R: binit() and count()
LTasM=binit(df$lm1t2, target=4)$bins
LCasM=binit(df$lm1t1, target=3)$bins

ATasM=binit(df$am1t2, target=4)$bins
ACasM=binit(df$am1t1, target=3)$bins


#reformat for addis:
ATasM$`Type of Mediation`=rep("via trans gene",48)
ACasM$`Type of Mediation`=rep("via cis gene",48)

ATCM=rbind.data.frame(ATasM, ACasM)
ATCM.final=subset(ATCM, ATCM$`Number of Genes`>0)
ATCM.final$`Number of Genes`=log(ATCM.final$`Number of Genes`+1, base = 10)


#reformat for long
LTasM$`Type of Mediation`=rep("via trans gene",48)
LCasM$`Type of Mediation`=rep("via cis gene",48)

LTCM=rbind.data.frame(LTasM, LCasM)
LTCM.final=subset(LTCM, LTCM$`Number of Genes`>0)
LTCM.final$`Number of Genes`=log(LTCM.final$`Number of Genes`+1, base = 10)

library(latex2exp)

png("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F5_Histograms of Tissue Sharing (ADDIS).png", units = "in",
    width=12, height=8, res=1000)
#tiff("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F5_Histograms of Tissue Sharing (ADDIS).tiff", units = "in",
#    width=12, height=8, res=800)

# plot for addis
ggbarplot(ATCM.final, x = "Number of Tissues Shared", y = "Number of Genes",
          fill = "Type of Mediation", color = "Type of Mediation", 
          palette = c("navyblue","salmon1"),
          label = FALSE, ylab = TeX("$\\log_{10}($Number of Genes $ + 1)$"), size = 0.2,
          position = position_dodge2(0.9, preserve = "single"))+
  theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7), 
        axis.text.y = element_text(size=10))+
  scale_x_continuous(breaks = seq(0, 48, 4), lim = c(0, 48))+
  scale_y_continuous(breaks = seq(0, round(max(ATCM.final$`Number of Genes`))+0.5, 0.5), 
                     lim = c(0, round(max(ATCM.final$`Number of Genes`))+0.5))+
  #coord_flip()+
  ggtitle("Number of Shared Tissues Among Trans and Cis Mediation Genes (ADDIS)")

dev.off()

ggsave("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F5_Histograms of Tissue Sharing (ADDIS).pdf", units = "in",
    width=12, height=8)





png("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F5_Histograms of Tissue Sharing (LOND).png", units = "in",
    width=12, height=8, res=1000)
#tiff("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F5_Histograms of Tissue Sharing (LOND).tiff", units = "in",
#    width=12, height=8, res=800)
#plot for Lond
ggbarplot(LTCM.final, x = "Number of Tissues Shared", y = "Number of Genes",
          fill = "Type of Mediation", color = "Type of Mediation", 
          palette = c("navyblue","salmon1"),
          label = FALSE, ylab = TeX("$\\log_{10}($Number of Genes $ + 1)$"), size = 0.2,
          position = position_dodge2(0.9, preserve = "single"))+
  theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7), 
        axis.text.y = element_text(size=10))+
  scale_x_continuous(breaks = seq(0, 48, 4), lim = c(0, 48))+
  scale_y_continuous(breaks = seq(0, round(max(ATCM.final$`Number of Genes`))+0.5, 0.5), 
                     lim = c(0, round(max(ATCM.final$`Number of Genes`))+0.5))+
  #coord_flip()+
  ggtitle("Number of Shared Tissues Among Trans and Cis Mediation Genes (LOND)")

dev.off()

ggsave("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F5_Histograms of Tissue Sharing (LOND).pdf", units = "in",
    width=12, height=8)



#===========================================================================================
#-----------------------Cis--Trans--Gene--Types--Breakdown----------------------------------
#===========================================================================================



source("C:/Users/Bruin/Documents/GitHub/MRPC_support/Reg_Net/AL_genetabV2.R")
tiss=read.csv("C:/Users/Bruin/Documents/GitHub/MRPC_support//ADDIS_ReRun/tissuenames.csv")

loadRData <- function(fileName=NULL){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#load the master summary table for M1 trios (found in ./Reg_Net)
df=loadRData(fileName="C:/Users/Bruin/Documents/GitHub/MRPC_support/Reg_Net/loaded_AL_datatables.Rdata")
BM=read.table(file="C:/Users/Bruin/Documents/GitHub/MRPC_support/Reg_Net/mart_export_merged_lncRNA_fixed.txt", sep="\t", header=T)
test.strs=c("pseudogene", "protein_coding", "lncRNA")


convert_cats=function(gene_types){
  test.strs=c("pseudogene", "protein_coding", "lncRNA")
  types.fixed=matrix(0, nrow = length(gene_types), ncol = 3)

  for(i in 1:length(test.strs)){

    logi=grepl(test.strs[i], gene_types, fixed = TRUE)
    gene_types[logi]=test.strs[i]
  }
  
  others=setdiff(c(1:length(gene_types)),
                 c(which(gene_types==test.strs[1]), 
                   which(gene_types==test.strs[2]), 
                   which(gene_types==test.strs[3])))
  
  gene_types[others]="other"
  
  
  return(gene_types)

}

#Addis convert gene type categories
df$am1t1$cis.Gene.type=convert_cats(df$am1t1$cis.Gene.type)
df$am1t2$trans.Gene.type=convert_cats(df$am1t2$trans.Gene.type)
#Lond convert gene type categories
df$lm1t1$cis.Gene.type=convert_cats(df$lm1t1$cis.Gene.type)
df$lm1t2$trans.Gene.type=convert_cats(df$lm1t2$trans.Gene.type)



#Addis

addis.data.t1=cbind.data.frame(df$am1t1$cis.Gene.type, df$am1t1$tissue, rep("via cis gene", dim(df$am1t1)[1]))
addis.data.t2=cbind.data.frame(df$am1t2$trans.Gene.type, df$am1t2$tissue, rep("via trans gene", dim(df$am1t2)[1]))

#Lond

lond.data.t1=cbind.data.frame(df$lm1t1$cis.Gene.type, df$lm1t1$tissue, rep("via cis gene", dim(df$lm1t1)[1]))
lond.data.t2=cbind.data.frame(df$lm1t2$trans.Gene.type, df$lm1t2$tissue, rep("via trans gene", dim(df$lm1t2)[1]))

colnames(lond.data.t1)=c("Gene Type", "Tissue", "Type of Mediation")
colnames(lond.data.t2)=c("Gene Type", "Tissue", "Type of Mediation")
colnames(addis.data.t1)=c("Gene Type", "Tissue", "Type of Mediation")
colnames(addis.data.t2)=c("Gene Type", "Tissue", "Type of Mediation")

#convert to a counts table:

#preallocate new tables
gene.types.list=c("pseudogene", "protein_coding", "lncRNA", "other")

count.at1=c()
count.at2=c()
count.lt1=c()
count.lt2=c()

#get the counts for type1 and type2 mediator gene types
for(i in 1:4){
  
  minict.at1=rep(0, 48)
  minict.at2=rep(0, 48)
  minict.lt1=rep(0, 48)
  minict.lt2=rep(0, 48)
  
  for(j in 1:48){
    
    #addis t1 and t2
    minict.at1[j]=dim(subset(addis.data.t1, addis.data.t1$`Gene Type`==gene.types.list[i] & addis.data.t1$Tissue==tiss[j,2]))[1]
    minict.at2[j]=dim(subset(addis.data.t2, addis.data.t2$`Gene Type`==gene.types.list[i] & addis.data.t2$Tissue==tiss[j,2]))[1]
    #lond t1 and t2
    minict.lt1[j]=dim(subset(lond.data.t1, lond.data.t1$`Gene Type`==gene.types.list[i] & lond.data.t1$Tissue==tiss[j,2]))[1]
    minict.lt2[j]=dim(subset(lond.data.t2, lond.data.t2$`Gene Type`==gene.types.list[i] & lond.data.t2$Tissue==tiss[j,2]))[1]
    
  }
  
  #append addis
  count.at1=c(count.at1, minict.at1)
  count.at2=c(count.at2, minict.at2)
  #append lond
  count.lt1=c(count.lt1, minict.lt1)
  count.lt2=c(count.lt2, minict.lt2)
  
}

#complete table Addis
addis.data.t1.new=cbind.data.frame(count.at1, expand.grid(tiss[,2], gene.types.list), rep("Mediation via cis gene", length(count.at1)))
addis.data.t2.new=cbind.data.frame(count.at2, expand.grid(tiss[,2], gene.types.list), rep("Mediation via trans gene", length(count.at2)))
colnames(addis.data.t1.new)=colnames(addis.data.t2.new)
addis.final=rbind(addis.data.t1.new, addis.data.t2.new)
#colnames
colnames(addis.final)=c("Number of Genes", "Tissue", "Gene Type", "Type_Mediation")
#complete table Lond
lond.data.t1.new=cbind.data.frame(count.lt1, expand.grid(tiss[,2], gene.types.list), rep("Mediation via cis gene", length(count.lt1)))
lond.data.t2.new=cbind.data.frame(count.lt2, expand.grid(tiss[,2], gene.types.list), rep("Mediation via trans gene", length(count.lt2)))
colnames(lond.data.t1.new)=colnames(lond.data.t2.new)
lond.final=rbind(lond.data.t1.new, lond.data.t2.new)
#colnames
colnames(lond.final)=c("Number of Genes", "Tissue", "Gene Type", "Type_Mediation")


png("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F6_Distribution of Mediator Gene Types For each Tissue (LOND).png", units = "in",
    width=12, height=8, res=1000)
#tiff("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F6_Distribution of Mediator Gene Types For each Tissue (LOND).tiff", units = "in",
#    width=12, height=8, res=800)
#plot for lond
ggbarplot(lond.final, x = "Tissue", y = "Number of Genes",
          fill = "Gene Type", color = "Gene Type", 
          palette = c("red4","salmon1","cadetblue1", "navyblue"),
          label = FALSE, ylab = FALSE, order = sort(tiss[,2], decreasing = T),
          facet.by = "Type_Mediation")+
  
  theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7), 
        axis.text.y = element_text(size=9))+
  
  scale_y_continuous(breaks = seq(0, 120, 10), lim = c(0, 120))+
  
  coord_flip()+
  ggtitle("Distribution of Mediator Gene Types For each Tissue (LOND)")
dev.off()

ggsave("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F6_Distribution of Mediator Gene Types For each Tissue (LOND).pdf", units = "in",
    width=12, height=8)



png("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F6_Distribution of Mediator Gene Types For each Tissue (ADDIS).png", units = "in",
    width=12, height=8, res=1000)
#tiff("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F6_Distribution of Mediator Gene Types For each Tissue (ADDIS).tiff", units = "in",
#    width=12, height=8, res=800)

#plot for addis
ggbarplot(addis.final, x = "Tissue", y = "Number of Genes",
          fill = "Gene Type", color = "Gene Type", 
          palette = c("red4","salmon1","cadetblue1", "navyblue"),
          label = FALSE, ylab = FALSE, order = sort(tiss[,2], decreasing = T),
          facet.by = "Type_Mediation")+
  
  theme(axis.text.x = element_text(size=10, hjust=0.5,vjust=0.7), 
        axis.text.y = element_text(size=9))+
  
  scale_y_continuous(breaks = seq(0, 120, 10), lim = c(0, 120))+
  
  coord_flip()+
  ggtitle("Distribution of Mediator Gene Types For each Tissue (ADDIS)")

dev.off()

ggsave("C:/Users/Bruin/Documents/GitHub/MRPC_support/Manuscript/F6_Distribution of Mediator Gene Types For each Tissue (ADDIS).pdf", units = "in",
    width=12, height=8)
