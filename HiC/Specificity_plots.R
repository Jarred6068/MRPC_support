source("/mnt/ceph/jarredk/HiC_Analyses/HiC_search.R")
library(gtools, lib = '/mnt/ceph/jarredk/Rpackages')
library(ggpubr, lib = '/mnt/ceph/jarredk/Rpackages')

fibro.sig = read.csv('/mnt/ceph/jarredk/HiC_Analyses/fibroblast_sig_trios.csv')


chrm.name = paste(rep('chr', 24), as.character(c(1:22, 'X', "Y")), sep = '')
chrm = as.character(c(1:22, 'X', "Y"))
coords = c(1:24)

combos = gtools::combinations(length(chrm), 2, coords, repeats.allowed = T)
combos2 = gtools::combinations(length(chrm), 2, chrm, repeats.allowed = T)
combos3 = gtools::combinations(length(chrm), 2, chrm.name, repeats.allowed = T)

get.count = function(chr1, chr2, fname="/mnt/ceph/jarredk/HiC_Analyses/Data/lymphoblastoid_cells/ENCFF355OWW.hic", 
                     resolution = 5000){
  
  hic1=extract_hic(fileName=fname, 
                   chrs=c(chr1,chr2), 
                   resol=resolution)
  
  ct = log10(sum(hic1$counts))
  
  return(ct)
  
}




int.check=function(trio.data, hic.filename=NULL, resolution=10000, search.size=100000, 
                   tiss="CellsEBVtransformedlymphocytes", verbose=TRUE, pack.path="/mnt/ceph/jarredk/Rpackages", 
                   plot.h=FALSE, FDR=NULL){
  #SYNTAX:
  #hic.filename -- the path to the desired .hic file
  #trio -- the trios index (or vector of indexes) typically obtained from ADDIS.M1.check() or ADDIS.PostProc() 
  #resolution -- the BP resolution desired and to be passed to "extract_hic() and Resample_interactions()"
  #search.size -- the size of the bin around gene positions: also passed to Resample_interactions()
  #tiss -- a string specifying the name of the tissue to which the trio belongs 
  #Verbose -- logical which if true prints checkpoints as trios are checked
  #pack.path -- the path to you Rpackages library for dependencies such as StrawR
  #plot.h -- logical passed to Resample_interactions()
  #FDR -- string specifying "LOND" or "ADDIS" passed to Resample_interactions()
  
  #pre-allocate
  num.trios = dim(trio.data)[1]
  reads = NULL
  averages = NULL
  p.values = NULL
  p.values2 = NULL
  total.nas = NULL
  resampled_dataset=list()
  
  hic.extent=as.data.frame(matrix(0, nrow = num.trios, ncol = 4))
  colnames(hic.extent)=c("lower_x", "upper_x", "lower_y", "upper_y")
  
  for(i in 1:num.trios){
    
    cis.chr = trio.data$Cis.Gene.Chr[i]
    var.pos = trio.data$Variant.Position[i]
    
    trans.chr = trio.data$Trans.Gene.Chr[i]
    trans.left = trio.data$Trans.Gene.Left.End[i]
    trans.right = trio.data$Trans.Gene.Right.End[i]
    
    
    trio.name = paste(trio.data$Cis.Gene.ID[i], trio.data$Trans.Gene.ID[i], sep=':')
    
    print(paste0('Preparing trio :', trio.name))
    
    #establish snp - gene search "box"
    left.pos=paste(cis.chr,":",
                   var.pos-search.size,":",
                   var.pos+search.size, sep = "")
    
    #establish trans gene search "box"
    right.pos=paste(trans.chr,":",
                    trans.left-search.size,":",
                    trans.right+search.size, sep = "")
    
    
    if(verbose==TRUE){
      print(paste('Start: ', trans.left, ' End: ', trans.right))
    }
    
    
    #save search upper and lower bounds
    hic.extent[i,]=c(var.pos-search.size,
                     var.pos+search.size,
                     trans.left-search.size,
                     trans.right+search.size)
    
    #extract data from .hic based on search parameters
    data.hic=extract_hic(fileName = hic.filename, 
                         chrs = c(left.pos, right.pos),
                         resol=resolution)
    
    #count interactions in observed
    reads[i]=ifelse(plyr::empty(as.data.frame(data.hic))==TRUE, NA, sum(data.hic$counts))
    
    
    #resample the data to obtain a probability distribution by MC integration
    if(is.na(reads[i])==TRUE){
      #if no reads for variant<-->gene return NA's
      averages[i]=NA
      p.values[i]=NA
      p.values2[i]=NA
      total.nas[i] = NA
    }else{
      #else return the resampled data and determine P(>=obs)
      RS=Resample_interactions(filePath = hic.filename,
                               chrs=c(paste(cis.chr),paste(trans.chr)),
                               res=10000,
                               search.size=search.size,
                               tiss = tiss,
                               trio = trio.name,
                               plot.hist = plot.h,
                               FDR = FDR)
      
      # obtain resampled data, extract NA's and non-NA data points and calc pvalue
      averages[i]=RS$confInterval[2]
      totals=as.vector(na.omit(RS$resampled.totals))
      num_nas=as.vector(attr(RS$resampled.totals,"na.action"))
      
      #MC-integration
      vec=ifelse(totals>=reads[i], 1, 0)
      #P(>obs) = P(>obs and NA) + P(>obs and !NA) = 0 + P(>obs)*P(!NA)
      p.values[i]=( sum(na.omit(vec))/length(totals) )
      p.values2[i]=( sum(na.omit(vec))/(length(totals)+length(num_nas)) )
      
      #return the resampled data set
      resampled_dataset[[i]]=list(sampled=totals, nas=num_nas)
      total.nas[i]=length(num_nas)
      
    }
    
    
    
    
    
  }
  
  
  
  #preallocate space for FDR/FWER adjustments
  bh.thresh=rep(0, length(p.values))
  qvals=rep(0, length(p.values))
  BY=rep(0, length(p.values))
  
  bh.thresh2=rep(0, length(p.values2))
  qvals2=rep(0, length(p.values2))
  BY2=rep(0, length(p.values2))
  
  #print(dim(hic.extent))
  print("made it this far1")
  
  #print(colnames(trio.attr$Attributes$trans))
  #print(head(trio.attr$Attributes$trans))
  #print(colnames(trio.attr$Attributes$trans))
  #print(head(trio.attr$Attributes$trans))
  #summarize
  info.list = cbind.data.frame(trio.data$SNP, 
                               trio.data$Trans.Gene.Name,
                               trio.data$Trans.Gene.ID,
                               trio.data$Trans.Gene.Chr,
                               reads, averages, total.nas, p.values, BY,
                               bh.thresh, qvals, 
                               trio.data$Variant.Position,
                               trio.data$Trans.Gene.Left.End,
                               trio.data$Trans.Gene.Right.End,
                               hic.extent, 
                               trio.data$Cis.Gene.Name,
                               trio.data$Cis.Gene.ID,
                               trio.data$Cis.Gene.Chr)
  
  #print(info.list)
  
  info.list2 = cbind.data.frame(trio.data$SNP, 
                               trio.data$Trans.Gene.Name,
                               trio.data$Trans.Gene.ID,
                               trio.data$Trans.Gene.Chr,
                               reads, averages, total.nas, p.values2, BY2,
                               bh.thresh2, qvals2, 
                               trio.data$Variant.Position,
                               trio.data$Trans.Gene.Left.End,
                               trio.data$Trans.Gene.Right.End,
                               hic.extent, 
                               trio.data$Cis.Gene.Name,
                               trio.data$Cis.Gene.ID,
                               trio.data$Cis.Gene.Chr)
  #name cols
  cnames=c("SNP", "trans.gene.name", "trans.gene.ID", "trans.chr", "obs.reads", "expected.reads", "total_NA's", "P(>obs)", 
           "BY","HB.Adjusted", "qvals","variant.pos", "trans.left","trans.right", "variant.lower.bound", 
           "variant.upper.bound", "trans.lower.bound", "trans.upper.bound","cis.gene.name", "cis.gene.ID", "cis.chr")
  
  colnames(info.list)=cnames
  colnames(info.list2)=cnames
  
  
  #Holm-Bonferroni correction at FWER alpha = 0.05
  #initialize
  alpha=0.05
  m=length(na.omit(p.values))
  sorted.p=sort(p.values, index.return=TRUE, na.last = TRUE)
  sorted.p2=sort(p.values2, index.return=TRUE, na.last = TRUE)
  print(sorted.p$x)
  HB.adjust=rep(0, length(p.values))
  
  #calculate the rejections using step-down procedure
  for(k in 1:m){
    HB.adjust[k]=sorted.p$x[k]*(m+1-k)
  }
  
  #store in output table
  info.list=info.list[sorted.p$ix,]
  info.list2=info.list2[sorted.p2$ix,]
  HB.adjust=ifelse(HB.adjust>1, 1, HB.adjust)
  HB.adjust=ifelse(HB.adjust==0, NA, HB.adjust)
  info.list$HB.Adjusted=HB.adjust
  
  info.list2$HB.Adjusted=p.adjust(sorted.p2$x, method = "holm")
  
  #get qvalues using BH step-up procedure
  library(qvalue, lib=pack.path)
  
  Q=qvalue(as.vector(na.omit(info.list$`P(>obs)`)), fdr.level = alpha, pi0 = 1)
  Q2=qvalue(as.vector(na.omit(info.list2$`P(>obs)`)), fdr.level = alpha, pi0 = 1)
  
  Qq=as.vector(Q$qvalues)
  Qq2=as.vector(Q2$qvalues)
  
  Qq=ifelse(is.na(info.list$`P(>obs)`)==TRUE, NA, Qq)
  Qq2=ifelse(is.na(info.list2$`P(>obs)`)==TRUE, NA, Qq2)
  info.list$qvals=Qq
  info.list2$qvals=Qq2
  
  #include the Hochberg step up correction
  
  info.list$BY=p.adjust(sorted.p$x, method="BY")
  info.list2$BY=p.adjust(sorted.p2$x, method = "BY")
  
  #return as list
  return(list(summary.table1=info.list, summary.table2=info.list2, data=resampled_dataset))
  
  
  
}





tisspath = c("/mnt/ceph/jarredk/HiC_Analyses/Data/lymphoblastoid_cells/ENCFF355OWW.hic",
             "/mnt/ceph/jarredk/HiC_Analyses/Data/Skin/ENCFF569RJM.hic",
             "/mnt/ceph/jarredk/HiC_Analyses/Data/Lung/ENCFF366ERB.hic",
             "/mnt/ceph/jarredk/HiC_Analyses/Data/fibroblast_cells/ENCFF768UBD.hic")


tissues = c("CellsEBVtransformedlymphocytes",
            "SkinNotSunExposed",
            "Lung",
            "CellsCulturedfibroblasts")

num.trios = dim(fibro.sig)[1]                         


ftab = data.frame(Tissue = character(),
                  mean = numeric(),
                  sd = numeric(),
                  stringsAsFactors = FALSE)
colnames(ftab) = c('Tissue', 'mean', 'sd')

final.list = vector('list', length = 4)
names(final.list) = tissues

for(p in 1:length(tisspath)){
  
  path = tisspath[p]
  tissue = tissues[p]
  
  
  # print(paste0('extracting interactions from: ', path))
  # output = int.check(trio.data = fibro.sig, 
  #                             hic.filename=path, 
  #                             resolution=10000, 
  #                             search.size=100000, 
  #                             tiss=tisue, 
  #                             verbose=TRUE, 
  #                             pack.path="/mnt/ceph/jarredk/Rpackages", 
  #                             plot.h=FALSE, 
  #                             FDR='LOND')
  # 
  #   final.list[[p]] = output    
  #   
  # 
  # 
  # for(tt in 1:num.trios){
  #   
  #   trio.ch1 = as.character(fibro.sig[tt, ]$Trans.Gene.Chr)
  #   trio.ch2 = as.character(fibro.sig[tt, ]$Cis.Gene.Chr)
  #   
  #   trio.name = paste(fibro.sig[tt, ]$Trans.Gene.ID, fibro.sig[tt, ]$Cis.Gene.ID, sep = ':')
  #   
  #   print(trio.name)
  #   print(paste0('Cis chr: ', trio.ch2, ', Trans chr: ', trio.ch1))
  #   
  #   output = Resample_interactions(filePath=path, 
  #                                  chrs=c(trio.ch1,trio.ch2), 
  #                                  res=10000, 
  #                                  search.size=100000, 
  #                                  resamples=10000, 
  #                                  verbose=TRUE, 
  #                                  plot.hist=TRUE, 
  #                                  trio=trio.name, 
  #                                  tiss=tissue, 
  #                                  FDR='LOND', 
  #                                  plot.title=NULL)
  #   
  #   
  #   #print(str(output))
  #   print(output$confInterval)
  #   
  #   ftab[nrow(ftab) + 1, ] = c(tissue, mean(output$resampled.totals), sd(output$resampled.totals))
  #   
  #   print(ftab)
  #   
  # print(paste0('finished checking interactions of trio '))
  # save(final.list, file = paste0('/mnt/ceph/jarredk/HiC_Analyses/sigfibro_tissue_interactions.RData'))
  #   
  #   
  #   #preallocate
  mat1 = as.data.frame(matrix(0, nrow = 24, ncol = 24))
  long.mat = as.data.frame(matrix(rep(0, dim(combos)[1]*3), nrow =dim(combos)[1], ncol = 3))
  colnames(long.mat) = c('Chromosome i', 'Chromosome j', 'Log10(# of interactions)')
    
  for(i in 1:dim(combos)[1]){
    
    r = combos[i,1]
    c = combos[i,2]
    ch1 = combos2[i,1]
    ch2 = combos2[i,2]
    
    
    xcoord = ifelse(ch1 %in% c('X', 'Y'), ifelse(ch1 == 'X', 23, 24), as.numeric(ch1))
    ycoord = ifelse(ch2 %in% c('X', 'Y'), ifelse(ch2 == 'X', 23, 24), as.numeric(ch2))
   
    val = get.count(chr1 = ch1, chr2 = ch2, fname = path)
    
    ct = ifelse(length(val) == 0, 0, val)
    
    print(paste0('found ', ct, ' interactions between ', ch1, ' and ', ch2))
    
    long.mat[i,] = c(combos3[i,], ct) 
    
    mat1[xcoord,ycoord] = ct
    
  }
  
  
  long.mat$`Log10(# of interactions)` = as.numeric(long.mat$`Log10(# of interactions)`)
  
  print('plotting interaction matrix')
  A = ggplot(long.mat, aes(`Chromosome i`, `Chromosome j`)) +
    geom_tile(aes(fill = `Log10(# of interactions)`)) +
    geom_text(aes(label = round(`Log10(# of interactions)`, 2)), size = 2) +
    scale_fill_gradient(low = "white", high = "red")+
    theme(legend.position = 'top',
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 90))
  
  
  
  pdf(paste0('/mnt/ceph/jarredk/HiC_Analyses/',tissue,'_histogram.pdf'))
  plot(A)
  dev.off()
  
}
    


  

  
  
  
  

