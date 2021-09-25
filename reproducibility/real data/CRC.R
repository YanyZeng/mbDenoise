setwd("~/mbDenoise/reproducibility")

##load packages/methods source
library(curatedMetagenomicData)
library("phyloseq")
library(ade4)
library(edgeR)
library(DESeq2)
library(readxl)
library(openxlsx)
library(ape)
library(metagenomeSeq)
library(coin)
library(reticulate)
library(mbImpute)
library(glmnet)
library(SAVER)
library(DirichletMultinomial)
library('Matrix')
library('parallel')
library('doParallel')
library('foreach')
library('optimx')
library('trustOptim')
library('ZIPFA')
library(mbDenoise)

source("other_software/proximal_gradient.R")
source("other_software/basic_functions.R")
source("other_software/tune_proximal_gradient.R")
source("other_software/naive_methods.R")
source_python('other_software/ZIFA_copy.py')

##Figures 10 and 11 data processing
{
##ZellerG_2014 dataset
{
  ZellerG_2014 <- curatedMetagenomicData("ZellerG_2014.metaphlan_bugs_list.stool",
                                         counts=TRUE, dryrun=FALSE, bugs.as.phyloseq=TRUE)
  ZellerG_2014_genus <- phyloseq::tax_glom(ZellerG_2014$ZellerG_2014.metaphlan_bugs_list.stool, taxrank="Species")
  ZellerG_2014_tax <- tax_table(ZellerG_2014_genus )
  ZellerG_2014_otu <- as.data.frame(t(otu_table(ZellerG_2014_genus)))
  ZellerG_2014_sample <- as.data.frame( sample_data(ZellerG_2014_genus))

  ZellerG_2014_sample <- ZellerG_2014_sample[ZellerG_2014_sample[,4]!="adenoma",]
  ZellerG_2014_sample <- ZellerG_2014_sample[complete.cases(ZellerG_2014_sample[,10]),]
  range(ZellerG_2014_sample[,10])
  ZellerG_2014_co_index <- ZellerG_2014_sample[18.49<ZellerG_2014_sample[,10] &ZellerG_2014_sample[,10]<25,]
  table(ZellerG_2014_co_index[,4])

  ZellerG_2014_co_index <- ZellerG_2014_co_index[order(ZellerG_2014_co_index$study_condition),]
  ZellerG_2014_otu_g <- ZellerG_2014_otu[rownames(ZellerG_2014_co_index),]
  zerocol <- which(colSums(ZellerG_2014_otu_g)==0)
  if(length(zerocol) >0 ){
    ZellerG_2014_otu_g <- ZellerG_2014_otu_g[,-zerocol];
  }
  dim(ZellerG_2014_otu_g)

  data_c <- t(ZellerG_2014_otu_g)
  rownames(data_c) <- NULL
  samp <- ZellerG_2014_co_index
  phenotypeData = AnnotatedDataFrame(samp)
  ZellerG_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)

}
##FengQ_2015 dataset
{
  FengQ_2015 <- curatedMetagenomicData("FengQ_2015.metaphlan_bugs_list.stool",
                                       counts=TRUE, dryrun=FALSE, bugs.as.phyloseq=TRUE)
  FengQ_2015_genus <- phyloseq::tax_glom(FengQ_2015$FengQ_2015.metaphlan_bugs_list.stool, taxrank="Species")
  FengQ_2015_tax <- tax_table(FengQ_2015_genus )
  FengQ_2015_otu <- as.data.frame(t(otu_table(FengQ_2015_genus)))
  FengQ_2015_sample <- as.data.frame( sample_data(FengQ_2015_genus))

  FengQ_2015_sample <- FengQ_2015_sample[FengQ_2015_sample[,4]!="adenoma",]
  FengQ_2015_sample <- FengQ_2015_sample[complete.cases(FengQ_2015_sample[,10]),]
  range(FengQ_2015_sample[,10])
  FengQ_2015_co_index <- FengQ_2015_sample[18.49<FengQ_2015_sample[,10] &FengQ_2015_sample[,10]<25,]
  #FengQ_2015_co_index <- FengQ_2015_sample
  table(FengQ_2015_co_index[,4])

  FengQ_2015_co_index <- FengQ_2015_co_index[order(FengQ_2015_co_index$study_condition),]
  FengQ_2015_otu_g <- FengQ_2015_otu[rownames(FengQ_2015_co_index),]
  zerocol <- which(colSums(FengQ_2015_otu_g)==0)
  if(length(zerocol) >0 ){
    FengQ_2015_otu_g <- FengQ_2015_otu_g[,-zerocol];
  }
  dim(FengQ_2015_otu_g)

  data_c <- t(FengQ_2015_otu_g)
  rownames(data_c) <- NULL
  samp <- FengQ_2015_co_index
  phenotypeData = AnnotatedDataFrame(samp)
  FengQ_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)

}
##YuJ_2015 dataset
{
  YuJ_2015 <- curatedMetagenomicData("YuJ_2015.metaphlan_bugs_list.stool",
                                     counts=TRUE, dryrun=FALSE, bugs.as.phyloseq=TRUE)
  YuJ_2015_genus <- phyloseq::tax_glom(YuJ_2015$YuJ_2015.metaphlan_bugs_list.stool, taxrank="Species")
  YuJ_2015_tax <- tax_table(YuJ_2015_genus )
  YuJ_2015_otu <- as.data.frame(t(otu_table(YuJ_2015_genus)))
  YuJ_2015_sample <- as.data.frame( sample_data(YuJ_2015_genus))

  YuJ_2015_sample <- YuJ_2015_sample[YuJ_2015_sample[,4]!="adenoma",]
  YuJ_2015_sample <- YuJ_2015_sample[order(YuJ_2015_sample$study_condition),]
  YuJ_2015_otu_g <- YuJ_2015_otu[rownames(YuJ_2015_sample),]
  #YuJ_2015_otu_g <- as.matrix(YuJ_2015_otu[,sapply("g__", grepl,colnames(YuJ_2015_otu)) == TRUE])
  dim(YuJ_2015_otu_g)

  data_c <- t(YuJ_2015_otu_g)
  rownames(data_c) <- NULL
  samp <- YuJ_2015_sample
  phenotypeData = AnnotatedDataFrame(samp)
  YuJ_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)

}
##VogtmannE_2016 dataset
{
  VogtmannE_2016 <- curatedMetagenomicData("VogtmannE_2016.metaphlan_bugs_list.stool",
                                           counts=TRUE, dryrun=FALSE, bugs.as.phyloseq=TRUE)
  VogtmannE_2016_genus <- phyloseq::tax_glom(VogtmannE_2016$VogtmannE_2016.metaphlan_bugs_list.stool, taxrank="Species")
  VogtmannE_2016_tax <- tax_table(VogtmannE_2016_genus )
  VogtmannE_2016_otu <- as.data.frame(t(otu_table(VogtmannE_2016_genus)))
  VogtmannE_2016_sample <- as.data.frame( sample_data(VogtmannE_2016_genus))

  VogtmannE_2016_sample <- VogtmannE_2016_sample[complete.cases(VogtmannE_2016_sample[,9]),]
  range(VogtmannE_2016_sample[,9])
  VogtmannE_2016_co_index <- VogtmannE_2016_sample[18.49<VogtmannE_2016_sample[,9] &VogtmannE_2016_sample[,9]<25,]
  table(VogtmannE_2016_co_index[,4])

  VogtmannE_2016_co_index <- VogtmannE_2016_co_index[complete.cases(VogtmannE_2016_co_index[,4]),]
  VogtmannE_2016_co_index <- VogtmannE_2016_co_index[order(VogtmannE_2016_co_index$study_condition),]
  VogtmannE_2016_otu_g <- VogtmannE_2016_otu[rownames(VogtmannE_2016_co_index),]
  #VogtmannE_2016_otu_g <- as.matrix(VogtmannE_2016_otu[,sapply("g__", grepl,colnames(VogtmannE_2016_otu)) == TRUE])
  zerocol <- which(colSums(VogtmannE_2016_otu_g)==0)
  if(length(zerocol) >0 ){
    VogtmannE_2016_otu_g <- VogtmannE_2016_otu_g[,-zerocol];
  }
  dim(VogtmannE_2016_otu_g)

  data_c <- t(VogtmannE_2016_otu_g)
  rownames(data_c) <- NULL
  samp <- VogtmannE_2016_co_index
  phenotypeData = AnnotatedDataFrame(samp)
  VogtmannE_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)

}
{
  HanniganGD_2017 <- curatedMetagenomicData("HanniganGD_2017.metaphlan_bugs_list.stool",
                                            counts=TRUE, dryrun=FALSE, bugs.as.phyloseq=TRUE)
  HanniganGD_2017_genus <- phyloseq::tax_glom(HanniganGD_2017$HanniganGD_2017.metaphlan_bugs_list.stool, taxrank="Species")
  HanniganGD_2017_tax <- tax_table(HanniganGD_2017_genus )
  HanniganGD_2017_otu <- as.data.frame(t(otu_table(HanniganGD_2017_genus)))
  HanniganGD_2017_sample <- as.data.frame( sample_data(HanniganGD_2017_genus))

  HanniganGD_2017_sample <- HanniganGD_2017_sample[HanniganGD_2017_sample[,4]!="adenoma",]
  HanniganGD_2017_sample <- HanniganGD_2017_sample[order(HanniganGD_2017_sample$study_condition),]
  HanniganGD_2017_otu_g <- HanniganGD_2017_otu[rownames(HanniganGD_2017_sample),]
  #HanniganGD_2017_otu_g <- as.matrix(HanniganGD_2017_otu[,sapply("g__", grepl,colnames(HanniganGD_2017_otu)) == TRUE])
  zerocol <- which(colSums(HanniganGD_2017_otu_g)==0)
  if(length(zerocol) >0 ){
    HanniganGD_2017_otu_g <- HanniganGD_2017_otu_g[,-zerocol];
  }
  dim(HanniganGD_2017_otu_g)

  data_c <- t(HanniganGD_2017_otu_g)
  rownames(data_c) <- NULL
  samp <- HanniganGD_2017_sample
  phenotypeData = AnnotatedDataFrame(samp)
  HanniganGD_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)

}
{
  ThomasAM_2018a <- curatedMetagenomicData("ThomasAM_2018a.metaphlan_bugs_list.stool",
                                           counts=TRUE, dryrun=FALSE, bugs.as.phyloseq=TRUE)
  ThomasAM_2018a_genus <- phyloseq::tax_glom(ThomasAM_2018a$ThomasAM_2018a.metaphlan_bugs_list.stool, taxrank="Species")
  ThomasAM_2018a_tax <- tax_table(ThomasAM_2018a_genus )
  ThomasAM_2018a_otu <- as.data.frame(t(otu_table(ThomasAM_2018a_genus)))
  ThomasAM_2018a_sample <- as.data.frame( sample_data(ThomasAM_2018a_genus))

  ThomasAM_2018a_sample <- ThomasAM_2018a_sample[complete.cases(ThomasAM_2018a_sample[,10]),]
  range(ThomasAM_2018a_sample[,10])
  ThomasAM_2018a_co_index <- ThomasAM_2018a_sample[18.49<ThomasAM_2018a_sample[,10] &ThomasAM_2018a_sample[,10]<25,]
  table(ThomasAM_2018a_co_index[,4])

  ThomasAM_2018a_co_index <- ThomasAM_2018a_co_index[ThomasAM_2018a_co_index[,4]!="adenoma",]
  ThomasAM_2018a_co_index <- ThomasAM_2018a_co_index[order(ThomasAM_2018a_co_index$study_condition),]
  ThomasAM_2018a_otu_g <- ThomasAM_2018a_otu[rownames(ThomasAM_2018a_co_index),]
  #ThomasAM_2018a_otu_g <- as.matrix(ThomasAM_2018a_otu[,sapply("g__", grepl,colnames(ThomasAM_2018a_otu)) == TRUE])
  zerocol <- which(colSums(ThomasAM_2018a_otu_g)==0)
  if(length(zerocol) >0 ){
    ThomasAM_2018a_otu_g <- ThomasAM_2018a_otu_g[,-zerocol];
  }
  dim(ThomasAM_2018a_otu_g)

  data_c <- t(ThomasAM_2018a_otu_g)
  rownames(data_c) <- NULL
  samp <- ThomasAM_2018a_co_index
  phenotypeData = AnnotatedDataFrame(samp)
  ThomasAM_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)

}

X1 <- ZellerG_2014_otu_g
group1 <-   ZellerG_2014_co_index$study_condition
Z1 <- ifelse(group1=="CRC",1,0)
X2 <- FengQ_2015_otu_g
group2 <- FengQ_2015_co_index$study_condition
Z2 <- ifelse(group2=="CRC",1,0)
X3 <- YuJ_2015_otu_g
group3 <- YuJ_2015_sample$study_condition
Z3 <- ifelse(group3=="CRC",1,0)
X4 <- VogtmannE_2016_otu_g
group4 <- VogtmannE_2016_co_index$study_condition
Z4 <- ifelse(group4=="CRC",1,0)
X5 <- HanniganGD_2017_otu_g
group5 <- HanniganGD_2017_sample$study_condition
Z5 <- ifelse(group5=="CRC",1,0)
X6 <- ThomasAM_2018a_otu_g
group6 <- ThomasAM_2018a_co_index$study_condition
Z6 <- ifelse(group6=="CRC",1,0)

group <- list(group1,group2,group3,group4);
X <- list(as.matrix(X1),as.matrix(X2),as.matrix(X3),as.matrix(X4))
Z <- list(Z1,Z2,Z3,Z4)
X_MRexperiment <- list(ZellerG_MRexperiment,FengQ_MRexperiment,YuJ_MRexperiment,
                       VogtmannE_MRexperiment)

n <- length(X)

}

## run each method
{
ZILNM <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)
ZINB <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
ZIP <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
ZILNM_cov <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
ZINB_cov <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
ZIP_cov <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
ZR <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
SVT <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
pg <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
pca_out_x<- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
pca_out_y<- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);
pcoa <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)
tsne_out_y <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)
tsne_out_x <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)
gllvm_nbva<- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL);gllvm_nbva_cov<- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)
ZIFA<- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)
mbi <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)
sav <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)
dmn <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)
zipfa <- list(data1=NULL,data2=NULL,data3=NULL,data4=NULL)

deseq_p_adj <- edger_p_adj<- meta_p_adj <- meta_p_adj2 <- vector("list", n)
zip_p_adj1 <- zip_p_adj2<- zip_p_adj3 <- zip_p_adj4 <- zip_p_adj5 <- vector("list", n)
zinb_p_adj1 <- zinb_p_adj2<- zinb_p_adj3 <- zinb_p_adj4 <- zinb_p_adj5 <- vector("list", n)
zilnm_p_adj1 <- zilnm_p_adj2<- zilnm_p_adj3 <- zilnm_p_adj4 <- zilnm_p_adj5 <- vector("list", n)
mbimpute_p_adj1 <- mbimpute_p_adj2<- mbimpute_p_adj3 <- mbimpute_p_adj4 <- mbimpute_p_adj5 <- vector("list", n)
saver_p_adj1 <- saver_p_adj2<- saver_p_adj3 <- saver_p_adj4 <- saver_p_adj5 <- vector("list", n)

zip_p <- zip_p2<- zip_p3 <- zip_p4 <- vector("list", n)
zinb_p <- zinb_p2<- zinb_p3 <- zinb_p4 <- vector("list", n)
zilnm_p <- zilnm_p2<- zilnm_p3 <- zilnm_p4 <-  vector("list", n)
mbimpute_p <- mbimpute_p2<- mbimpute_p3 <- mbimpute_p4 <-  vector("list", n)
saver_p <- saver_p2<- saver_p3 <- saver_p4 <-  vector("list", n)

zip_p_adj_cov1 <- zip_p_adj_cov2<- zip_p_adj_cov3 <- zip_p_adj_cov4 <- zip_p_adj_cov5 <- vector("list", n)
zinb_p_adj_cov1 <- zinb_p_adj_cov2<- zinb_p_adj_cov3 <- zinb_p_adj_cov4 <- zinb_p_adj_cov5 <- vector("list", n)
zilnm_p_adj_cov1 <- zilnm_p_adj_cov2<- zilnm_p_adj_cov3 <- zilnm_p_adj_cov4 <- zilnm_p_adj_cov5 <- vector("list", n)

zip_p_cov <- zip_p_cov2<- zip_p_cov3 <- zip_p_cov4 <- vector("list", n)
zinb_p_cov <- zinb_p_cov2<- zinb_p_cov3 <- zinb_p_cov4 <- vector("list", n)
zilnm_p_cov <- zilnm_p_cov2<- zilnm_p_cov3 <- zilnm_p_cov4 <-  vector("list", n)

t_test_p <- t_test_p_adj <- vector("list",n)

for(i in 1:n){
  ZILNM[[i]] <- tryCatch({ZIPPCAlnm(X[[i]],rank = T)},
                         error=function(e){ NULL})
  ZINB[[i]] <- tryCatch({ZIPPCApn(X[[i]],rank = T)},
                        error=function(e){ NULL})
  ZIP[[i]] <- tryCatch({ZIPPCApn(X[[i]],family = "poisson",rank = T)},
                       error=function(e){ NULL})
  ZILNM_cov[[i]] <- tryCatch({ZIPPCAlnm(X[[i]],Z[[i]],rank = T)},
                             error=function(e){ NULL})
  ZINB_cov[[i]] <- tryCatch({ZIPPCApn(X[[i]],Z[[i]],rank = T)},
                            error=function(e){ NULL})
  ZIP_cov[[i]] <- tryCatch({ZIPPCApn(X[[i]],Z[[i]],family = "poisson",rank = T)},
                           error=function(e){ NULL})
  ZR[[i]] <- zr(W =X[[i]] , alpha = 0.5) # replaced by 0.5
  SVT[[i]] <- tryCatch({svt(W = X[[i]], r = 20, alpha = 0.5)},
                       error=function(e){NULL})
  pg[[i]] <- tryCatch({autoTuneProxGradient(W = X[[i]])$X_hat},
                      error=function(e){ NULL})
  dmn[[i]] <- tryCatch({DirichletMultinomial::dmn(X[[i]],k=10)},
                  error=function(e){NaN})
  zipfa[[i]] <- tryCatch({ZIPFA(X[[i]], k = 2, display = F,cut=1)},
                    error=function(e){NaN})

  pca_out_x[[i]] <- prcomp(X[[i]])
  pca_out_y[[i]] <- prcomp(log2(1+X[[i]]))
  dist_bray <- vegan::vegdist(log2(1+X[[i]]), method="bray")
  pcoa[[i]] <- ape::pcoa(dist_bray)
  tsne_out_y[[i]] <- Rtsne::Rtsne(log2(1+X[[i]]),initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))
  tsne_out_x[[i]] <- Rtsne::Rtsne(X[[i]],initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))

  gllvm_nbva[[i]] <- tryCatch({gllvm::gllvm(X[[i]], family = "negative.binomial",method = "VA",Lambda.struc = "diagonal",row.eff= "fixed")},
                              error=function(e){NULL})
  gllvm_nbva_cov[[i]] <-tryCatch({gllvm::gllvm(X[[i]],as.data.frame(Z[[i]]),family = "negative.binomial",method = "VA",Lambda.struc = "diagonal",row.eff= "fixed")},
                                 error=function(e){NULL})

  source_python('/lustre/home/acct-clswt/clswt/project/zyy/empirical/ZIFA_copy.py')
  ZIFA[[i]] <- tryCatch({fitModel(log2(1+X[[i]]),2)},error=function(e){NULL})

  ##mbimpute
  mbi[[i]] <- tryCatch({mbImpute(condition = Z[[i]], otu_tab = X[[i]],parallel = TRUE, ncores = 4)},
                       error=function(e){ NULL})

  ##saver
  sav[[i]] <- tryCatch({saver(t(X[[i]]), ncores = 12)},
                       error=function(e){ NULL})

  ##DESeq2
  group_gile <- data.frame(cbind(rownames(c(1:nrow(X[[i]]))),Z[[i]]))
  colnames(group_gile) <- c("X_cov")
  dds <- tryCatch({DESeqDataSetFromMatrix(countData = t(X[[i]]+1), colData = group_gile, design = ~X_cov)},
                  error=function(e){ NULL})
  if(is.null(dds)){
    deseq_p_adj[[i]] <- NaN;
  }else{
    dds <- DESeq(dds,quiet=T)
    res <- results(dds)#提取分析结果
    deseq_p_adj[[i]] <-res$padj}

  ##edgeR
  X_tr <- data.frame(t(X[[i]]))
  y<-DGEList(counts=X_tr,genes=rownames(X_tr))
  y<-edgeR::calcNormFactors(y)
  design <- model.matrix(~Z[[i]])
  Disp <- estimateDisp(y, design, robust = TRUE)
  fit <- glmFit(Disp, design)
  qlf <- glmLRT(fit, coef=2)
  edger_p_adj[[i]] <- p.adjust(qlf$table$PValue, method = "BH")

  ##metagenomeSeq:fitFeatureModel
  XX <- cumNorm(X_MRexperiment[[i]], p = 0.5)
  pd <- pData(XX)
  mod <- model.matrix(~1+Z[[i]], data = pd)
  fit <- tryCatch({fitFeatureModel(XX, mod)},
                  error=function(e){ NULL})
  if(is.null(fit)){
    meta_p_adj[[i]] <- NaN;
  }else{
    meta <- fit@pvalues
    meta_p_adj[[i]] <- p.adjust(meta,method = "fdr")
  }

}

}

## save results
save.image("real data/CRC.RData")

