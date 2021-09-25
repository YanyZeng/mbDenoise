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

##Figure 6 data processing
{
  DhakanDB <- curatedMetagenomicData("DhakanDB_2019.metaphlan_bugs_list.stool",
                                                 counts=TRUE, dryrun=FALSE, bugs.as.phyloseq=TRUE)
  DhakanDB_species <- phyloseq::tax_glom(DhakanDB$DhakanDB_2019.metaphlan_bugs_list.stool, taxrank="Species")
  DhakanDB_species_tax <- tax_table(DhakanDB_species )
  DhakanDB_species_otu <- as.data.frame(t(otu_table(DhakanDB_species)))
  DhakanDB_species_sample <- as.data.frame( sample_data(DhakanDB_species))
  DhakanDB_species_sample <- DhakanDB_species_sample[18.49<DhakanDB_species_sample$BMI &DhakanDB_species_sample$BMI<25,]
  DhakanDB_species_sample <- DhakanDB_species_sample[complete.cases(DhakanDB_species_sample$BMI),]
  DhakanDB_species_sample <- DhakanDB_species_sample[order(DhakanDB_species_sample$location),]
  DhakanDB_species_otu <- DhakanDB_species_otu[rownames(DhakanDB_species_sample),]
  
  zerocol <- which(colSums(DhakanDB_species_otu)==0)
  if(length(zerocol) >0 ){
    DhakanDB_species_otu <- DhakanDB_species_otu[,-zerocol];
  }
  dim(DhakanDB_species_otu)
  
  data_c <- t(DhakanDB_species_otu)
  rownames(data_c) <- NULL
  samp <- DhakanDB_species_sample
  phenotypeData = AnnotatedDataFrame(samp)
  DhakanDB_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)
  
  X1 <- DhakanDB_species_otu
  group1 <-   DhakanDB_species_sample$location
  Z1 <- ifelse(group1=="Bhopal",1,0)
  
  DhakanDB_species_sample2 <- DhakanDB_species_sample[DhakanDB_species_sample$location=="Bhopal",]
  DhakanDB_species_otu2 <- DhakanDB_species_otu[rownames(DhakanDB_species_sample2),]
  
  zerocol <- which(colSums(DhakanDB_species_otu2)==0)
  if(length(zerocol) >0 ){
    DhakanDB_species_otu2 <- DhakanDB_species_otu2[,-zerocol];
  }
  dim(DhakanDB_species_otu2)
  
  data_c <- t(DhakanDB_species_otu2)
  rownames(data_c) <- NULL
  samp <- DhakanDB_species_sample2
  phenotypeData = AnnotatedDataFrame(samp)
  DhakanDB_MRexperiment2 <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)
  
  X2 <- DhakanDB_species_otu2
  Z2 <- c(rep(0,round(length(DhakanDB_species_sample2$location)/2)),rep(1,length(DhakanDB_species_sample2$location)-round(length(DhakanDB_species_sample2$location)/2)))
  group2 <- ifelse(Z2==0,"Bhopal","Kerala")
  

group <- list(group1,group2);
X <- list(as.matrix(X1),as.matrix(X2))
Z <- list(Z1,Z2)
X_MRexperiment <- list(DhakanDB_MRexperiment,DhakanDB_MRexperiment2)
n <- length(X)

}

## run each method
{
ZILNM <- list(data1=NULL,data2=NULL)
ZINB <- list(data1=NULL,data2=NULL);
ZIP <- list(data1=NULL,data2=NULL);
ZILNM_cov <- list(data1=NULL,data2=NULL);
ZINB_cov <- list(data1=NULL,data2=NULL);
ZIP_cov <- list(data1=NULL,data2=NULL);
ZR <- list(data1=NULL,data2=NULL);
SVT <- list(data1=NULL,data2=NULL);
pg <- list(data1=NULL,data2=NULL);
pca_out_x<- list(data1=NULL,data2=NULL);
pca_out_y<- list(data1=NULL,data2=NULL);
pcoa <- list(data1=NULL,data2=NULL)
tsne_out_y <- list(data1=NULL,data2=NULL)
tsne_out_x <- list(data1=NULL,data2=NULL)
gllvm_nbva<- list(data1=NULL,data2=NULL);
gllvm_nbva_cov<- list(data1=NULL,data2=NULL)
ZIFA<- list(data1=NULL,data2=NULL)
mbi <- list(data1=NULL,data2=NULL)
sav <- list(data1=NULL,data2=NULL)
dmn <- list(data1=NULL,data2=NULL)
zipfa <- list(data1=NULL,data2=NULL)


deseq_p_adj <- edger_p_adj<- meta_p_adj <- meta_p_adj2 <- vector("list", n)
zip_p_adj1 <- zip_p_adj2<- zip_p_adj3 <- zip_p_adj4 <- zip_p_adj5 <- vector("list", n)
zinb_p_adj1 <- zinb_p_adj2<- zinb_p_adj3 <- zinb_p_adj4 <- zinb_p_adj5 <- vector("list", n)
zilnm_p_adj1 <- zilnm_p_adj2<- zilnm_p_adj3 <- zilnm_p_adj4 <- zilnm_p_adj5 <- vector("list", n)

zip_p <- zip_p2<- zip_p3 <- zip_p4 <- vector("list", n)
zinb_p <- zinb_p2<- zinb_p3 <- zinb_p4 <- vector("list", n)
zilnm_p <- zilnm_p2<- zilnm_p3 <- zilnm_p4 <-  vector("list", n)

zip_p_adj_cov1 <- zip_p_adj_cov2<- zip_p_adj_cov3 <- zip_p_adj_cov4 <- zip_p_adj_cov5 <-  vector("list", n)
zinb_p_adj_cov1 <- zinb_p_adj_cov2<- zinb_p_adj_cov3 <- zinb_p_adj_cov4 <- zinb_p_adj_cov5 <-  vector("list", n)
zilnm_p_adj_cov1 <- zilnm_p_adj_cov2<- zilnm_p_adj_cov3 <- zilnm_p_adj_cov4 <- zilnm_p_adj_cov5 <-  vector("list", n)
mbimpute_p_adj1 <- mbimpute_p_adj2<- mbimpute_p_adj3 <- mbimpute_p_adj4 <- mbimpute_p_adj5 <- vector("list", n)
saver_p_adj1 <- saver_p_adj2<- saver_p_adj3 <- saver_p_adj4 <- saver_p_adj5 <- vector("list", n)

zip_p_cov <- zip_p_cov2<- zip_p_cov3 <- zip_p_cov4 <- vector("list", n)
zinb_p_cov <- zinb_p_cov2<- zinb_p_cov3 <- zinb_p_cov4 <- vector("list", n)
zilnm_p_cov <- zilnm_p_cov2<- zilnm_p_cov3 <- zilnm_p_cov4 <-  vector("list", n)
mbimpute_p <- mbimpute_p2<- mbimpute_p3 <- mbimpute_p4 <-  vector("list", n)
saver_p <- saver_p2<- saver_p3 <- saver_p4 <-  vector("list", n)
t_test_p <- t_test_p_adj <- vector("list",n)

for(i in 1:n){

  ##mbimpute
  mbi[[i]] <- tryCatch({mbImpute(condition = Z[[i]], otu_tab = X[[i]],parallel = TRUE, ncores = 4)},
                       error=function(e){ NULL})
 
  ##saver
  sav[[i]] <- tryCatch({saver(t(X[[i]]), ncores = 12)},
                       error=function(e){ NULL})

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



  ##DESeq2
  group_file <- data.frame(cbind(rownames(c(1:nrow(X[[i]]))),Z[[i]]))
  colnames(group_file) <- c("X_cov")
  dds <- tryCatch({DESeqDataSetFromMatrix(countData = t(X[[i]]+1), colData = group_file, design = ~X_cov)},
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
save.image("real data/India_diet.RData")

