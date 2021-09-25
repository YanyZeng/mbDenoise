setwd("~/mbDenoise/reproducibility")

##load packages/methods source
library(readxl)
library(openxlsx)
library(curatedMetagenomicData)
library("phyloseq")
library(ade4)
library(edgeR)
library(DESeq2)
library(ape)
library(metagenomeSeq)
library(coin)
library(reticulate)
library(mbImpute)
library(glmnet)
library(SAVER)
library(stringr)
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

##Figures 8b-c and 9 data processing
{
pool_index <- read.xlsx("real data/oral.xlsx",sheet=1,startRow = 2,rowNames = T)
pool_meta <- read.xlsx("real data/meta.xlsx",sheet=1,startRow = 2,rowNames = T)
Subp_meta <- pool_meta[pool_meta$Sample.Region!="Tongue",];

##pool
Subp <- (pool_index[,1:72]);
Subp_tax <- (pool_index$Consensus.Lineage %>% str_split(";",  simplify = TRUE))
colnames(Subp_tax) <- c( "Kindom","Phylum", "Class", "Order", "Family", "Genus","Species")
rownames(Subp_tax) <- rownames(Subp)
OTU = otu_table(Subp, taxa_are_rows = TRUE)
TAX = tax_table(Subp_tax)
SAMPLE = sample_data(Subp_meta)
Subp_physeq = phyloseq(OTU, TAX,SAMPLE)
Subp.genus <- phyloseq::tax_glom(Subp_physeq, taxrank="Species")
Subp.genus_otu <- as.data.frame(t(otu_table(Subp.genus)))
Subp.genus_sample <- as.data.frame(sample_data(Subp.genus))
Subp.genus_tax <- tax_table(Subp.genus)
colnames(Subp.genus_otu) <- Subp.genus_tax[,7]

X1 <- X2 <- Subp.genus_otu
Z1 <- Subp.genus_sample$`Healthy?`
group1 <- ifelse(Z1==1,"Healthy","Disease")
group2 <- Subp.genus_sample$Sample.Region
Z2<- ifelse(group2=="Supra",1,0)

data_c <- t(Subp.genus_otu)
rownames(data_c) <- NULL
samp <- Subp_meta
phenotypeData = AnnotatedDataFrame(samp)
Subp_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)


group <- list(group1,group2);
X <- list(as.matrix(X1),as.matrix(X2))
Z <- list(Z1,Z2)
X_MRexperiment <- list(Subp_MRexperiment,Subp_MRexperiment)
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

ZINB_pcoa <- ZINB_pca <- ZINB_tsne <- ZIP_pcoa <- ZIP_pca <- ZIP_tsne <-ZILNM_pcoa <- ZILNM_pca <- ZILNM_tsne <- vector("list", n)

t_test_p <- t_test_p_adj <- vector("list", n)

for(i in 1:1){
  ZILNM[[i]] <- tryCatch({ZIPPCAlnm(X[[i]],rank = T)},
                         error=function(e){ NULL})
  ZINB[[i]] <- tryCatch({ZIPPCApn(X[[i]],rank = T)},
                        error=function(e){ NULL})
  ZIP[[i]] <- tryCatch({ZIPPCApn(X[[i]],family = "poisson",rank = T)},
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

  ZIFA[[i]] <- tryCatch({fitModel(log2(1+X[[i]]),2)},error=function(e){NULL})

  ##saver
  sav[[i]] <- tryCatch({saver(t(X[[i]]), ncores = 12)},
                       error=function(e){ NULL})



}
for(i in 2:2){
  ZILNM[[i]] <- ZILNM[[1]]
  ZINB[[i]] <- ZINB[[1]]
  ZIP[[i]] <- ZIP[[1]]
  dmn[[i]] <- dmn[[1]]
  zipfa[[i]] <- zipfa[[1]]
  ZR[[i]] <- ZR[[1]]
  SVT[[i]] <- SVT[[1]]
  pg[[i]] <- pg[[1]]
  pca_out_x[[i]] <- pca_out_x[[1]]
  pca_out_y[[i]] <- pca_out_y[[1]]
  pcoa[[i]] <- pcoa[[1]]
  tsne_out_y[[i]] <- tsne_out_y[[1]]
  tsne_out_x[[i]] <- tsne_out_x[[1]]
  gllvm_nbva[[i]] <- gllvm_nbva[[1]]
  ZIFA[[i]] <- ZIFA[[1]]

  ##saver
  sav[[i]] <-sav[[1]]
}
for(i in 1:n){

  ZILNM_cov[[i]] <- tryCatch({ZIPPCAlnm(X[[i]],Z[[i]],rank = T)},
                             error=function(e){ NULL})
  ZINB_cov[[i]] <- tryCatch({ZIPPCApn(X[[i]],Z[[i]],rank = T)},
                            error=function(e){ NULL})
  ZIP_cov[[i]] <- tryCatch({ZIPPCApn(X[[i]],Z[[i]],family = "poisson",rank = T)},
                           error=function(e){ NULL})
  gllvm_nbva_cov[[i]] <-tryCatch({gllvm::gllvm(X[[i]],as.data.frame(Z[[i]]),family = "negative.binomial",method = "VA",Lambda.struc = "diagonal",row.eff= "fixed")},
                                 error=function(e){NULL})

  ##mbimpute
  mbi[[i]] <- tryCatch({mbImpute(condition = Z[[i]], otu_tab = X[[i]],parallel = TRUE, ncores = 4)},
                       error=function(e){ NULL})

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
save.image("real data/CP2.RData")




