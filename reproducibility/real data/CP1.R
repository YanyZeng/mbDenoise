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

##Figure 8a data processing
{
pool_index <- read.xlsx("real data/oral.xlsx",sheet=1,startRow = 2,rowNames = T)
Subg_index <- read.xlsx("real data/oral.xlsx",sheet=2,startRow = 2,rowNames = T)
Supra_index <- read.xlsx("real data/oral.xlsx",sheet=3,startRow = 2,rowNames = T)
Tongue_index <- read.xlsx("real data/oral.xlsx",sheet=4,startRow = 2,rowNames = T)
pool_meta <- read.xlsx("real data/meta.xlsx",sheet=1,startRow = 2,rowNames = T)
Subg_meta <- pool_meta[pool_meta$Sample.Region=="Subg",];rownames(Subg_meta) <- Subg_meta$SubG.SampleID
Supra_meta <- pool_meta[pool_meta$Sample.Region=="Supra",];rownames(Supra_meta) <- Supra_meta$Supra.SampleID
Tongue_meta <- pool_meta[pool_meta$Sample.Region=="Tongue",];rownames(Tongue_meta) <- Tongue_meta$Tongue.SampleID

##pool
pool <- (pool_index[,1:72]);
pool_tax <- (pool_index$Consensus.Lineage %>% str_split(";",  simplify = TRUE))
colnames(pool_tax) <- c( "Kindom","Phylum", "Class", "Order", "Family", "Genus","Species")
rownames(pool_tax) <- rownames(pool)
OTU = otu_table(pool, taxa_are_rows = TRUE)
TAX = tax_table(pool_tax)
SAMPLE = sample_data(pool_meta)
pool_physeq = phyloseq(OTU, TAX,SAMPLE)
pool.genus <- phyloseq::tax_glom(pool_physeq, taxrank="Species")
pool.genus_otu <- as.data.frame(t(otu_table(pool.genus)))
pool.genus_sample <- as.data.frame(sample_data(pool.genus))
pool.genus_tax <- tax_table(pool.genus)
colnames(pool.genus_otu) <- pool.genus_tax[,7]

X1 <- pool.genus_otu
Z1 <- pool.genus_sample$`Healthy?`
group1 <- ifelse(Z1==1,"Healthy","Disease")
pool_site <- pool_meta$Sample.Region
group11 <- pool_site
Z11 <- as.factor(group11)
data_c <- t(pool.genus_otu)
rownames(data_c) <- NULL
samp <- pool_meta
phenotypeData = AnnotatedDataFrame(samp)
pool_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)

##Subg
Subg <- (Subg_index[,1:24]);
Subg_tax <- (Subg_index$Consensus.Lineage %>% str_split(";",  simplify = TRUE))
colnames(Subg_tax) <- c( "Kindom","Phylum", "Class", "Order", "Family", "Genus","Species")
rownames(Subg_tax) <- rownames(Subg)
OTU = otu_table(Subg, taxa_are_rows = TRUE)
TAX = tax_table(Subg_tax)
SAMPLE = sample_data(Subg_meta)
Subg_physeq = phyloseq(OTU, TAX,SAMPLE)
Subg.genus <- phyloseq::tax_glom(Subg_physeq, taxrank="Species")
Subg.genus_otu <- as.data.frame(t(otu_table(Subg.genus)))
Subg.genus_sample <- as.data.frame(sample_data(Subg.genus))
Subg.genus_tax <- tax_table(Subg.genus)
colnames(Subg.genus_otu) <- Subg.genus_tax[,7]

X2 <- Subg.genus_otu
Z2 <- Subg.genus_sample$`Healthy?`
group2 <- ifelse(Z2==1,"Healthy","Disease")
data_c <- t(Subg.genus_otu)
rownames(data_c) <- NULL
samp <- Subg_meta
phenotypeData = AnnotatedDataFrame(samp)
Subg_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)


##Supra
Supra <- (Supra_index[,1:24]);
Supra_tax <- (Supra_index$Consensus.Lineage %>% str_split(";",  simplify = TRUE))
colnames(Supra_tax) <- c( "Kindom","Phylum", "Class", "Order", "Family", "Genus","Species")
rownames(Supra_tax) <- rownames(Supra)
OTU = otu_table(Supra, taxa_are_rows = TRUE)
TAX = tax_table(Supra_tax)
SAMPLE = sample_data(Supra_meta)
Supra_physeq = phyloseq(OTU, TAX,SAMPLE)
Supra.genus <- phyloseq::tax_glom(Supra_physeq, taxrank="Species")
Supra.genus_otu <- as.data.frame(t(otu_table(Supra.genus)))
Supra.genus_sample <- as.data.frame(sample_data(Supra.genus))
Supra.genus_tax <- tax_table(Supra.genus)
colnames(Supra.genus_otu) <- Supra.genus_tax[,7]

X3 <- Supra.genus_otu
Z3 <- Supra.genus_sample$`Healthy?`
group3 <- ifelse(Z3==1,"Healthy","Disease")
data_c <- t(Supra.genus_otu)
rownames(data_c) <- NULL
samp <- Supra_meta
phenotypeData = AnnotatedDataFrame(samp)
Supra_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)


##Tongue
Tongue <- (Tongue_index[,1:24]);
Tongue_tax <- (Tongue_index$Consensus.Lineage %>% str_split(";",  simplify = TRUE))
colnames(Tongue_tax) <- c( "Kindom","Phylum", "Class", "Order", "Family", "Genus","Species")
rownames(Tongue_tax) <- rownames(Tongue)
OTU = otu_table(Tongue, taxa_are_rows = TRUE)
TAX = tax_table(Tongue_tax)
SAMPLE = sample_data(Tongue_meta)
Tongue_physeq = phyloseq(OTU, TAX,SAMPLE)
Tongue.genus <- phyloseq::tax_glom(Tongue_physeq, taxrank="Species")
Tongue.genus_otu <- as.data.frame(t(otu_table(Tongue.genus)))
Tongue.genus_sample <- as.data.frame(sample_data(Tongue.genus))
Tongue.genus_tax <- tax_table(Tongue.genus)
colnames(Tongue.genus_otu) <- Tongue.genus_tax[,7]

X4 <- Tongue.genus_otu
Z4 <- Tongue.genus_sample$`Healthy?`
group4 <- ifelse(Z4==1,"Healthy","Disease")
data_c <- t(Tongue.genus_otu)
rownames(data_c) <- NULL
samp <- Tongue_meta
phenotypeData = AnnotatedDataFrame(samp)
Tongue_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)

group <- list(group1,group2,group3,group4);
X <- list(as.matrix(X1),as.matrix(X2),as.matrix(X3),as.matrix(X4))
Z <- list(Z1,Z2,Z3,Z4)
X_MRexperiment <- list(pool_MRexperiment,Subg_MRexperiment,Supra_MRexperiment,Tongue_MRexperiment)

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

ZINB_pcoa <- ZINB_pca <- ZINB_tsne <- ZIP_pcoa <- ZIP_pca <- ZIP_tsne <-ZILNM_pcoa <- ZILNM_pca <- ZILNM_tsne <- vector("list", n)

t_test_p <- t_test_p_adj <- vector("list", n)

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

  ZIFA[[i]] <- tryCatch({fitModel(log2(1+X[[i]]),2)},error=function(e){NULL})

  ##mbimpute
  mbi[[i]] <- tryCatch({mbImpute(condition = Z[[i]], otu_tab = X[[i]],parallel = TRUE, ncores = 4)},
                       error=function(e){ NULL})

  ##saver
  sav[[i]] <- tryCatch({saver(t(X[[i]]), ncores = 12)},
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
save.image("real data/CP1.RData")




