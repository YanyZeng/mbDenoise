# mbDenoise
Microbiome data denoising framework (mbDenoise) with zero-inflated probabilistic PCA, which can be used for downstream high-level statistical analysis, including ordination, compositional normalization, differential abundance analysis, etc.

# Installation
You can also install the released version of mbDenoise from github with:

``` r
install.packages("devtools")  
devtools::install_github("YanyZeng/mbDenoise")  
library(mbDenoise)
```

# Description
Marker gene and metagenomic sequencing studies have illustrated the importance of microbial communities for human and environmental health. However, the analysis of microbiome data has several caveats and technical challenges. One of the major issues is that count matrices contain a large proportion of zeros, some of which are biological zeros, whereas others are technical zeros. In addition, the measurements suffer from unequal sequencing depth, overdispersion, and data redundancy. These nuisance factors introduce substantial noise, which distorts the biological signal and hinders downstream analyses. To address these challenges, we propose mbDenoise, an accurate and robust method for denoising microbiome data. mbDenoise assumes a zero-inflated probabilistic PCA (ZIPPCA) model consisting of a point mass at zero and a negative binomial component. It uses variational approximation to learn the latent structure, borrowing information across samples and taxa, and then recovers the true underlying abundance levels using the posterior mean. We evaluate the performance of mbDenoise and compare it to state-of-the-art methods, using both simulated and real datasets, in terms of unconstrained ordination, alpha and beta diversity analysis, and differential abundance testing. mbDenoise outperforms existing methods to extract the signal for high-level analyses. 

# Usage
mbDenoise is based on zero-inflated probabilistic PCA with logistical normal multinomial (ZIPPCA-LNM),
Poisson (ZIPPCA-Poi) and negative-binomial model (ZIPPCA-NB). And mbDenoise with ZIPPCA-NB model is recommended for empirical data analysis.

```r
ZIPPCApn(X,V = NULL,family = "negative.binomial",n.factors = 2,rank = FALSE,trace = FALSE,
          maxit = 100,parallel = TRUE)

ZIPPCAlnm(X,V = NULL,n.factors = 2,rank = FALSE,trace = FALSE,maxit = 100,parallel = TRUE)
```
* X: matrix of observations.
* V: vector of the sample covariate.
* family: distribution of models. Two options are "poisson" and "negative.binomial". Defaults to "negative.binomial".
* n.factors: the rank or number of factors, after dimensional reduction. Defaults to 2.
* rank: logical, if TRUE, the rank or number of factors, is chosen from 1 to 5 by HIC (hybrid information criterion) in ZIPPCApn function and BIC (Bayesian information criterion) in ZIPPCAlnm function. Defaults to FALSE.
* trace: logical, defaults to FALSE. if TRUE each current iteration step information will be printed.
* maxit: maximum number of iterations within optim and constrOptim function, defaults to 100.
* parallel: logical, if TRUE, use parallel toolbox to accelerate.

# Example
We use a microbiome dataset of Dhakan et al (2019) as a basic example which shows you how to solve a common problem and demonstrates the use of mbDenoise.

``` r
## Data preparing
  DhakanDB <- curatedMetagenomicData::curatedMetagenomicData("DhakanDB_2019.metaphlan_bugs_list.stool",
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

  X2 <- DhakanDB_species_otu2
  Z2 <- c(rep(0,round(length(DhakanDB_species_sample2$location)/2)),rep(1,length(DhakanDB_species_sample2$location)-round(length(DhakanDB_species_sample2$location)/2)))
  group2 <- ifelse(Z2==0,"Bhopal","Kerala")
  
  group <- list(group1,group2);
  X <- list(as.matrix(X1),as.matrix(X2))
  Z <- list(Z1,Z2)
  n <- length(X)
  
## Fitting models  
  ZINB <- list(data1=NULL,data2=NULL); ZINB_cov <- list(data1=NULL,data2=NULL)
  for(i in 1:n){
  ZINB[[i]] <- tryCatch({ZIPPCApn(X[[i]],rank = T)},
                        error=function(e){ NULL})
  ZINB_cov[[i]] <- tryCatch({ZIPPCApn(X[[i]],Z[[i]],rank = T)},
                            error=function(e){ NULL})
}

## Ordination analysis (beta diversity)
#mbDenoise-zinb
   adonis_result_dis = adonis2(dist(ZINB[[i]]$lvs$factor_scores2)~group[[i]],method = "euclidean")
    p1 <- ggplot(data.frame(ZINB[[i]]$lvs$factor_scores2),aes(x=ZINB[[i]]$lvs$factor_scores2[,1], y=ZINB[[i]]$lvs$factor_scores2[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_brewer(palette = "Set1")+
      xlab("mbDenoise-zinb F1")+ylab("mbDenoise-zinb F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    X_new <- log2(ZINB[[i]]$muz+1)
    ZINB_pca[[i]] <- prcomp(X_new)
    dist_bray <- vegan::vegdist(X_new, method="bray")
    ZINB_pcoa[[i]] <- ape::pcoa(dist_bray)
    set.seed(4)
    ZINB_tsne[[i]] <- Rtsne::Rtsne(X_new,initial_dims = ncol(X_new),perplexity=round((nrow(X_new)-2)/3))
    
#mbDenoise-zinb_pca   
    adonis_result_dis = adonis2(dist(ZINB_pca[[i]]$x[,1:2])~group[[i]],method = "euclidean")
    p2 <- ggplot(data.frame(ZINB_pca[[i]]$x[,1:2]),aes(ZINB_pca[[i]]$x[,1], y=ZINB_pca[[i]]$x[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_brewer(palette = "Set1")+
      xlab("mbDenoise-zinb_pca F1")+ylab("mbDenoise-zinb_pca F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
#mbDenoise-zinb_pcoa
    adonis_result_dis = adonis2(dist(ZINB_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    p3 <- ggplot(data.frame(ZINB_pcoa[[i]]$vectors[,1:2]),aes(ZINB_pcoa[[i]]$vectors[,1], y=ZINB_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_brewer(palette = "Set1")+
      xlab("mbDenoise-zinb_pcoa F1")+ylab("mbDenoise-zinb_pcoa F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
#mbDenoise-zinb_tsne
    adonis_result_dis = adonis2(dist(ZINB_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    p4 <- ggplot(data.frame(ZINB_tsne[[i]]$Y),aes(ZINB_tsne[[i]]$Y[,1], y=ZINB_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_brewer(palette = "Set1")+
      xlab("mbDenoise-zinb_tsne F1")+ylab("mbDenoise-zinb_tsne F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))

## Compositional normalization analysis (alpha diversity)

   diversity <- function(X){
    # shannon index
    sh_ele <- -X * log(X)
    #sh_ele[is.nan(sh_ele)] <- 0
    sh <- rowSums(sh_ele)
    
    # simpson index
    sp  <- 1-rowSums(X * X)
    
    return(list(Shannon  = sh, Simpson  = sp))
  }

  zinb_div <-  list()
  
  for(i in 1:n){
    #mbDenoise-zinb
    zinb_div[[i]] <- diversity(ZINB[[i]]$Q)
  }  

## Differential abundance analysis
  zinb <- zinb_p <- zinb_p_adj <- list()
  for(i in 1:n){ 
  if(is.null(ZINB_cov[[i]])){
    zinb_p[[i]] <- NaN;
  }else{
    X_new <- log2(ZINB_cov[[i]]$muz+1)

    for(l in 1:ncol(X[[i]])){
      zinb_p[[i]][l] <- tryCatch({t.test(X_new[,l]~Z[[i]])$p.value},error=function(e){ NaN})
    }}
  
  zinb_p_adj[[i]] <- p.adjust(zinb_p[[i]],method = "BH")
  zinb[[i]] <- na.omit(colnames(X[[i]])[(zinb_p_adj[[i]]) <0.05])

}

```
Reference:
Dhakan, D. et al. The unique composition of indian gut microbiome, gene catalogue, and
associated fecal metabolome deciphered using multi-omics approaches. Gigascience 8, giz004
(2019).
