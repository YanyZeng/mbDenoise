setwd("~/mbDenoise/reproducibility")

library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(reshape)
library(vegan)
library(VennDiagram)
library(UpSetR)
library(gridExtra)

load("real data/India_diet.RData")

## Figure 6a: ordination and beta diversity
{
  ZINB_pcoa <- ZINB_pca <- ZINB_tsne <- ZIP_pcoa <- ZIP_pca <- ZIP_tsne <-list()
  ppca_pcoa <- ppca_pca <- ppca_tsne <- sav_pcoa <- sav_pca <-sav_tsne <- list()
  grid_all <- list()
  for(i in 1:1){ 
    set.seed(4)
    tsne_out_y[[i]] <- Rtsne::Rtsne(log2(1+X[[i]]),initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))
    tsne_out_x[[i]] <- Rtsne::Rtsne(X[[i]],initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))

    adonis_result_dis = adonis2(dist(pca_out_y[[i]]$x[,1:2])~group[[i]],method = "euclidean")
    pp2 <- ggplot(data.frame(pca_out_y[[i]]$x[,1:2]),aes(pca_out_y[[i]]$x[,1], y=pca_out_y[[i]]$x[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PCA F1")+ylab("PCA F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp3 <- ggplot(data.frame(pcoa[[i]]$vectors[,1:2]),aes(pcoa[[i]]$vectors[,1], y=pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PCoA F1")+ylab("PCoA F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
 
    adonis_result_dis = adonis2(dist(tsne_out_y[[i]]$Y)~group[[i]],method = "euclidean")
    pp5 <- ggplot(data.frame(tsne_out_y[[i]]$Y),aes(tsne_out_y[[i]]$Y[,1], y=tsne_out_y[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("t-SNE F1")+ylab("t-SNE F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = tryCatch({adonis2(dist(ZIFA[[i]][[1]])~group[[i]],method = "euclidean")},error=function(e){ NaN})
    pp6 <- tryCatch({ggplot(data.frame(ZIFA[[i]][[1]]),aes(x=ZIFA[[i]][[1]][,1], y=ZIFA[[i]][[1]][,2],colour=as.factor(group[[i]]))) +
        geom_point(size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
        scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
        xlab("ZIFA F1")+ylab("ZIFA F2")+labs(colour="Location") +
        theme_bw()+
        stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
        annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))},error=function(e){ NaN})
    
    adonis_result_dis = adonis2(dist(gllvm_nbva[[i]]$lvs)~group[[i]],method = "euclidean")
    pp7 <- ggplot(data.frame( gllvm_nbva[[i]]$lvs),aes(x= gllvm_nbva[[i]]$lvs[,1], y= gllvm_nbva[[i]]$lvs[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB F1")+ylab("PPCA-NB F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))

    adonis_result_dis = adonis2(dist(ZIP[[i]]$lvs$factor_scores2)~group[[i]],method = "euclidean")
    pp10 <- ggplot(data.frame(ZIP[[i]]$lvs$factor_score2),aes(x=ZIP[[i]]$lvs$factor_scores2[,1], y=ZIP[[i]]$lvs$factor_scores2[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("ZIPPCA-Poi F1")+ylab("ZIPPCA-Poi F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB[[i]]$lvs$factor_scores2)~group[[i]],method = "euclidean")
    pp11 <- ggplot(data.frame(ZINB[[i]]$lvs$factor_scores2),aes(x=ZINB[[i]]$lvs$factor_scores2[,1], y=ZINB[[i]]$lvs$factor_scores2[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
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
    
    adonis_result_dis = adonis2(dist(ZINB_pca[[i]]$x[,1:2])~group[[i]],method = "euclidean")
    pp12 <- ggplot(data.frame(ZINB_pca[[i]]$x[,1:2]),aes(ZINB_pca[[i]]$x[,1], y=ZINB_pca[[i]]$x[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb_pca F1")+ylab("mbDenoise-zinb_pca F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp13 <- ggplot(data.frame(ZINB_pcoa[[i]]$vectors[,1:2]),aes(ZINB_pcoa[[i]]$vectors[,1], y=ZINB_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb_pcoa F1")+ylab("mbDenoise-zinb_pcoa F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    pp14 <- ggplot(data.frame(ZINB_tsne[[i]]$Y),aes(ZINB_tsne[[i]]$Y[,1], y=ZINB_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb_tsne F1")+ylab("mbDenoise-zinb_tsne F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    n.nn <- nrow(X[[i]])
    n.mm <- ncol(X[[i]])
    ppcanb <- tryCatch({exp(gllvm_nbva[[i]]$lvs %*% t(gllvm_nbva[[i]]$params$theta) +matrix(gllvm_nbva[[i]]$params$beta0,n.nn,n.mm,byrow=TRUE))},
                       error=function(e){NaN})
    X_new <- log2(ppcanb+1)
    ppca_pca[[i]] <- prcomp(X_new)
    dist_bray <- vegan::vegdist(X_new, method="bray")
    ppca_pcoa[[i]] <- ape::pcoa(dist_bray)
    set.seed(4)
    ppca_tsne[[i]] <- Rtsne::Rtsne(X_new,initial_dims = ncol(X_new),perplexity=round((nrow(X_new)-2)/3))
    
    adonis_result_dis = adonis2(dist(ppca_pca[[i]]$x[,1:2])~group[[i]],method = "euclidean")
    pp15 <- ggplot(data.frame(ppca_pca[[i]]$x[,1:2]),aes(ppca_pca[[i]]$x[,1], y=ppca_pca[[i]]$x[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB_pca F1")+ylab("PPCA-NB_pca F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ppca_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp16 <- ggplot(data.frame(ppca_pcoa[[i]]$vectors[,1:2]),aes(ppca_pcoa[[i]]$vectors[,1], y=ppca_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB_pcoa F1")+ylab("PPCA-NB_pcoa F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ppca_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    pp17 <- ggplot(data.frame(ppca_tsne[[i]]$Y),aes(ppca_tsne[[i]]$Y[,1], y=ppca_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB_tsne F1")+ylab("PPCA-NB_tsne F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    X_new <-  log2(t(sav[[i]]$estimate)+1)
    sav_pca[[i]] <- prcomp(X_new)
    dist_bray <- vegan::vegdist(X_new, method="bray")
    sav_pcoa[[i]] <- ape::pcoa(dist_bray)
    set.seed(4)
    sav_tsne[[i]] <- Rtsne::Rtsne(X_new,initial_dims = ncol(X_new),perplexity=round((nrow(X_new)-2)/3))
    
    adonis_result_dis = adonis2(dist(sav_pca[[i]]$x[,1:2])~group[[i]],method = "euclidean")
    pp18 <- ggplot(data.frame(sav_pca[[i]]$x[,1:2]),aes(sav_pca[[i]]$x[,1], y=sav_pca[[i]]$x[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("SAVER_pca F1")+ylab("SAVER_pca F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(sav_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp19 <- ggplot(data.frame(sav_pcoa[[i]]$vectors[,1:2]),aes(sav_pcoa[[i]]$vectors[,1], y=sav_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("SAVER_pcoa F1")+ylab("SAVER_pcoa F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(sav_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    pp20 <- ggplot(data.frame(sav_tsne[[i]]$Y),aes(sav_tsne[[i]]$Y[,1], y=sav_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("SAVER_tsne F1")+ylab("SAVER_tsne F2")+ labs(colour="Location") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("Bhopal","Kerala"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
  }
  
  ggarrange(pp2,pp6,pp15,pp18,pp12,pp3,pp7,pp16,pp19,pp13,pp5,pp11,pp17,pp20,pp14,nrow=3,ncol = 5, common.legend = TRUE,
            legend = "right")
  
}
## Figure 6b: composition and alpha diversity
{
  { 
    diversity <- function(X){
      sh_ele <- -X * log(X)
      sh <- rowSums(sh_ele)
      sp  <- 1-rowSums(X * X)
      
      return(list(Shannon  = sh, Simpson  = sp))
    }
    zip_div <- zinb_div <- zr_div <- svt_div <- pg_div <-  sav_div <-  dmn_div <-ppcanbva_div <- list()
    data1 <- data2 <- data3 <-data4<-data5<-data6 <-data7<-data8<-data9 <- data10<-data11 <-data_all <- list()
    for(i in 1:n){
      
      n.nn <- nrow(X[[i]])
      n.mm <- ncol(X[[i]])

      ppcanb_q <- tryCatch({exp(gllvm_nbva[[i]]$lvs %*% t(gllvm_nbva[[i]]$params$theta) + matrix(gllvm_nbva[[i]]$params$beta0,n.nn,n.mm,byrow=TRUE))/rowSums(exp(gllvm_nbva[[i]]$lvs %*% t(gllvm_nbva[[i]]$params$theta) + matrix(gllvm_nbva[[i]]$params$beta0,n.nn,n.mm,byrow=TRUE)))},
                           error=function(e){NaN})
      ppcanbva_div[[i]] <- diversity(ppcanb_q)
      
      alpha3 <- dmn[[i]]@fit$Estimate
      comesti3_s <- list()
      for (j in 1:ncol(alpha3)) {
        comesti3_s[[j]] <- t(apply(X[[i]],1,function(x){ (x+alpha3[,j])/(sum(x+alpha3[,j]))}))*dmn[[i]]@group[,j]
      }
      QQQ_dmn <- Reduce('+', comesti3_s)
      dmn_div[[i]] <- diversity(QQQ_dmn)

      zip_div[[i]] <- tryCatch({diversity(ZIP_cov[[i]]$Q)},
                               error=function(e){ NaN})
      zinb_div[[i]] <- diversity(ZINB[[i]]$Q)
      zr_div[[i]] <- diversity(ZR[[i]])
      svt_div[[i]] <- diversity(SVT[[i]])
      pg_div[[i]] <- diversity(pg[[i]])
      sav_div[[i]] <- diversity(zr(W =t(sav[[i]]$estimate), alpha = 0.5) )
      
      data1[[i]] <- data.frame(melt(zr_div[[i]]),group[[i]],rep("zr",nrow(X[[i]])));
      data2[[i]] <- data.frame(melt(svt_div[[i]]),group[[i]],rep("svt",nrow(X[[i]])));
      data3[[i]] <- data.frame(melt(pg_div[[i]]),group[[i]],rep("pmr",nrow(X[[i]])));
      data6[[i]] <- data.frame(melt(zinb_div[[i]]),group[[i]],rep("mbDenoise-zinb",nrow(X[[i]])));
      data8[[i]] <- data.frame(melt(ppcanbla_div[[i]]),group[[i]],rep("PPCA-NB",nrow(X[[i]])));
      data9[[i]] <- data.frame(melt(sav_div[[i]]),group[[i]],rep("SAVER_zr",nrow(X[[i]])));
      data10[[i]] <- data.frame(melt(dmn_div[[i]]),group[[i]],rep("dmm",nrow(X[[i]])));
 
      colnames(data1[[i]])[1] <-colnames(data2[[i]])[1] <-colnames(data3[[i]])[1]<-colnames(data6[[i]])[1] <- colnames(data8[[i]])[1]<- colnames(data9[[i]])[1] <- colnames(data10[[i]])[1] <-"value";
      colnames(data1[[i]])[2] <-colnames(data2[[i]])[2] <-colnames(data3[[i]])[2] <- colnames(data6[[i]])[2] <-  colnames(data8[[i]])[2]<-colnames(data9[[i]])[2] <-  colnames(data10[[i]])[2] <-"index";
      colnames(data1[[i]])[3] <-colnames(data2[[i]])[3] <-colnames(data3[[i]])[3] <- colnames(data6[[i]])[3] <-  colnames(data8[[i]])[3]<-colnames(data9[[i]])[3] <- colnames(data10[[i]])[3] <-"State";
      colnames(data1[[i]])[4] <-colnames(data2[[i]])[4] <-colnames(data3[[i]])[4] <- colnames(data6[[i]])[4] <-  colnames(data8[[i]])[4]<-colnames(data9[[i]])[4] <- colnames(data10[[i]])[4] <-"Method";
      data_all[[i]] <- rbind(data1[[i]],data2[[i]],data3[[i]],data10[[i]],data8[[i]],data9[[i]],data6[[i]])
      
    }
  }
  shp_ggplots <- function(data){
    ggplot(data,aes(x = factor(State,level=c("Bhopal","Kerala")), y = value,fill=factor(State))) +
      stat_boxplot(geom = "errorbar",width=0.15,color="darkgrey")+ 
      geom_boxplot(outlier.alpha = 0.1,lwd=0.5,color="darkgrey")+ 
      geom_jitter(aes(color=factor(State)),width =0.2,shape = 21,size=2,fill="white")+ 
      scale_fill_manual(values =  brewer.pal(12, "Paired")[9:10])+
      scale_colour_manual(values =  brewer.pal(12, "Paired")[9:10])+
      scale_y_continuous(name = "alpha diversity")+ scale_x_discrete(name = "") +
      guides(fill=guide_legend(title="Location"),color=guide_legend(title="Location"))+
      theme_bw()+
      theme(axis.text.x =element_blank(),axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=13),
            legend.text= element_text(size=11)) +
      stat_compare_means(label = "p.format")+
      facet_grid(index ~ Method,scales="free")+
      theme(strip.text.x = element_text(size = 11),strip.text.y = element_text(size = 13))
    
    
  }
  shp_ggplots(data_all[[1]])
  
  
}
## Figure 6c: DA-upsetR
{
  for(i in 1:n){ 
    if(is.null(ZINB_cov[[i]])){
      zinb_p[[i]] <- NaN;
    }else{
      X_new <- log2(ZINB_cov[[i]]$muz+1)
      n.nn <- nrow(X[[i]])
      n.mm <- ncol(X[[i]])
      ppcanb_cov <- tryCatch({exp(gllvm_nbva_cov[[i]]$lvs %*% t(gllvm_nbva_cov[[i]]$params$theta) + matrix(gllvm_nbva_cov[[i]]$params$Xcoef,n.nn,n.mm,byrow=TRUE)*Z[[i]]+matrix(gllvm_nbva_cov[[i]]$params$beta0,n.nn,n.mm,byrow=TRUE))},
                             error=function(e){NaN})
      X_new2 <- log2(ppcanb_cov+1)
      
      for(l in 1:ncol(X[[i]])){
        zinb_p[[i]][l] <- tryCatch({t.test(X_new[,l]~Z[[i]])$p.value},error=function(e){ NaN})
        nbva_p[[i]][l] <- tryCatch({t.test(X_new2[,l]~Z[[i]])$p.value},error=function(e){ NaN})
        t_test_p[[i]][l] <- tryCatch({t.test(log2(X[[i]][,l]+1)~Z[[i]])$p.value},error=function(e){ NaN})
        saver_p[[i]][l] <- tryCatch({t.test((t(log2(sav[[i]]$estimate+1)))[,l]~Z[[i]])$p.value},error=function(e){ NaN})
        mbimpute_p[[i]][l] <- tryCatch({t.test((log2(mbi[[i]]$imp_count_mat_norm+1))[,l]~Z[[i]])$p.value},error=function(e){ NaN})
        
      }}
    
    nbva_p_adj[[i]] <- p.adjust(nbva_p[[i]],method = "BH")
    t_test_p_adj[[i]] <-  p.adjust(t_test_p[[i]],method = "BH")
    zinb_p_adj1[[i]] <- p.adjust(zinb_p[[i]],method = "BH")
    saver_p_adj1[[i]] <- p.adjust(saver_p[[i]],method = "BH")
    mbimpute_p_adj1[[i]] <- p.adjust(mbimpute_p[[i]],method = "BH")
    
  }
  deseq <- edger<- meta<- zinb1 <- sav1 <- mbimpute1 <- tt <-set_diff <- ppca <- list()
  for(i in 1:n){
    deseq[[i]] <- na.omit(colnames(X[[i]])[(deseq_p_adj[[i]]) <0.05])
    edger[[i]] <- na.omit(colnames(X[[i]])[(edger_p_adj[[i]]) <0.05])
    meta[[i]] <- na.omit(colnames(X[[i]])[(meta_p_adj[[i]]) <0.05])
    zinb1[[i]] <- na.omit(colnames(X[[i]])[(zinb_p_adj1[[i]]) <0.05])
    sav1[[i]] <- na.omit(colnames(X[[i]])[(saver_p_adj1[[i]]) <0.05])
    mbimpute1[[i]] <- na.omit(colnames(X[[i]])[(mbimpute_p_adj1[[i]]) <0.05])
    tt[[i]]<- na.omit(colnames(X[[i]])[(t_test_p_adj[[i]]) <0.05])
    ppca[[i]] <- na.omit(colnames(X[[i]])[(nbva_p_adj[[i]]) <0.05])
    
    names_diff <- c(rep("DESeq2",length(deseq[[i]])),rep("edgeR",length(edger[[i]])),
                    rep("metagenomeSeq",length(meta[[i]])),rep("mbDenoise-zinb-cov",length(zinb1[[i]])),
                    rep("SAVER",length(sav1[[i]])),rep("mbImpute",length(mbimpute1[[i]])),
                    rep("t test",length(tt[[i]])),rep("PPCA-NB",length(ppca[[i]])))
    
    set_diff[[i]] <- data.frame(c(deseq[[i]],edger[[i]],meta[[i]],zinb1[[i]],sav1[[i]],mbimpute1[[i]],tt[[i]],ppca[[i]]),names_diff)
    colnames(set_diff[[i]]) <- c("OTU","Compartments")                         
    
  }

  input <- list()
  for(i in 1:n){
    input[[i]] <- fromList(list(DESeq2=na.omit(deseq[[i]]),edgeR=na.omit(edger[[i]]),metagenomeSeq=na.omit(meta[[i]]),
                                t.test=na.omit(tt[[i]]),mbDenoise_zinb_cov=na.omit(zinb1[[i]]),
                                SAVER=na.omit(sav1[[i]]),mbImpute=na.omit(mbimpute1[[i]])))
  }
  colnames(input[[1]]) <- colnames(input[[2]]) <- c("DESeq2","edgeR","metagenomeSeq","t-test","mbDenoise-zinb-cov","SAVER","mbImpute")
  col <- c(brewer.pal(7, "Set2")[1:7])
  
  upset(input[[1]], nsets = 7, number.angles = 0, point.size = 2.5, line.size = 1,
        mainbar.y.label = "Intersection size", sets.x.label = "No. of DA species", 
        sets.bar.color = brewer.pal(7, "Set2")[1:7],main.bar.color = rev(colorRampPalette(brewer.pal(7, "RdPu"))(18)),
        order.by = c("degree", "freq"),decreasing = c(F,T),
        text.scale = c(1.8, 1.5, 1.3, 1, 1.5, 2),
        queries = list(list(query = intersects, color="#4DAF4A",params = list("mbDenoise-zinb-cov"),active=T),
                       list(query = intersects, color="#377EB8",active=T,params = list("DESeq2","metagenomeSeq","mbDenoise-zinb-cov","mbImpute"))))
  
  
}

