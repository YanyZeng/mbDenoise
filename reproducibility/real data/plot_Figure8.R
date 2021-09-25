setwd("~/mbDenoise/reproducibility")

library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(reshape)
library(vegan)
library(VennDiagram)
library(UpSetR)
library(gridExtra)

## Figure 8a: ordination and beta diversity
load("real data/CP1.RData")
{
  ZINB_pcoa <- ZINB_pca <- ZINB_tsne <- ZIP_pcoa <- ZIP_pca <- ZIP_tsne <- list()
  ppca_pcoa <- ppca_pca <- ppca_tsne <- sav_pcoa <- sav_pca <-sav_tsne <- list()
  grid_all <- list()
  group[[1]] <- ifelse( group[[1]]=="Healthy","control","CP")
  levels(Z11) <-  c("subg","supra","tongue")
  for(i in 1:1){ 
    set.seed(4)
    tsne_out_y[[i]] <- Rtsne::Rtsne(log2(1+X[[i]]),initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))
    tsne_out_x[[i]] <- Rtsne::Rtsne(X[[i]],initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))
    Z12 <- factor(Z11)
    
    adonis_result_dis = adonis2(dist(pca_out_y[[i]]$x[,1:2])~Z12*group[[i]], by = NULL,method = "euclidean")
    pp2 <- ggplot(data.frame(pca_out_y[[i]]$x[,1:2]),aes(pca_out_y[[i]]$x[,1], y=pca_out_y[[i]]$x[,2],colour=as.factor(Z11))) +
      geom_point( size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("PCA F1")+ylab("PCA F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(pcoa[[i]]$vectors[,1:2])~Z12*group[[i]], by = NULL,method = "euclidean")
    pp3 <- ggplot(data.frame(pcoa[[i]]$vectors[,1:2]),aes(pcoa[[i]]$vectors[,1], y=pcoa[[i]]$vectors[,2],colour=as.factor(Z11))) +
      geom_point( size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("PCoA F1")+ylab("PCoA F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
  
    adonis_result_dis = adonis2(dist(tsne_out_y[[i]]$Y)~Z12*group[[i]], by = NULL,method = "euclidean")
    pp5 <- ggplot(data.frame(tsne_out_y[[i]]$Y),aes(tsne_out_y[[i]]$Y[,1], y=tsne_out_y[[i]]$Y[,2],colour=as.factor(Z11))) +
      geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("t-SNE F1")+ylab("t-SNE F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = tryCatch({adonis2(dist(ZIFA[[i]][[1]])~Z12*group[[i]], by = NULL,method = "euclidean")},error=function(e){ NaN})
    pp6 <- tryCatch({ggplot(data.frame(ZIFA[[i]][[1]]),aes(x=ZIFA[[i]][[1]][,1], y=ZIFA[[i]][[1]][,2],colour=as.factor(Z11))) +
        geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) +  scale_colour_brewer(palette = "Paired")+
        xlab("ZIFA F1")+ylab("ZIFA F2")+ labs(shape="Disease", colour="Site") +
        theme_bw()+
        annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))},error=function(e){ NaN})
    
    
    adonis_result_dis = adonis2(dist(gllvm_nbva[[i]]$lvs)~Z12*group[[i]], by = NULL,method = "euclidean")
    pp7 <- ggplot(data.frame( gllvm_nbva[[i]]$lvs),aes(x= gllvm_nbva[[i]]$lvs[,1], y= gllvm_nbva[[i]]$lvs[,2],colour=as.factor(Z11))) +
      geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("PPCA-NB F1")+ylab("PPCA-NB F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))

    adonis_result_dis = adonis2(dist(ZIP[[i]]$lvs$factor_scores2)~Z12*group[[i]], by = NULL,method = "euclidean")
    pp10 <- ggplot(data.frame(ZIP[[i]]$lvs$factor_scores2),aes(x=ZIP[[i]]$lvs$factor_scores2[,1], y=ZIP[[i]]$lvs$factor_scores2[,2],colour=as.factor(Z11))) +
      geom_point( size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("ZIPPCA-Poi F1")+ylab("ZIPPCA-Poi F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB[[i]]$lvs$factor_scores2)~Z12*group[[i]], by = NULL,method = "euclidean")
    pp11 <- ggplot(data.frame(ZINB[[i]]$lvs$factor_scores2),aes(x=ZINB[[i]]$lvs$factor_scores2[,1], y=ZINB[[i]]$lvs$factor_scores2[,2],colour=as.factor(Z11))) +
      geom_point( size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("mbDenoise-zinb F1")+ylab("mbDenoise-zinb F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    X_new <- log2(ZINB[[i]]$muz+1)
    ZINB_pca[[i]] <- prcomp(X_new)
    dist_bray <- vegan::vegdist(X_new, method="bray")
    ZINB_pcoa[[i]] <- ape::pcoa(dist_bray)
    set.seed(4)
    ZINB_tsne[[i]] <- Rtsne::Rtsne(X_new,initial_dims = ncol(X_new),perplexity=round((nrow(X_new)-2)/3))
    
    adonis_result_dis = adonis2(dist(ZINB_pca[[i]]$x[,1:2])~Z12*group[[i]], by = NULL,method = "euclidean")
    pp12 <- ggplot(data.frame(ZINB_pca[[i]]$x[,1:2]),aes(ZINB_pca[[i]]$x[,1], y=ZINB_pca[[i]]$x[,2],colour=as.factor(Z11))) +
      geom_point( size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("mbDenoise-zinb_pca F1")+ylab("mbDenoise-zinb_pca F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB_pcoa[[i]]$vectors[,1:2])~Z12*group[[i]], by = NULL,method = "euclidean")
    pp13 <- ggplot(data.frame(ZINB_pcoa[[i]]$vectors[,1:2]),aes(ZINB_pcoa[[i]]$vectors[,1], y=ZINB_pcoa[[i]]$vectors[,2],colour=as.factor(Z11))) +
      geom_point( size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("mbDenoise-zinb_pcoa F1")+ylab("mbDenoise-zinb_pcoa F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB_tsne[[i]]$Y)~Z12*group[[i]], by = NULL,method = "euclidean")
    pp14 <- ggplot(data.frame(ZINB_tsne[[i]]$Y),aes(ZINB_tsne[[i]]$Y[,1], y=ZINB_tsne[[i]]$Y[,2],colour=as.factor(Z11))) +
      geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("mbDenoise-zinb_tsne F1")+ylab("mbDenoise-zinb_tsne F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    X_new <-  log2(t(sav[[i]]$estimate)+1)
    sav_pca[[i]] <- prcomp(X_new)
    dist_bray <- vegan::vegdist(X_new, method="bray")
    sav_pcoa[[i]] <- ape::pcoa(dist_bray)
    set.seed(4)
    sav_tsne[[i]] <- Rtsne::Rtsne(X_new,initial_dims = ncol(X_new),perplexity=round((nrow(X_new)-2)/3))
    
    adonis_result_dis = adonis2(dist(sav_pca[[i]]$x[,1:2])~Z12*group[[i]],by = NULL,method = "euclidean")
    pp18 <- ggplot(data.frame(sav_pca[[i]]$x[,1:2]),aes(sav_pca[[i]]$x[,1], y=sav_pca[[i]]$x[,2],colour=as.factor(Z11))) +
      geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("SAVER_pca F1")+ylab("SAVER_pca F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(sav_pcoa[[i]]$vectors[,1:2])~Z12*group[[i]],by = NULL,method = "euclidean")
    pp19 <- ggplot(data.frame(sav_pcoa[[i]]$vectors[,1:2]),aes(sav_pcoa[[i]]$vectors[,1], y=sav_pcoa[[i]]$vectors[,2],colour=as.factor(Z11))) +
      geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("SAVER_pcoa F1")+ylab("SAVER_pcoa F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(sav_tsne[[i]]$Y)~Z12*group[[i]],by = NULL,method = "euclidean")
    pp20 <- ggplot(data.frame(sav_tsne[[i]]$Y),aes(sav_tsne[[i]]$Y[,1], y=sav_tsne[[i]]$Y[,2],colour=as.factor(Z11))) +
      geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("SAVER_tsne F1")+ylab("SAVER_tsne F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
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
    
    adonis_result_dis = adonis2(dist(ppca_pca[[i]]$x[,1:2])~Z12*group[[i]],by = NULL,method = "euclidean")
    pp15 <- ggplot(data.frame(ppca_pca[[i]]$x[,1:2]),aes(ppca_pca[[i]]$x[,1], y=ppca_pca[[i]]$x[,2],colour=as.factor(Z11))) +
      geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("PPCA-NB_pca F1")+ylab("PPCA-NB_pca F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ppca_pcoa[[i]]$vectors[,1:2])~Z12*group[[i]],by = NULL,method = "euclidean")
    pp16 <- ggplot(data.frame(ppca_pcoa[[i]]$vectors[,1:2]),aes(ppca_pcoa[[i]]$vectors[,1], y=ppca_pcoa[[i]]$vectors[,2],colour=as.factor(Z11))) +
      geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("PPCA-NB_pcoa F1")+ylab("PPCA-NB_pcoa F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
    adonis_result_dis = adonis2(dist(ppca_tsne[[i]]$Y)~Z12*group[[i]],by = NULL,method = "euclidean")
    pp17 <- ggplot(data.frame(ppca_tsne[[i]]$Y),aes(ppca_tsne[[i]]$Y[,1], y=ppca_tsne[[i]]$Y[,2],colour=as.factor(Z11))) +
      geom_point(size = 2.5,aes(shape=as.factor(group[[i]]))) + scale_colour_brewer(palette = "Paired")+
      xlab("PPCA-NB_tsne F1")+ylab("PPCA-NB_tsne F2")+ labs(shape="Disease", colour="Site") +
      theme_bw()+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
  }
  
  ggarrange(pp2,pp6,pp15,pp18,pp12,pp3,pp7,pp16,pp19,pp13,pp5,pp11,pp17,pp20,pp14,nrow=3,ncol = 5, common.legend = TRUE,
            legend = "right")
  
}

## Figure 8b-c: ordination and beta diversity
load("real data/CP2.RData")
{
  ZINB_pcoa <- ZINB_pca <- ZINB_tsne <- ZIP_pcoa <- ZIP_pca <- ZIP_tsne <- list()
  ppca_pcoa <- ppca_pca <- ppca_tsne <- sav_pcoa <- sav_pca <-sav_tsne <- list()
  grid_all <- list()
  ## Figure 8b
  group[[1]] <- ifelse( group[[1]]=="Healthy","control","CP")
  for(i in 1:1){ 
    set.seed(4)
    tsne_out_y[[i]] <- Rtsne::Rtsne(log2(1+X[[i]]),initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))
    tsne_out_x[[i]] <- Rtsne::Rtsne(X[[i]],initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))

    adonis_result_dis = adonis2(dist(pca_out_y[[i]]$x[,1:2])~group[[i]],method = "euclidean")
    pp2 <- ggplot(data.frame(pca_out_y[[i]]$x[,1:2]),aes(pca_out_y[[i]]$x[,1], y=pca_out_y[[i]]$x[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PCA2 F1")+ylab("PCA2 F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp3 <- ggplot(data.frame(pcoa[[i]]$vectors[,1:2]),aes(pcoa[[i]]$vectors[,1], y=pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PCoA(Bray-curtis) F1")+ylab("PCoA(Bray-curtis) F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(tsne_out_y[[i]]$Y)~group[[i]],method = "euclidean")
    pp5 <- ggplot(data.frame(tsne_out_y[[i]]$Y),aes(tsne_out_y[[i]]$Y[,1], y=tsne_out_y[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2)+ scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("t-SNE2 F1")+ylab("t-SNE2 F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = tryCatch({adonis2(dist(ZIFA[[i]][[1]])~group[[i]],method = "euclidean")},error=function(e){ NaN})
    pp6 <- tryCatch({ggplot(data.frame(ZIFA[[i]][[1]]),aes(x=ZIFA[[i]][[1]][,1], y=ZIFA[[i]][[1]][,2],colour=as.factor(group[[i]]))) +
        geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
        scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
        xlab("ZIFA F1")+ylab("ZIFA F2")+labs(colour="Disease") +
        theme_bw()+
        stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
        #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
        annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))},error=function(e){ NaN})
    
    adonis_result_dis = adonis2(dist(gllvm_nbva[[i]]$lvs)~group[[i]],method = "euclidean")
    pp7 <- ggplot(data.frame( gllvm_nbva[[i]]$lvs),aes(x= gllvm_nbva[[i]]$lvs[,1], y= gllvm_nbva[[i]]$lvs[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB F1")+ylab("PPCA-NB F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))

    adonis_result_dis = adonis2(dist(ZIP[[i]]$lvs$factor_scores2)~group[[i]],method = "euclidean")
    pp10 <- ggplot(data.frame(ZIP[[i]]$lvs$factor_scores[[2]]),aes(x=ZIP[[i]]$lvs$factor_scores[[2]][,1], y=ZIP[[i]]$lvs$factor_scores[[2]][,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("ZIPPCA-Poi F1")+ylab("ZIPPCA-Poi F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB[[i]]$lvs$factor_scores2)~group[[i]],method = "euclidean")
    pp11 <- ggplot(data.frame(ZINB[[i]]$lvs$factor_scores2),aes(x=ZINB[[i]]$lvs$factor_scores2[,1], y=ZINB[[i]]$lvs$factor_scores2[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb F1")+ylab("mbDenoise-zinb F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
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
      xlab("mbDenoise-zinb_pca F1")+ylab("mbDenoise-zinb_pca F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp13 <- ggplot(data.frame(ZINB_pcoa[[i]]$vectors[,1:2]),aes(ZINB_pcoa[[i]]$vectors[,1], y=ZINB_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb_pcoa F1")+ylab("mbDenoise-zinb_pcoa F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    pp14 <- ggplot(data.frame(ZINB_tsne[[i]]$Y),aes(ZINB_tsne[[i]]$Y[,1], y=ZINB_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb_tsne F1")+ylab("mbDenoise-zinb_tsne F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
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
      xlab("SAVER_pca F1")+ylab("SAVER_pca F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
    adonis_result_dis = adonis2(dist(sav_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp19 <- ggplot(data.frame(sav_pcoa[[i]]$vectors[,1:2]),aes(sav_pcoa[[i]]$vectors[,1], y=sav_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("SAVER_pcoa F1")+ylab("SAVER_pcoa F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
    adonis_result_dis = adonis2(dist(sav_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    pp20 <- ggplot(data.frame(sav_tsne[[i]]$Y),aes(sav_tsne[[i]]$Y[,1], y=sav_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("SAVER_tsne F1")+ylab("SAVER_tsne F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
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
      xlab("PPCA-NB_pca F1")+ylab("PPCA-NB_pca F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ppca_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp16 <- ggplot(data.frame(ppca_pcoa[[i]]$vectors[,1:2]),aes(ppca_pcoa[[i]]$vectors[,1], y=ppca_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB_pcoa F1")+ylab("PPCA-NB_pcoa F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ppca_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    pp17 <- ggplot(data.frame(ppca_tsne[[i]]$Y),aes(ppca_tsne[[i]]$Y[,1], y=ppca_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB_tsne F1")+ylab("PPCA-NB_tsne F2")+ labs(colour="Disease") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("control","CP"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
  }
  ggarrange(pp2,pp6,pp15,pp18,pp12,pp3,pp7,pp16,pp19,pp13,pp5,pp11,pp17,pp20,pp14,nrow=3,ncol = 5, common.legend = TRUE,
            legend = "right")
  
  ## Figure 8c
  group[[2]] <- ifelse( group[[2]]=="Supra","supra","subg")
  for(i in 2:2){ 
    set.seed(4)
    tsne_out_y[[i]] <- Rtsne::Rtsne(log2(1+X[[i]]),initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))
    tsne_out_x[[i]] <- Rtsne::Rtsne(X[[i]],initial_dims = ncol(X[[i]]),perplexity=round((nrow(X[[i]])-2)/3))
    
    adonis_result_dis = adonis2(dist(pca_out_y[[i]]$x[,1:2])~group[[i]],method = "euclidean")
    pp2 <- ggplot(data.frame(pca_out_y[[i]]$x[,1:2]),aes(pca_out_y[[i]]$x[,1], y=pca_out_y[[i]]$x[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PCA2 F1")+ylab("PCA2 F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
    adonis_result_dis = adonis2(dist(pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp3 <- ggplot(data.frame(pcoa[[i]]$vectors[,1:2]),aes(pcoa[[i]]$vectors[,1], y=pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PCoA(Bray-curtis) F1")+ylab("PCoA(Bray-curtis) F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
   
    adonis_result_dis = adonis2(dist(tsne_out_y[[i]]$Y)~group[[i]],method = "euclidean")
    pp5 <- ggplot(data.frame(tsne_out_y[[i]]$Y),aes(tsne_out_y[[i]]$Y[,1], y=tsne_out_y[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("t-SNE2 F1")+ylab("t-SNE2 F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = tryCatch({adonis2(dist(ZIFA[[i]][[1]])~group[[i]],method = "euclidean")},error=function(e){ NaN})
    pp6 <- tryCatch({ggplot(data.frame(ZIFA[[i]][[1]]),aes(x=ZIFA[[i]][[1]][,1], y=ZIFA[[i]][[1]][,2],colour=as.factor(group[[i]]))) +
        geom_point(size = 2) +   scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
        scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
        xlab("ZIFA F1")+ylab("ZIFA F2")+labs(colour="Site") +
        theme_bw()+
        stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
        #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
        annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))},error=function(e){ NaN})
    
    
    adonis_result_dis = adonis2(dist(gllvm_nbva[[i]]$lvs)~group[[i]],method = "euclidean")
    pp7 <- ggplot(data.frame( gllvm_nbva[[i]]$lvs),aes(x= gllvm_nbva[[i]]$lvs[,1], y= gllvm_nbva[[i]]$lvs[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB F1")+ylab("PPCA-NB F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZIP[[i]]$lvs$factor_scores2)~group[[i]],method = "euclidean")
    pp10 <- ggplot(data.frame(ZIP[[i]]$lvs$factor_scores2),aes(x=ZIP[[i]]$lvs$factor_scores2[,1], y=ZIP[[i]]$lvs$factor_scores2[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("ZIPPCA-Poi F1")+ylab("ZIPPCA-Poi F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB[[i]]$lvs$factor_scores2)~group[[i]],method = "euclidean")
    pp11 <- ggplot(data.frame(ZINB[[i]]$lvs$factor_scores2),aes(x=ZINB[[i]]$lvs$factor_scores2[,1], y=ZINB[[i]]$lvs$factor_scores2[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb F1")+ylab("mbDenoise-zinb F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    X_new <- log2(ZINB[[i]]$muz+1)
    ZINB_pca[[i]] <- prcomp(X_new)
    dist_bray <- vegan::vegdist(X_new, method="bray")
    ZINB_pcoa[[i]] <- ape::pcoa(dist_bray)
    set.seed(4)
    ZINB_tsne[[i]] <- Rtsne::Rtsne(X_new,initial_dims = ncol(X_new),perplexity=round((nrow(X_new)-2)/3))
    
    adonis_result_dis = adonis2(dist(ZINB_pca[[i]]$x[,1:2])~group[[i]],method = "euclidean")
    pp12 <- ggplot(data.frame(ZINB_pca[[i]]$x[,1:2]),aes(ZINB_pca[[i]]$x[,1], y=ZINB_pca[[i]]$x[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb_pca F1")+ylab("mbDenoise-zinb_pca F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ZINB_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp13 <- ggplot(data.frame(ZINB_pcoa[[i]]$vectors[,1:2]),aes(ZINB_pcoa[[i]]$vectors[,1], y=ZINB_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb_pcoa F1")+ylab("mbDenoise-zinb_pcoa F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
    adonis_result_dis = adonis2(dist(ZINB_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    pp14 <- ggplot(data.frame(ZINB_tsne[[i]]$Y),aes(ZINB_tsne[[i]]$Y[,1], y=ZINB_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("mbDenoise-zinb_tsne F1")+ylab("mbDenoise-zinb_tsne F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
    X_new <-  log2(t(sav[[i]]$estimate)+1)
    sav_pca[[i]] <- prcomp(X_new)
    dist_bray <- vegan::vegdist(X_new, method="bray")
    sav_pcoa[[i]] <- ape::pcoa(dist_bray)
    set.seed(4)
    sav_tsne[[i]] <- Rtsne::Rtsne(X_new,initial_dims = ncol(X_new),perplexity=round((nrow(X_new)-2)/3))
    
    adonis_result_dis = adonis2(dist(sav_pca[[i]]$x[,1:2])~group[[i]],method = "euclidean")
    pp18 <- ggplot(data.frame(sav_pca[[i]]$x[,1:2]),aes(sav_pca[[i]]$x[,1], y=sav_pca[[i]]$x[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("SAVER_pca F1")+ylab("SAVER_pca F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(sav_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp19 <- ggplot(data.frame(sav_pcoa[[i]]$vectors[,1:2]),aes(sav_pcoa[[i]]$vectors[,1], y=sav_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) +  scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("SAVER_pcoa F1")+ylab("SAVER_pcoa F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
    adonis_result_dis = adonis2(dist(sav_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    pp20 <- ggplot(data.frame(sav_tsne[[i]]$Y),aes(sav_tsne[[i]]$Y[,1], y=sav_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("SAVER_tsne F1")+ylab("SAVER_tsne F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
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
      xlab("PPCA-NB_pca F1")+ylab("PPCA-NB_pca F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    
    adonis_result_dis = adonis2(dist(ppca_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
    pp16 <- ggplot(data.frame(ppca_pcoa[[i]]$vectors[,1:2]),aes(ppca_pcoa[[i]]$vectors[,1], y=ppca_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
      geom_point( size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB_pcoa F1")+ylab("PPCA-NB_pcoa F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
    adonis_result_dis = adonis2(dist(ppca_tsne[[i]]$Y)~group[[i]],method = "euclidean")
    pp17 <- ggplot(data.frame(ppca_tsne[[i]]$Y),aes(ppca_tsne[[i]]$Y[,1], y=ppca_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
      geom_point(size = 2) + scale_colour_manual(values =  brewer.pal(8, "Set2")[3:4])+
      scale_fill_manual(values =  brewer.pal(8,  "Set2")[3:4])+
      xlab("PPCA-NB_tsne F1")+ylab("PPCA-NB_tsne F2")+ labs(colour="Site") +
      theme_bw()+
      stat_ellipse(aes(fill=factor(group[[i]],level=c("subg","supra"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)+
      #stat_ellipse(type = "norm", show.legend = F,linetype=2)+
      annotate("text_npc", npcx = "left", npcy = "top", parse=T,label= paste0('atop(p ==', adonis_result_dis$`Pr(>F)`[1], ', R^2 ==', round(adonis_result_dis$R2[1],2), ')'))
    
  }
  ggarrange(pp2,pp6,pp15,pp18,pp12,pp3,pp7,pp16,pp19,pp13,pp5,pp11,pp17,pp20,pp14,nrow=3,ncol = 5, common.legend = TRUE,
            legend = "right")
  
}


