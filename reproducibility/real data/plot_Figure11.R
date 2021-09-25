setwd("~/mbDenoise/reproducibility")

library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(reshape)
library(vegan)
library(VennDiagram)
library(UpSetR)
library(gridExtra)
library(viridis)
library(hrbrthemes)
library(pheatmap)

load("real data/CRC.RData")

## Figure 11a: DA-upsetR
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

  #venn-DESeq2
  venn.diagram(x =list(ZellerG=na.omit(deseq[[1]]),FengQ=na.omit(deseq[[2]]),YuJ=na.omit(deseq[[3]]),VogtmannE=na.omit(deseq[[4]])), main="DESeq2",filename = "real data/DeSeq2.png", height = 450, width= 450, resolution =300, imagetype="png", 
               lwd=0.6, fill =c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"),cex=0.45, col="white",cat.cex=0.45,main.cex = 0.5,main.pos= c(0.5, 0.08))
  #venn-edgeR
  venn.diagram(x =list(ZellerG=na.omit(edger[[1]]),FengQ=na.omit(edger[[2]]),YuJ=na.omit(edger[[3]]),VogtmannE=na.omit(edger[[4]])), main="edgeR",filename = "real data/edgeR.png", height = 450, width= 450, resolution =300, imagetype="png", 
               lwd=0.6, fill =c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"),cex=0.45, col="white",cat.cex=0.45,main.cex = 0.5,main.pos= c(0.5, 0.08))
  #venn-metagenomeSeq
  venn.diagram(x =list(ZellerG=na.omit(meta[[1]]),FengQ=na.omit(meta[[2]]),YuJ=na.omit(meta[[3]]),VogtmannE=na.omit(meta[[4]])), main="metagenomeSeq",filename = "real data/metagenomeSeq.png", height = 450, width= 450, resolution =300, imagetype="png", 
               lwd=0.6, fill =c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"),cex=0.45, col="white",cat.cex=0.45,main.cex = 0.5,main.pos= c(0.5, 0.08))
  #venn-SAVER
  venn.diagram(x =list(ZellerG=sav1[[1]],FengQ=sav1[[2]],YuJ=na.omit(sav1[[3]]),VogtmannE=sav1[[4]]), main="SAVER",filename = "real data/SAVER.png", height = 450, width= 450, resolution =300, imagetype="png", 
               lwd=0.6, fill =c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"),cex=0.45, col="white",cat.cex=0.45,main.cex = 0.5,main.pos= c(0.5, 0.08))
  #venn-mbImpute
  venn.diagram(x =list(ZellerG=mbimpute1[[1]],FengQ=mbimpute1[[2]],YuJ=na.omit(mbimpute1[[3]]),VogtmannE=mbimpute1[[4]]), main="mbImpute",filename = "real data/mbImpute.png", height = 450, width= 450, resolution =300, imagetype="png", 
               lwd=0.6, fill =c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"),cex=0.45, col="white",cat.cex=0.45,main.cex = 0.5,main.pos= c(0.5, 0.08))
  #venn-mbDenoise-zinb-cov
  venn.diagram(x =list(ZellerG=na.omit(zinb1[[1]]),FengQ=na.omit(zinb1[[2]]),YuJ=na.omit(zinb1[[3]]),VogtmannE=na.omit(zinb1[[4]])), main="mbDenoise-zinb-cov",filename = "real data/mbDenoise-zinb-cov.png", height = 450, width= 450, resolution =300, imagetype="png", 
               lwd=0.6, fill =c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"),cex=0.45, col="white",cat.cex=0.45,main.cex = 0.5,main.pos= c(0.5, 0.08))
  ##t test
  venn.diagram(x =list(ZellerG=na.omit(tt[[1]]),FengQ=na.omit(tt[[2]]),YuJ=na.omit(tt[[3]]),VogtmannE=na.omit(tt[[4]])), main=c(expression(paste(italic(t), "-test"))),filename = "real data/t.test.png", height = 450, width= 450, resolution =300, imagetype="png", 
               lwd=0.6, fill =c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"),cex=0.45, col="white",cat.cex=0.45,main.cex = 0.5,main.pos= c(0.5, 0.08))
  
  
}

## Figure 11b
{
  
  zinb_inter_4 <- 23
  deseq_inter_4 <- 15
  edger_inter_4 <-2
  mbimpute_inter_4 <-9
  sav_inter_4 <- sav_inter_3 <-0
  sav_inter_2 <-3
  sav_inter_1 <- 36
  t_inter_4 <- t_inter_3 <- 0
  meta_inter_4 <- meta_inter_3 <- 0
  meta_inter_2 <- 4
  t_inter_2 <- 3
  zinb_inter_3 <- 32;zinb_inter_2 <- 70
  deseq_inter_3 <- 42;deseq_inter_2 <- 93
  edger_inter_3 <- 19;edger_inter_2 <- 91
  mbimpute_inter_3 <- 15;mbimpute_inter_2 <- 52
  zinb_inter_1 <- 238
  deseq_inter_1 <- 177
  edger_inter_1 <-276
  mbimpute_inter_1 <-161
  meta_inter_1 <- 15
  t_inter_1 <- 36

  # create a dataset
  method <- c(rep("DESeq2" , 4) , rep("edgeR" , 4) , rep("metagenomeSeq" , 4) , rep("t-test" , 4),
              rep("mbDenoise-zinb-cov" , 4) , rep("mbImpute" , 4), rep("SAVER" , 4) )
  Condition <- rep(c("1" , "2" , "3","4") , 7)
  value <- c(deseq_inter_1,deseq_inter_2,deseq_inter_3,deseq_inter_4,
             edger_inter_1,edger_inter_2,edger_inter_3,edger_inter_4,
             meta_inter_1,meta_inter_2,meta_inter_3,meta_inter_4,
             t_inter_1,t_inter_2,t_inter_3,t_inter_4,
             zinb_inter_1,zinb_inter_2,zinb_inter_3,zinb_inter_4,
             mbimpute_inter_1,mbimpute_inter_2,mbimpute_inter_3,mbimpute_inter_4,
             sav_inter_1,sav_inter_2,sav_inter_3,sav_inter_4)
  data <- data.frame(method,Condition,value)
  
  # Small multiple
  library(viridis)
  library(hrbrthemes)
  ggplot(data, aes(fill=Condition, y=value, x=method)) + 
    geom_bar(position="fill",stat="identity") +
    scale_x_discrete(labels=c("DESeq2","edgeR","mbDenoise-zinb-cov","mbImpute","metagenomeSeq",
                              "SAVER",expression(paste(italic(t), "-test"))))+
    scale_fill_viridis(discrete = T,option = "G",direction=-1)+
    ggtitle("") +
    theme_ipsum() +
    xlab("")+
    coord_flip()+
    geom_text(position='fill',aes(label=value), hjust=1,color="red")
  
  
}

## Figure 11c
{
  inter_name <- intersect(intersect(intersect(zinb1[[1]],zinb1[[2]]),zinb1[[3]]),zinb1[[4]])
  inter <- matrix(0,n,length(inter_name),byrow = T)
  for(i in 1:n){
    names(zinb_p_adj1[[i]]) <- colnames(X[[i]])
    inter[i,] <- zinb_p_adj1[[i]][inter_name]
  }
  colnames(inter) <- inter_name
  rownames(inter) <-   c("ZellerG","FengQ","YuJ","VogtmannE")
  inter
  pheatmap((-log10(inter)),cluster_row = F,cluster_col = F,
           color = colorRampPalette((brewer.pal(n = 7, name ="YlGn")))(15),angle_col = "315")
  
}

## Figure 11d
{
  ## 4
  {
  length(intersect(intersect(intersect(intersect(zinb1[[1]],zinb1[[2]]),zinb1[[3]]),zinb1[[4]]),net))
  length(intersect(intersect(intersect(intersect(deseq[[1]],deseq[[2]]),deseq[[3]]),deseq[[4]]),net))
  length(intersect(intersect(intersect(intersect(edger[[1]],edger[[2]]),edger[[3]]),edger[[4]]),net))
  length(intersect(intersect(intersect(intersect(mbimpute1[[1]],mbimpute1[[2]]),mbimpute1[[3]]),mbimpute1[[4]]),net))
  length(intersect(intersect(intersect(intersect(meta[[1]],meta[[2]]),meta[[3]]),meta[[4]]),net))
  length(intersect(intersect(intersect(intersect(tt[[1]],tt[[2]]),tt[[3]]),tt[[4]]),net))
  length(intersect(intersect(intersect(intersect(sav1[[1]],sav1[[2]]),sav1[[3]]),sav1[[4]]),net))
  }
  ## 3
  {
    #zinb: 6
    (intersect(intersect(intersect(zinb1[[1]],zinb1[[2]]),zinb1[[3]]),net))
    (intersect(intersect(intersect(zinb1[[1]],zinb1[[2]]),zinb1[[4]]),net))
    (intersect(intersect(intersect(zinb1[[1]],zinb1[[4]]),zinb1[[3]]),net))
    (intersect(intersect(intersect(zinb1[[4]],zinb1[[2]]),zinb1[[3]]),net))
    
    #deseq: 4
    (intersect(intersect(intersect(deseq[[1]],deseq[[2]]),deseq[[3]]),net))
    (intersect(intersect(intersect(deseq[[1]],deseq[[2]]),deseq[[4]]),net))
    (intersect(intersect(intersect(deseq[[1]],deseq[[4]]),deseq[[3]]),net))
    (intersect(intersect(intersect(deseq[[4]],deseq[[2]]),deseq[[3]]),net))
    #edger:2
    (intersect(intersect(intersect(edger[[1]],edger[[2]]),edger[[3]]),net))
    (intersect(intersect(intersect(edger[[1]],edger[[2]]),edger[[4]]),net))
    (intersect(intersect(intersect(edger[[1]],edger[[4]]),edger[[3]]),net))
    (intersect(intersect(intersect(edger[[4]],edger[[2]]),edger[[3]]),net))
    #mbimpute:3
    (intersect(intersect(intersect(mbimpute1[[1]],mbimpute1[[2]]),mbimpute1[[3]]),net))
    (intersect(intersect(intersect(mbimpute1[[1]],mbimpute1[[2]]),mbimpute1[[4]]),net))
    (intersect(intersect(intersect(mbimpute1[[1]],mbimpute1[[4]]),mbimpute1[[3]]),net))
    (intersect(intersect(intersect(mbimpute1[[4]],mbimpute1[[2]]),mbimpute1[[3]]),net))
    #meta:0
    (intersect(intersect(intersect(meta[[1]],meta[[2]]),meta[[3]]),net))
    (intersect(intersect(intersect(meta[[1]],meta[[2]]),meta[[4]]),net))
    (intersect(intersect(intersect(meta[[1]],meta[[4]]),meta[[3]]),net))
    (intersect(intersect(intersect(meta[[4]],meta[[2]]),meta[[3]]),net))
    #t : 0
    (intersect(intersect(intersect(tt[[1]],tt[[2]]),tt[[3]]),net))
    (intersect(intersect(intersect(tt[[1]],tt[[2]]),tt[[4]]),net))
    (intersect(intersect(intersect(tt[[1]],tt[[4]]),tt[[3]]),net))
    (intersect(intersect(intersect(tt[[4]],tt[[2]]),tt[[3]]),net))
    #sav:0
    (intersect(intersect(intersect(sav1[[1]],sav1[[2]]),sav1[[3]]),net))
    (intersect(intersect(intersect(sav1[[1]],sav1[[2]]),sav1[[4]]),net))
    (intersect(intersect(intersect(sav1[[1]],sav1[[4]]),sav1[[3]]),net))
    (intersect(intersect(intersect(sav1[[4]],sav1[[2]]),sav1[[3]]),net))
  }
  ## 2
  {
    #zinb : 7
    intersect(intersect(zinb1[[1]],zinb1[[2]]),net)
    intersect(intersect(zinb1[[1]],zinb1[[3]]),net)
    intersect(intersect(zinb1[[1]],zinb1[[4]]),net)
    intersect(intersect(zinb1[[3]],zinb1[[2]]),net)
    intersect(intersect(zinb1[[4]],zinb1[[2]]),net)
    intersect(intersect(zinb1[[4]],zinb1[[3]]),net)
    
    ##5
    intersect(intersect(deseq[[1]],deseq[[2]]),net)
    intersect(intersect(deseq[[1]],deseq[[3]]),net)
    intersect(intersect(deseq[[1]],deseq[[4]]),net)
    intersect(intersect(deseq[[3]],deseq[[2]]),net)
    intersect(intersect(deseq[[4]],deseq[[2]]),net)
    intersect(intersect(deseq[[4]],deseq[[3]]),net)
    
    #4
    intersect(intersect(edger[[1]],edger[[2]]),net)
    intersect(intersect(edger[[1]],edger[[3]]),net)
    intersect(intersect(edger[[1]],edger[[4]]),net)
    intersect(intersect(edger[[3]],edger[[2]]),net)
    intersect(intersect(edger[[4]],edger[[2]]),net)
    intersect(intersect(edger[[4]],edger[[3]]),net)
    
    #mbimpute:4
    intersect(intersect(mbimpute1[[1]],mbimpute1[[2]]),net)
    intersect(intersect(mbimpute1[[1]],mbimpute1[[3]]),net)
    intersect(intersect(mbimpute1[[1]],mbimpute1[[4]]),net)
    intersect(intersect(mbimpute1[[3]],mbimpute1[[2]]),net)
    intersect(intersect(mbimpute1[[4]],mbimpute1[[2]]),net)
    intersect(intersect(mbimpute1[[4]],mbimpute1[[3]]),net)
    
    ##2
    intersect(intersect(meta[[1]],meta[[2]]),net)
    intersect(intersect(meta[[1]],meta[[3]]),net)
    intersect(intersect(meta[[1]],meta[[4]]),net)
    intersect(intersect(meta[[3]],meta[[2]]),net)
    intersect(intersect(meta[[4]],meta[[2]]),net)
    intersect(intersect(meta[[4]],meta[[3]]),net)
    
    #t test :2
    intersect(intersect(tt[[1]],tt[[2]]),net)
    intersect(intersect(tt[[1]],tt[[3]]),net)
    intersect(intersect(tt[[1]],tt[[4]]),net)
    intersect(intersect(tt[[3]],tt[[2]]),net)
    intersect(intersect(tt[[4]],tt[[2]]),net)
    intersect(intersect(tt[[4]],tt[[3]]),net)
    
    #saver: 2
    intersect(intersect(sav1[[1]],sav1[[2]]),net)
    intersect(intersect(sav1[[1]],sav1[[3]]),net)
    intersect(intersect(sav1[[1]],sav1[[4]]),net)
    intersect(intersect(sav1[[3]],sav1[[2]]),net)
    intersect(intersect(sav1[[4]],sav1[[2]]),net)
    intersect(intersect(sav1[[4]],sav1[[3]]),net)
    
  }
  
  data_all <- matrix(0,7,3,byrow = T)
  rownames(data_all) <- c("t-test","DESeq2","edgeR","metagenomeSeq","mbDenoise-zinb-cov","mbImpute","SAVER")
  colnames(data_all) <- c("four","three","two")
  data_all[,1] <- c(0,3,0,0,3,3,0)
  data_all[,2] <- c(0,4,2,0,6,3,0)
  data_all[,3] <- c(2,5,4,2,7,4,2)
  
  df <- melt(data_all)
  colnames(df) <- c("method","Condition","value")
  
  levels(df$method) <- 
    c("DESeq2","edgeR",expression(paste("mbDenoise","-zinb-cov")),"mbImpute","metagenomeSeq",
      "SAVER",expression(paste(italic(t), "-test")))
  
  df$Condition<-factor(df$Condition,
                       levels = c("four","three","two"),
                       labels = c("4","3","2"))
  
  ggplot(data=df, aes(x=factor(Condition,level=c("four","three","two")),y=value))+
    geom_bar(stat="identity", position="dodge", aes(fill=Condition))+
    labs(x = "",y = "", title = "")+
    scale_fill_brewer(palette = "Accent",breaks=c("4","3","2"))+
    theme_bw()+
    facet_grid(~factor(method,c(expression(paste(italic(t), "-test")),"DESeq2","edgeR","metagenomeSeq","mbImpute","SAVER",expression(paste("mbDenoise","-zinb-cov")))),
               labeller = label_parsed ,scales="free")+
    coord_cartesian(ylim = c(0, 7.1))+scale_y_continuous(expand=c(0,0),breaks=seq(0,10,1))+
    theme(
      strip.text=element_text(size = 11),
      axis.text.x =element_blank(),
      axis.ticks.x =element_blank(),
      legend.text = element_text(size = 11),
      axis.text.y.left = element_text(size = 12) #修改坐标轴文本大小
    )
}