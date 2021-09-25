setwd("~/mbDenoise/reproducibility")

library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(reshape)
library(vegan)
library(VennDiagram)
library(UpSetR)
library(gridExtra)
library(igraph)
library(Hmisc)
library(psych)
library(pheatmap)

load("real data/CP2.RData")

## Figure 9a-b: DA-upsetR
#DA upsetR
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
  colnames(input[[1]]) <- c("DESeq2","edgeR","metagenomeSeq","t-test","mbDenoise-zinb-cov","SAVER","mbImpute")
  upset(input[[1]], nsets = 7, number.angles = 0, point.size = 2.5, line.size = 1,
        mainbar.y.label = "Intersection size", sets.x.label = "No. of DA species", 
        sets.bar.color = brewer.pal(7, "Set2")[1:7],main.bar.color = rev(colorRampPalette(brewer.pal(7, "RdPu"))(25)),
        order.by = c("degree", "freq"),decreasing = c(F,T),
        text.scale = c(1.8, 1.5, 1.3, 1, 1.5, 2),
        queries = list(list(query = intersects, color="#4DAF4A",params = list("mbDenoise-zinb-cov"),active=T),
                       list(query = intersects, color="#7570B3",params = list("DESeq2"),active=T),
                       list(query = intersects, color="#377EB8",params = list("edgeR"),active=T)))
  
  colnames(input[[2]]) <- c("DESeq2","edgeR","metagenomeSeq","t-test","mbDenoise-zinb-cov","SAVER","mbImpute")
  upset(input[[2]], nsets = 7, number.angles = 0, point.size = 2.5, line.size = 1,
        mainbar.y.label = "Intersection size", sets.x.label = "No. of DA species", 
        sets.bar.color = brewer.pal(7, "Set2")[1:7],main.bar.color = rev(colorRampPalette(brewer.pal(7, "RdPu"))(4)),
        order.by = c("degree", "freq"),decreasing = c(F,T),
        text.scale = c(1.8, 1.5, 1.3, 1, 1.5, 2))
  
}

## Figure 9c: obtain edgelist data as NetShift input
{
  X_case <- X1[Z1==0,]
  X_control <- X1[Z1==1,]
  zerocol <- which(colSums(X_case)==0 | colSums(X_control)==0)
  if(length(zerocol) >0 ){
    X_case <- X_case[,-zerocol];X_control <- X_control[,-zerocol];
  }
  dim(X_case);dim(X_control)
  
  ## get edgelist of X_case 
  {
    occor = corr.test(X_case,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
    occor.r = occor$r # 取相关性矩阵R值
    occor.p = occor$p # 取相关性矩阵p值
    
    # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
    occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 
    dim(occor.r)
    zerocol <- which((colSums(occor.r==0))==(nrow(occor.r)-1))
    if(length(zerocol) >0 ){
      occor.r <- occor.r[,-zerocol];occor.r <- occor.r[-zerocol,];
    }
    dim(occor.r)
    
    # 构建igraph对象
    occor.r[occor.r!=0] = 1
    igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)
    igraph
    bad.vs = V(igraph)[degree(igraph) == 0]
    igraph = delete.vertices(igraph, bad.vs)
    igraph
    # 将igraph weight属性赋值到igraph.weight
    igraph.weight = E(igraph)$weight
    
    # 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
    E(igraph)$weight = NA
    e_case <- get.edgelist(igraph)
  }
  ## get edgelist of X_control 
  {
    occor = corr.test(X_control,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
    occor.r = occor$r # 取相关性矩阵R值
    occor.p = occor$p # 取相关性矩阵p值
    
    # 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
    occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0 
    dim(occor.r)
    zerocol <- which((colSums(occor.r==0))==(nrow(occor.r)-1))
    if(length(zerocol) >0 ){
      occor.r <- occor.r[,-zerocol];occor.r <- occor.r[-zerocol,];
    }
    dim(occor.r)
    
    # 构建igraph对象
    occor.r[occor.r!=0] = 1
    igraph = graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)
    igraph
    bad.vs = V(igraph)[degree(igraph) == 0]
    igraph = delete.vertices(igraph, bad.vs)
    igraph
    # 将igraph weight属性赋值到igraph.weight
    igraph.weight = E(igraph)$weight
    
    # 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
    E(igraph)$weight = NA
    e_control <- get.edgelist(igraph)
  }
  
  ## save 
  write.csv(e_case,file="real data/CP_case.csv" ,  row.names =FALSE, quote =FALSE)
  write.csv(e_control,file="real data/CP_control.csv" ,  row.names =FALSE, quote =FALSE)
  
  
  
}

## Figure 9d
{

  ## three drivers
  net <- c(" s__oral taxon 362"," s__maltophilum"," s__pneumosintes")
  
  intersect(zinb1[[1]],net)
  intersect(deseq[[1]],net)
  intersect(edger[[1]],net)
  intersect(mbimpute1[[1]],net)
  intersect(meta[[1]],net)
  intersect(tt[[1]],net)
  intersect(sav1[[1]],net)

  df <- matrix("absence",3,7)
  #rownames(df) <- net
  colnames(df) <- c("mbDenoise-zinb-cov","DESeq2","edgeR","mbImpute","metagenomeSeq","t-test","SAVER")
  df[,1] <- c("presence","presence","presence")
  df[,2] <- c("presence","presence","presence")
  df[,3] <- c("presence","presence","presence")
  df[,5] <- c("presence","absence","absence")
  df[,7] <- c("presence","absence","absence")
  df[,6] <- c("presence","presence","absence")
  
  df <- data.frame(net,df)
  colnames(df) <- c("net","mbDenoise-zinb-cov","DESeq2","edgeR","mbImpute","metagenomeSeq","t-test","SAVER")
  df <- melt(df,id="net")
  cols=c(
    "presence"="#F781BF","absence"="lightgrey"
  )
  ggplot(df,aes(x=variable,y=net))+
    geom_tile(aes(fill=value),color="white",size=0.5)+ #color和size分别指定方块边线的颜色和粗细
    scale_x_discrete(labels=c("mbDenoise-zinb-cov","DESeq2","edgeR","mbImpute","metagenomeSeq",expression(paste(italic(t), "-test")),
                              "SAVER"),expand = c(0,0))+
    scale_y_discrete("",expand = c(0,0))+
    scale_fill_manual(values=cols)+ #设置填充的颜色
    theme_bw()+
    theme(
      axis.text.x.bottom = element_text(size=13,angle = 45, hjust = 1),
      legend.text = element_text(size = 11),
      axis.text.y.left = element_text(size = 12), #修改坐标轴文本大小
      axis.ticks = element_blank(), #不显示坐标轴刻度
      legend.title = element_blank() #不显示图例title
    )
  
}

## Figure 9e
{
  setdiff(zinb1[[1]], union(union(union(union(union(deseq[[1]],edger[[1]]),meta[[1]]),tt[[1]]),sav1[[1]]),mbimpute1[[1]]))
  setdiff(deseq[[1]], union(union(union(union(union(zinb1[[1]],edger[[1]]),meta[[1]]),tt[[1]]),sav1[[1]]),mbimpute1[[1]]))
  setdiff(edger[[1]], union(union(union(union(union(zinb1[[1]],deseq[[1]]),meta[[1]]),tt[[1]]),sav1[[1]]),mbimpute1[[1]]))
  
  net <- c(" s__oral taxon 362"," s__maltophilum"," s__pneumosintes")
  net_a <- c(" s__endodontalis"," s__FT050"," s__sanguinis"," s__II:C:T1")
  net_all <- c(net,net_a)
  X_net <- X1[,net]
  X_neta <- X1[,net_a]
  
  colnames(X_net) <- c("NetShift: s__oral taxon 362","NetShift: s__maltophilum","NetShift: s__pneumosintes")
  colnames(X_neta) <- c("mbDenoise-zinb-cov: s__endodontalis","mbDenoise-zinb-cov: s__FT050",
                        "DESeq2: s__sanguinis","edgeR: s__II:C:T1")

  cor_zinb <- psych::corr.test(X_net,X_neta,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
  cor_zinb.r = round(cor_zinb$r,2) # 取相关性矩阵R值
  cor_zinb.p = cor_zinb$p
  
  bk <- c(seq(-1,-0.1,by=0.1),seq(0,1,by=0.1))
  p_m <- matrix(0,nrow(cor_zinb.p),ncol(cor_zinb.p))
  p_m <- ifelse(0.01<cor_zinb.p & cor_zinb.p<0.05,"*","")
  p_m <- ifelse(0.01>cor_zinb.p,"**",p_m)
  
  pheatmap(t(cor_zinb.r),color = rev(colorRampPalette((brewer.pal(n = 11, name ="RdYlBu")))(20)),cluster_rows = F,cluster_cols = F,
           legend_breaks=seq(-1,1,0.2),border_color =NA,
           breaks=bk,display_numbers =t(matrix(paste(cor_zinb.r,p_m),nrow(cor_zinb.r),ncol(cor_zinb.r))),fontsize = 14,
           number_color="black",fontsize_number = 1.2* 10,
           angle_col = "315")
  
}


