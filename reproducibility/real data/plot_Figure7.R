setwd("~/mbDenoise/reproducibility")

library(stringr)
library(curatedMetagenomicData)
library(pheatmap)

load("real data/India_diet.RData")
# ##get DhakanDB_path.csv
# {
# DhakanDB <- curatedMetagenomicData("DhakanDB_2019.pathabundance_relab.stool", counts=TRUE, dryrun=FALSE)
# DhakanDB_path  <- exprs(DhakanDB[[1]])
# DhakanDB_path <- DhakanDB_path[,rownames(X1)]
# 
# zerorow <- which(rowSums(DhakanDB_path)==0 )
# if(length(zerocol) >0 ){
#   DhakanDB_path <- DhakanDB_path[-zerorow,];
# }
# dim(DhakanDB_path);
# DhakanDB_path2 <- zr(DhakanDB_path)
# write.csv(DhakanDB_path2,file="real data/DhakanDB_path.csv" ,  row.names =TRUE, quote =TRUE)
# }

lefse_diff_path <- read_csv("real data/diff_path_lefse.csv",col_names=F)
DhakanDB_path_new <- read_csv("real data/DhakanDB_path.csv")
rownames(DhakanDB_path_new) <- DhakanDB_path_new$X1
diff_path <- lefse_diff_path[lefse_diff_path$X5!="-",]
diff_path_new <- diff_path[order(as.numeric(diff_path$X5),decreasing = F),]
diff_path_top <- diff_path_new[1:20,]

DhakanDB_path_diff <- DhakanDB_path_new[diff_path_top$X1,]
path_diff <- (DhakanDB_path_diff[,-(1:2)])
rownames(path_diff) <- DhakanDB_path_diff$X1
zinb_diff <- t(X1[,zinb1[[1]]])

##Figure 7a
cor_zinb = corr.test(t(zinb_diff),t(path_diff),use="pairwise",method="spearman",adjust="BH",alpha=.05)
cor_zinb.r = round(cor_zinb$r,2) # 取相关性矩阵R值
cor_zinb.p = cor_zinb$p
bk <- c(seq(-1,-0.1,by=0.1),seq(0,1,by=0.1))
p_m <- matrix(0,nrow(cor_zinb.p),ncol(cor_zinb.p))
p_m <- ifelse(0.01<cor_zinb.p & cor_zinb.p<0.05,"*","")
p_m <- ifelse(0.01>cor_zinb.p,"**",p_m)
pheatmap(cor_zinb.r,color = rev(colorRampPalette((brewer.pal(n = 11, name ="PiYG")))(20)),
         legend_breaks=seq(-1,1,0.2),border_color =NA,
         breaks=bk,display_numbers =matrix(paste(cor_zinb.r,p_m),nrow(cor_zinb.r),ncol(cor_zinb.r)),fontsize = 14, 
         number_color="white",fontsize_number = 1.2 * 10,
         angle_col = "315")

##Figure 7b
set <- c("s__Prevotella_copri","s__Lactobacillus_ruminis", "s__Veillonella_unclassified")
se1 <- zinb_diff[set,]
se2 <- zinb_diff[-which(rownames(zinb_diff)==set),];
se2 <- se2[-18,]
cor_zinb = corr.test(t(se1),t(se2),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
cor_zinb.r = round(cor_zinb$r,2) # 取相关性矩阵R值
cor_zinb.p = cor_zinb$p
bk <- c(seq(-1,-0.1,by=0.1),seq(0,1,by=0.1))
p_m <- matrix(0,nrow(cor_zinb.p),ncol(cor_zinb.p))
p_m <- ifelse(0.01<cor_zinb.p & cor_zinb.p<0.05,"*","")
p_m <- ifelse(0.01>cor_zinb.p,"**",p_m)
pheatmap(t(cor_zinb.r),color = rev(colorRampPalette((brewer.pal(n = 11, name ="RdYlBu")))(20)),
         legend_breaks=seq(-1,1,0.2),border_color =NA,
         breaks=bk,display_numbers =t(matrix(paste(cor_zinb.r,p_m),nrow(cor_zinb.r),ncol(cor_zinb.r))),fontsize = 14, 
         number_color="white",fontsize_number = 1.2* 10,
         angle_col = "315")




