library(reshape)
library(ggplot2)
library(ggpubr)
library(grid)
library(cowplot)

###denoise
#n<p main text eta=0.25
{
  {
    load("~/mbDenoise/reproducibility/simulation/sim_nb1.RData")
    MSE <- log10(cbind(noimpute_mse2,zinb_mse2,zip_mse2,unlist(ppcanb_mse2),sav_mse2,mbi_mse2))
    distance <- log10(cbind(noimpute_wasserstein1d,zinb_wasserstein1d,zip_wasserstein1d,
                            unlist(ppcanb_wasserstein1d),sav_wasserstein1d,mbi_wasserstein1d))
    cor <- cbind(rowMeans(do.call(rbind,noimpute_cor),na.rm=T),rowMeans(do.call(rbind,zinb_cor),na.rm=T),
                 rowMeans(do.call(rbind,zip_cor),na.rm=T),
                 rowMeans(do.call(rbind,ppcanb_cor),na.rm=T),
                 rowMeans(do.call(rbind,sav_cor),na.rm=T),
                 rowMeans(do.call(rbind,mbi_cor),na.rm=T))
    time <- log10(cbind(NaN,difftime_zinb,difftime_zip,difftime_nbva,difftime_sav,difftime_mbi))
    
    colnames(MSE) <- colnames(time) <-colnames(distance) <- colnames(cor) <-
      c("No imputation","mbDenoise-zinb","mbDenoise-zip","PPCA-NB","SAVER","mbImpute")
    
    mse <- rep("MSE",sum(nrow(melt(MSE))))
    mse_melt <- cbind(melt(MSE),mse)
    mse <- rep("time",sum(nrow(melt(time))))
    time_melt <- cbind(melt(time),mse)
    mse <- rep("correlation",sum(nrow(melt(cor))))
    cor_melt <- cbind(melt(cor),mse)
    mse <- rep("Wasserstein",sum(nrow(melt(distance))))
    dis_melt <- cbind(melt(distance),mse)
    
    all_melt <- rbind(mse_melt,dis_melt,cor_melt,time_melt)
    setting <- rep("1",sum(nrow(all_melt)))
    all_melt <- cbind(all_melt,setting)
    colnames(all_melt)[2] <- "method"
  }
  {
    load("~/mbDenoise/reproducibility/simulation/sim_poi1.RData")
    MSE <- log10(cbind(noimpute_mse2,zinb_mse2,zip_mse2,unlist(ppcanb_mse2),sav_mse2,mbi_mse2))
    distance <- log10(cbind(noimpute_wasserstein1d,zinb_wasserstein1d,zip_wasserstein1d,
                            unlist(ppcanb_wasserstein1d),sav_wasserstein1d,mbi_wasserstein1d))
    cor <- cbind(rowMeans(do.call(rbind,noimpute_cor),na.rm=T),rowMeans(do.call(rbind,zinb_cor),na.rm=T),
                 rowMeans(do.call(rbind,zip_cor),na.rm=T),
                 rowMeans(do.call(rbind,ppcanb_cor),na.rm=T),
                 rowMeans(do.call(rbind,sav_cor),na.rm=T),
                 rowMeans(do.call(rbind,mbi_cor),na.rm=T))
    time <- log10(cbind(NaN,difftime_zinb,difftime_zip,difftime_nbva,difftime_sav,difftime_mbi))
    
    colnames(MSE) <- colnames(time) <-colnames(distance) <- colnames(cor) <-
      c("No imputation","mbDenoise-zinb","mbDenoise-zip","PPCA-NB","SAVER","mbImpute")
    mse <- rep("MSE",sum(nrow(melt(MSE))))
    mse_melt <- cbind(melt(MSE),mse)
    mse <- rep("time",sum(nrow(melt(time))))
    time_melt <- cbind(melt(time),mse)
    mse <- rep("correlation",sum(nrow(melt(cor))))
    cor_melt <- cbind(melt(cor),mse)
    mse <- rep("Wasserstein",sum(nrow(melt(distance))))
    dis_melt <- cbind(melt(distance),mse)
    
    all_melt2 <- rbind(mse_melt,dis_melt,cor_melt,time_melt)
    setting <- rep("2",sum(nrow(all_melt2)))
    all_melt2 <- cbind(all_melt2,setting)
    colnames(all_melt2)[2] <- "method"
  }
  {
    load("~/mbDenoise/reproducibility/simulation/sim_lnm1.RData")
    MSE <- log10(cbind(noimpute_mse2,zinb_mse2,zip_mse2,unlist(ppcanb_mse2),sav_mse2,mbi_mse2))
    distance <- log10(cbind(noimpute_wasserstein1d,zinb_wasserstein1d,zip_wasserstein1d,
                            unlist(ppcanb_wasserstein1d),sav_wasserstein1d,mbi_wasserstein1d))
    cor <- cbind(rowMeans(do.call(rbind,noimpute_cor),na.rm=T),rowMeans(do.call(rbind,zinb_cor),na.rm=T),
                 rowMeans(do.call(rbind,zip_cor),na.rm=T),
                 rowMeans(do.call(rbind,ppcanb_cor),na.rm=T),
                 rowMeans(do.call(rbind,sav_cor),na.rm=T),
                 rowMeans(do.call(rbind,mbi_cor),na.rm=T))
    time <- log10(cbind(NaN,difftime_zinb,difftime_zip,difftime_nbva,difftime_sav,difftime_mbi))
    
    colnames(MSE) <- colnames(time) <-colnames(distance) <- colnames(cor) <-
      c("No imputation","mbDenoise-zinb","mbDenoise-zip","PPCA-NB","SAVER","mbImpute")
    mse <- rep("MSE",sum(nrow(melt(MSE))))
    mse_melt <- cbind(melt(MSE),mse)
    mse <- rep("time",sum(nrow(melt(time))))
    time_melt <- cbind(melt(time),mse)
    mse <- rep("correlation",sum(nrow(melt(cor))))
    cor_melt <- cbind(melt(cor),mse)
    mse <- rep("Wasserstein",sum(nrow(melt(distance))))
    dis_melt <- cbind(melt(distance),mse)
    
    all_melt3 <- rbind(mse_melt,dis_melt,cor_melt,time_melt)
    setting <- rep("3",sum(nrow(all_melt3)))
    all_melt3 <- cbind(all_melt3,setting)
    colnames(all_melt3)[2] <- "method"
  }
  {
    load("~/mbDenoise/reproducibility/simulation/sim_zifa1.RData")
    MSE <- log10(cbind(noimpute_mse2,zinb_mse2,zip_mse2,unlist(ppcanb_mse2),sav_mse2,mbi_mse2))
    distance <- log10(cbind(noimpute_wasserstein1d,zinb_wasserstein1d,zip_wasserstein1d,
                            unlist(ppcanb_wasserstein1d),sav_wasserstein1d,mbi_wasserstein1d))
    cor <- cbind(rowMeans(do.call(rbind,noimpute_cor),na.rm=T),rowMeans(do.call(rbind,zinb_cor),na.rm=T),
                 rowMeans(do.call(rbind,zip_cor),na.rm=T),
                 rowMeans(do.call(rbind,ppcanb_cor),na.rm=T),
                 rowMeans(do.call(rbind,sav_cor),na.rm=T),
                 rowMeans(do.call(rbind,mbi_cor),na.rm=T))
    time <- log10(cbind(NaN,difftime_zinb,difftime_zip,difftime_nbva,difftime_sav,difftime_mbi))
    
    colnames(MSE) <- colnames(time) <-colnames(distance) <- colnames(cor) <-
      c("No imputation","mbDenoise-zinb","mbDenoise-zip","PPCA-NB","SAVER","mbImpute")
    
    mse <- rep("MSE",sum(nrow(melt(MSE))))
    mse_melt <- cbind(melt(MSE),mse)
    mse <- rep("time",sum(nrow(melt(time))))
    time_melt <- cbind(melt(time),mse)
    mse <- rep("correlation",sum(nrow(melt(cor))))
    cor_melt <- cbind(melt(cor),mse)
    mse <- rep("Wasserstein",sum(nrow(melt(distance))))
    dis_melt <- cbind(melt(distance),mse)
    
    all_melt7 <- rbind(mse_melt,dis_melt,cor_melt,time_melt)
    setting <- rep("7",sum(nrow(all_melt7)))
    all_melt7 <- cbind(all_melt7,setting)
    colnames(all_melt7)[2] <- "method"
  }
  {
    load("~/mbDenoise/reproducibility/simulation/sim_gllvm1.RData")
    MSE <- log10(cbind(noimpute_mse2,zinb_mse2,zip_mse2,unlist(ppcanb_mse2),sav_mse2,mbi_mse2))
    distance <- log10(cbind(noimpute_wasserstein1d,zinb_wasserstein1d,zip_wasserstein1d,
                            unlist(ppcanb_wasserstein1d),sav_wasserstein1d,mbi_wasserstein1d))
    cor <- cbind(rowMeans(do.call(rbind,noimpute_cor),na.rm=T),rowMeans(do.call(rbind,zinb_cor),na.rm=T),
                 rowMeans(do.call(rbind,zip_cor),na.rm=T),
                 rowMeans(do.call(rbind,ppcanb_cor),na.rm=T),
                 rowMeans(do.call(rbind,sav_cor),na.rm=T),
                 rowMeans(do.call(rbind,mbi_cor),na.rm=T))
    time <- log10(cbind(NaN,difftime_zinb,difftime_zip,difftime_nbva,difftime_sav,difftime_mbi))
    
    colnames(MSE) <- colnames(time) <-colnames(distance) <- colnames(cor) <-
      c("No imputation","mbDenoise-zinb","mbDenoise-zip","PPCA-NB","SAVER","mbImpute")
    mse <- rep("MSE",sum(nrow(melt(MSE))))
    mse_melt <- cbind(melt(MSE),mse)
    mse <- rep("time",sum(nrow(melt(time))))
    time_melt <- cbind(melt(time),mse)
    mse <- rep("correlation",sum(nrow(melt(cor))))
    cor_melt <- cbind(melt(cor),mse)
    mse <- rep("Wasserstein",sum(nrow(melt(distance))))
    dis_melt <- cbind(melt(distance),mse)
    
    all_melt8 <- rbind(mse_melt,dis_melt,cor_melt,time_melt)
    setting <- rep("8",sum(nrow(all_melt8)))
    all_melt8 <- cbind(all_melt8,setting)
    colnames(all_melt8)[2] <- "method"
    
  }
  
}

df <- rbind(all_melt,all_melt2,all_melt3,all_melt7,all_melt8)
df$mse <- factor(df$mse,labels = c("MSE","Wasserstein","Pearson","time"))
df$setting <- factor(df$setting,labels = c("M1","M2","M3","M4","M5"))
colourCount = length(unique(df$method))

ggplot(data=df, aes(x=factor(method,
                             level=c("No imputation","PPCA-NB","mbImpute",
                                     "SAVER","mbDenoise-zip","mbDenoise-zinb")),y=value,
                    fill=factor(method,level=c("No imputation","PPCA-NB","mbImpute","SAVER",
                      "mbDenoise-zip","mbDenoise-zinb")),
                    colour=factor(method,level=c("No imputation","PPCA-NB","mbImpute","SAVER",
                                                 "mbDenoise-zip","mbDenoise-zinb"))))+
  geom_boxplot(outlier.alpha = 0.18,outlier.size = 0.8,lwd=0.2)+
  stat_boxplot(geom = "errorbar",width=0.35)+ 
  scale_fill_manual(values =  c((colorRampPalette(brewer.pal(9, "Set1")[-6])(colourCount))[-1],(colorRampPalette(brewer.pal(9, "Set1")[-6])(colourCount))[1]))+
  scale_colour_manual(values =  c((colorRampPalette(brewer.pal(9, "Set1")[-6])(colourCount))[-1],(colorRampPalette(brewer.pal(9, "Set1")[-6])(colourCount))[1]))+
  guides(fill = guide_legend(title = "Method"),color=guide_legend(title = "Method"))+
  labs(x = "",y = "", title = "")+
  theme_bw()+theme(strip.text =element_text(size = 10),axis.text=element_text(size = 6),legend.title=element_text(size = 8),
                   legend.text = element_text(size = 7),axis.ticks.x =element_blank()) +
  theme(axis.text.x =element_blank()) +
  facet_grid(mse ~ setting,scales="free")


