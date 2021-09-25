library(reshape)
library(ggplot2)
library(ggpubr)
library(grid)
library(cowplot)

#n<p main text eta=0.25 Q
{
  load("~/mbDenoise/reproducibility/simulation/sim_nb1.RData")
  {
    a_bind <- cbind(a_0.5,a_svt,a_multi,a_dmn,a_poi,a_nb,a_ppcanb,a_sav)
    b_bind <- cbind(b_0.5,b_svt,b_multi,b_dmn,b_poi,b_nb,b_ppcanb,b_sav)
    c_bind <- cbind(c_0.5,c_svt,c_multi,c_dmn,c_poi,c_nb,c_ppcanb,c_sav)
    d_bind <- cbind(d_0.5,d_svt,d_multi,d_dmn,d_poi,d_nb,d_ppcanb,d_sav)
    time_bind <- cbind(difftime_zr,difftime_svt,difftime_pg,difftime_dmn,difftime_zip,difftime_zinb,
                       difftime_nbva,difftime_sav)
    
    colnames(a_bind) <- colnames(b_bind) <-colnames(time_bind) <- colnames(c_bind) <- colnames(d_bind) <-
      c("zr","svt","pmr","dmm","ZIPPCA-Poi","ZIPPCA-NB","PPCA-NB","SAVER")
    
    mse <- rep("A ",sum(nrow(melt(a_bind))))
    a_melt <- cbind(melt(a_bind),mse)
    mse <- rep("B",sum(nrow(melt(b_bind))))
    b_melt <- cbind(melt(b_bind),mse)
    mse <- rep("C",sum(nrow(melt(c_bind))))
    c_melt <- cbind(melt(c_bind),mse)
    mse <- rep("D",sum(nrow(melt(d_bind))))
    d_melt <- cbind(melt(d_bind),mse)
    mse <- rep("time",sum(nrow(melt(time_bind))))
    time_melt <- cbind(melt(time_bind),mse)
    
    all_melt <- rbind(a_melt,b_melt,c_melt,d_melt,time_melt)
    setting <- rep("1",sum(nrow(all_melt)))
    all_melt <- cbind(all_melt,setting)
    colnames(all_melt)[2] <- "method"
  }
  load("~/mbDenoise/reproducibility/simulation/sim_poi1.RData")
  {
    a_bind <- cbind(a_0.5,a_svt,a_multi,a_dmn,a_poi,a_nb,a_ppcanb,a_sav)
    b_bind <- cbind(b_0.5,b_svt,b_multi,b_dmn,b_poi,b_nb,b_ppcanb,b_sav)
    c_bind <- cbind(c_0.5,c_svt,c_multi,c_dmn,c_poi,c_nb,c_ppcanb,c_sav)
    d_bind <- cbind(d_0.5,d_svt,d_multi,d_dmn,d_poi,d_nb,d_ppcanb,d_sav)
    time_bind <- cbind(difftime_zr,difftime_svt,difftime_pg,difftime_dmn,difftime_zip,difftime_zinb,
                       difftime_nbva,difftime_sav)
    
    colnames(a_bind) <- colnames(b_bind) <-colnames(time_bind) <- colnames(c_bind) <- colnames(d_bind) <-
      c("zr","svt","pmr","dmm","ZIPPCA-Poi","ZIPPCA-NB","PPCA-NB","SAVER")
    
    mse <- rep("A ",sum(nrow(melt(a_bind))))
    a_melt <- cbind(melt(a_bind),mse)
    mse <- rep("B",sum(nrow(melt(b_bind))))
    b_melt <- cbind(melt(b_bind),mse)
    mse <- rep("C",sum(nrow(melt(c_bind))))
    c_melt <- cbind(melt(c_bind),mse)
    mse <- rep("D",sum(nrow(melt(d_bind))))
    d_melt <- cbind(melt(d_bind),mse)
    mse <- rep("time",sum(nrow(melt(time_bind))))
    time_melt <- cbind(melt(time_bind),mse)
    
    all_melt2 <- rbind(a_melt,b_melt,c_melt,d_melt,time_melt)
    setting <- rep("2",sum(nrow(all_melt2)))
    all_melt2 <- cbind(all_melt2,setting)
    colnames(all_melt2)[2] <- "method"
  }
  load("~/mbDenoise/reproducibility/simulation/sim_lnm1.RData")
  {
    a_bind <- cbind(a_0.5,a_svt,a_multi,a_dmn,a_poi,a_nb,a_ppcanb,a_sav)
    b_bind <- cbind(b_0.5,b_svt,b_multi,b_dmn,b_poi,b_nb,b_ppcanb,b_sav)
    c_bind <- cbind(c_0.5,c_svt,c_multi,c_dmn,c_poi,c_nb,c_ppcanb,c_sav)
    d_bind <- cbind(d_0.5,d_svt,d_multi,d_dmn,d_poi,d_nb,d_ppcanb,d_sav)
    time_bind <- cbind(difftime_zr,difftime_svt,difftime_pg,difftime_dmn,difftime_zip,difftime_zinb,
                       difftime_nbva,difftime_sav)
    
    colnames(a_bind) <- colnames(b_bind) <-colnames(time_bind) <- colnames(c_bind) <- colnames(d_bind) <-
      c("zr","svt","pmr","dmm","ZIPPCA-Poi","ZIPPCA-NB","PPCA-NB","SAVER")
    
    mse <- rep("A ",sum(nrow(melt(a_bind))))
    a_melt <- cbind(melt(a_bind),mse)
    mse <- rep("B",sum(nrow(melt(b_bind))))
    b_melt <- cbind(melt(b_bind),mse)
    mse <- rep("C",sum(nrow(melt(c_bind))))
    c_melt <- cbind(melt(c_bind),mse)
    mse <- rep("D",sum(nrow(melt(d_bind))))
    d_melt <- cbind(melt(d_bind),mse)
    mse <- rep("time",sum(nrow(melt(time_bind))))
    time_melt <- cbind(melt(time_bind),mse)
    
    all_melt3 <- rbind(a_melt,b_melt,c_melt,d_melt,time_melt)
    setting <- rep("3",sum(nrow(all_melt3)))
    all_melt3 <- cbind(all_melt3,setting)
    colnames(all_melt3)[2] <- "method"
  }
  load("~/mbDenoise/reproducibility/simulation/sim_pmr1.RData")
  {
    a_bind <- cbind(a_0.5,a_svt,a_multi,a_dmn,a_poi,a_nb,a_ppcanb,a_sav)
    b_bind <- cbind(b_0.5,b_svt,b_multi,b_dmn,b_poi,b_nb,b_ppcanb,b_sav)
    c_bind <- cbind(c_0.5,c_svt,c_multi,c_dmn,c_poi,c_nb,c_ppcanb,c_sav)
    d_bind <- cbind(d_0.5,d_svt,d_multi,d_dmn,d_poi,d_nb,d_ppcanb,d_sav)
    time_bind <- cbind(difftime_zr,difftime_svt,difftime_pg,difftime_dmn,difftime_zip,difftime_zinb,
                       difftime_nbva,difftime_sav)
    
    colnames(a_bind) <- colnames(b_bind) <-colnames(time_bind) <- colnames(c_bind) <- colnames(d_bind) <-
      c("zr","svt","pmr","dmm","ZIPPCA-Poi","ZIPPCA-NB","PPCA-NB","SAVER")
    
    mse <- rep("A ",sum(nrow(melt(a_bind))))
    a_melt <- cbind(melt(a_bind),mse)
    mse <- rep("B",sum(nrow(melt(b_bind))))
    b_melt <- cbind(melt(b_bind),mse)
    mse <- rep("C",sum(nrow(melt(c_bind))))
    c_melt <- cbind(melt(c_bind),mse)
    mse <- rep("D",sum(nrow(melt(d_bind))))
    d_melt <- cbind(melt(d_bind),mse)
    mse <- rep("time",sum(nrow(melt(time_bind))))
    time_melt <- cbind(melt(time_bind),mse)
    
    all_melt13 <- rbind(a_melt,b_melt,c_melt,d_melt,time_melt)
    setting <- rep("16",sum(nrow(all_melt13)))
    all_melt13 <- cbind(all_melt13,setting)
    colnames(all_melt13)[2] <- "method"
  }
}

df <- rbind(all_melt,all_melt2,all_melt3,all_melt13)
df$value <- log10(df$value)
df$mse <- factor(df$mse,labels =c("C1","C2","C3","C4","time"))
df$setting <- factor(df$setting,labels = c("M1","M2","M3","M6"))
df$method <- factor(df$method,labels = c("dmm","pmr","PPCA-NB","SAVER_zr","svt","mbDenoise-zinb","mbDenoise-zip", "zr"))
colourCount = length(unique(df$method))

ggplot(data=df, aes(x=factor(method,
                             level=c("zr","svt","pmr","dmm","PPCA-NB","SAVER_zr","mbDenoise-zip","mbDenoise-zinb")),
                    y=value,fill=factor(method,
                                        level=c("zr","svt","pmr","dmm","PPCA-NB","SAVER_zr","mbDenoise-zip","mbDenoise-zinb")),
                    colour=factor(method,
                                  level=c("zr","svt","pmr","dmm","PPCA-NB","SAVER_zr","mbDenoise-zip","mbDenoise-zinb"))))+
  geom_boxplot(outlier.alpha = 0.18,outlier.size = 0.8,lwd=0.2)+
  stat_boxplot(geom = "errorbar",width=0.35)+ 
  scale_fill_manual(values =  c(colorRampPalette(brewer.pal(7, "Paired"))(colourCount)[-7],colorRampPalette(brewer.pal(7, "Paired"))(colourCount)[7]))+
  scale_colour_manual(values =  c(colorRampPalette(brewer.pal(7, "Paired"))(colourCount)[-7],colorRampPalette(brewer.pal(7, "Paired"))(colourCount)[7]))+
  guides(fill = guide_legend(title = "Method"),color=guide_legend(title = "Method"))+
  labs(x = "",y = "", title = "")+
  theme_bw()+theme(strip.text =element_text(size = 10),axis.text=element_text(size = 6),legend.title=element_text(size = 8),
                   legend.text = element_text(size = 7),axis.ticks.x =element_blank()) +
  theme(axis.text.x =element_blank()) +
  facet_grid(mse~setting,scales="free")
