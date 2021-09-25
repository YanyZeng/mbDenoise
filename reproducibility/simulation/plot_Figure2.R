library(reshape)
library(ggplot2)
library(ggpubr)
library(grid)
library(cowplot)

##estimation and prediction
##n<p:f+b main text eta=0.25
{
  load("~/mbDenoise/reproducibility/simulation/sim_nb1.RData")
  {
    m_bind <- cbind(m_pr_y,m_pcoa_y,m_tsne_y,m_zifa,m_nbva,m_zip,m_zinb)
    m_bind2 <- cbind(m_pr_y2,m_pcoa_y2,m_tsne_y2,m_zifa2,m_nbva2,m_zip2,m_zinb2)
    
    beta_bind <- cbind(beta_pr_y,NaN,NaN,beta_zifa,beta_nbva,beta_zip,
                       beta_zinb)
    beta_bind2 <- cbind(beta_pr_y2,NaN,NaN,beta_zifa2,beta_nbva2,beta_zip2,
                        beta_zinb2)
    gamma_bind <- cbind(NaN,NaN,NaN,NaN,NaN,gamma_nbva_cov,NaN,gamma_zip_cov,NaN,
                        gamma_zinb_cov,NaN,gamma_zilnm_cov)
    
    time_bind <- log10(cbind(difftime_pr_y,difftime_pcoa_y,difftime_tsne_y,difftime_zifa,difftime_nbva,difftime_zip,difftime_zinb))
    
    colnames(m_bind) <- colnames(m_bind2) <-colnames(time_bind) <- colnames(beta_bind) <- colnames(beta_bind2) <-
      c("PCA","PCoA","tsne","ZIFA","GLLVM-NB","ZIPPCA-Poi","ZIPPCA-NB")
    
    # mse <- rep(paste("MSE1",expression(f_i),sep = ":"),sum(nrow(melt(m_bind))))
    # m_melt <- cbind(melt(m_bind),mse)
    # mse <- rep(paste("MSE2",expression(f_i),sep = ":"),sum(nrow(melt(m_bind2))))
    # m_melt2 <- cbind(melt(m_bind2),mse)
    # mse <- rep(paste("MSE1",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind))))
    # beta_melt <- cbind(melt(beta_bind),mse)
    # mse <- rep(paste("MSE2",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind2))))
    # beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("F_MSE1 ",sum(nrow(melt(m_bind))))
    m_melt <- cbind(melt(m_bind),mse)
    mse <- rep("F_MSE2",sum(nrow(melt(m_bind2))))
    m_melt2 <- cbind(melt(m_bind2),mse)
    mse <- rep("B_MSE1",sum(nrow(melt(beta_bind))))
    beta_melt <- cbind(melt(beta_bind),mse)
    mse <- rep("B_MSE2",sum(nrow(melt(beta_bind2))))
    beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("time",sum(nrow(melt(time_bind))))
    time_melt <- cbind(melt(time_bind),mse)
    
    all_melt <- rbind(m_melt,m_melt2,beta_melt,beta_melt2,time_melt)
    setting <- rep("1",sum(nrow(all_melt)))
    all_melt <- cbind(all_melt,setting)
    colnames(all_melt)[2] <- "method"
  }
  load("~/mbDenoise/reproducibility/simulation/sim_poi1.RData")
  {
    m_bind <- cbind(m_pr_y,m_pcoa_y,m_tsne_y,m_zifa,m_nbva,m_zip,m_zinb)
    m_bind2 <- cbind(m_pr_y2,m_pcoa_y2,m_tsne_y2,m_zifa2,m_nbva2,m_zip2,m_zinb2)
    
    beta_bind <- cbind(beta_pr_y,NaN,NaN,beta_zifa,beta_nbva,beta_zip,
                       beta_zinb)
    beta_bind2 <- cbind(beta_pr_y2,NaN,NaN,beta_zifa2,beta_nbva2,beta_zip2,
                        beta_zinb2)
    gamma_bind <- cbind(NaN,NaN,NaN,NaN,NaN,gamma_nbva_cov,NaN,gamma_zip_cov,NaN,
                        gamma_zinb_cov,NaN,gamma_zilnm_cov)
    
    time_bind <- log10(cbind(difftime_pr_y,difftime_pcoa_y,difftime_tsne_y,difftime_zifa,difftime_nbva,difftime_zip,difftime_zinb))
    
    colnames(m_bind) <- colnames(m_bind2) <-colnames(time_bind) <- colnames(beta_bind) <- colnames(beta_bind2) <-
      c("PCA","PCoA","tsne","ZIFA","GLLVM-NB","ZIPPCA-Poi","ZIPPCA-NB")
    
    # mse <- rep(paste("MSE1",expression(f_i),sep = ":"),sum(nrow(melt(m_bind))))
    # m_melt <- cbind(melt(m_bind),mse)
    # mse <- rep(paste("MSE2",expression(f_i),sep = ":"),sum(nrow(melt(m_bind2))))
    # m_melt2 <- cbind(melt(m_bind2),mse)
    # mse <- rep(paste("MSE1",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind))))
    # beta_melt <- cbind(melt(beta_bind),mse)
    # mse <- rep(paste("MSE2",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind2))))
    # beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("F_MSE1 ",sum(nrow(melt(m_bind))))
    m_melt <- cbind(melt(m_bind),mse)
    mse <- rep("F_MSE2",sum(nrow(melt(m_bind2))))
    m_melt2 <- cbind(melt(m_bind2),mse)
    mse <- rep("B_MSE1",sum(nrow(melt(beta_bind))))
    beta_melt <- cbind(melt(beta_bind),mse)
    mse <- rep("B_MSE2",sum(nrow(melt(beta_bind2))))
    beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("G_MSE1 ",sum(nrow(melt(gamma_bind))))
    g_melt <- cbind(melt(gamma_bind),mse)
    mse <- rep("time",sum(nrow(melt(time_bind))))
    time_melt <- cbind(melt(time_bind),mse)
    
    all_melt2 <- rbind(m_melt,m_melt2,beta_melt,beta_melt2,time_melt)
    setting <- rep("2",sum(nrow(all_melt2)))
    all_melt2 <- cbind(all_melt2,setting)
    colnames(all_melt2)[2] <- "method"
  }
  load("~/mbDenoise/reproducibility/simulation/sim_lnm1.RData")
  {
    m_bind <- cbind(m_pr_y,m_pcoa_y,m_tsne_y,m_zifa,m_nbva,m_zip,m_zinb)
    m_bind2 <- cbind(m_pr_y2,m_pcoa_y2,m_tsne_y2,m_zifa2,m_nbva2,m_zip2,m_zinb2)
    
    beta_bind <- cbind(beta_pr_y,NaN,NaN,beta_zifa,beta_nbva,beta_zip,
                       beta_zinb)
    beta_bind2 <- cbind(beta_pr_y2,NaN,NaN,beta_zifa2,beta_nbva2,beta_zip2,
                        beta_zinb2)
    gamma_bind <- cbind(NaN,NaN,NaN,NaN,NaN,gamma_nbva_cov,NaN,gamma_zip_cov,NaN,
                        gamma_zinb_cov,NaN,gamma_zilnm_cov)
    
    time_bind <- log10(cbind(difftime_pr_y,difftime_pcoa_y,difftime_tsne_y,difftime_zifa,difftime_nbva,difftime_zip,difftime_zinb))
    
    colnames(m_bind) <- colnames(m_bind2) <-colnames(time_bind) <- colnames(beta_bind) <- colnames(beta_bind2) <-
      c("PCA","PCoA","tsne","ZIFA","GLLVM-NB","ZIPPCA-Poi","ZIPPCA-NB")
    
    # mse <- rep(paste("MSE1",expression(f_i),sep = ":"),sum(nrow(melt(m_bind))))
    # m_melt <- cbind(melt(m_bind),mse)
    # mse <- rep(paste("MSE2",expression(f_i),sep = ":"),sum(nrow(melt(m_bind2))))
    # m_melt2 <- cbind(melt(m_bind2),mse)
    # mse <- rep(paste("MSE1",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind))))
    # beta_melt <- cbind(melt(beta_bind),mse)
    # mse <- rep(paste("MSE2",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind2))))
    # beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("F_MSE1 ",sum(nrow(melt(m_bind))))
    m_melt <- cbind(melt(m_bind),mse)
    mse <- rep("F_MSE2",sum(nrow(melt(m_bind2))))
    m_melt2 <- cbind(melt(m_bind2),mse)
    mse <- rep("B_MSE1",sum(nrow(melt(beta_bind))))
    beta_melt <- cbind(melt(beta_bind),mse)
    mse <- rep("B_MSE2",sum(nrow(melt(beta_bind2))))
    beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("G_MSE1 ",sum(nrow(melt(gamma_bind))))
    g_melt <- cbind(melt(gamma_bind),mse)
    mse <- rep("time",sum(nrow(melt(time_bind))))
    time_melt <- cbind(melt(time_bind),mse)
    
    all_melt3 <- rbind(m_melt,m_melt2,beta_melt,beta_melt2,time_melt)
    setting <- rep("3",sum(nrow(all_melt3)))
    all_melt3 <- cbind(all_melt3,setting)
    colnames(all_melt3)[2] <- "method"
  }
  load("~/mbDenoise/reproducibility/simulation/sim_zifa1.RData")
  {
    m_bind <- cbind(m_pr_y,m_pcoa_y,m_tsne_y,m_zifa,m_nbva,m_zip,m_zinb)
    m_bind2 <- cbind(m_pr_y2,m_pcoa_y2,m_tsne_y2,m_zifa2,m_nbva2,m_zip2,m_zinb2)
    
    beta_bind <- cbind(beta_pr_y,NaN,NaN,beta_zifa,beta_nbva,beta_zip,
                       beta_zinb)
    beta_bind2 <- cbind(beta_pr_y2,NaN,NaN,beta_zifa2,beta_nbva2,beta_zip2,
                        beta_zinb2)
    gamma_bind <- cbind(NaN,NaN,NaN,NaN,NaN,gamma_nbva_cov,NaN,gamma_zip_cov,NaN,
                        gamma_zinb_cov,NaN,gamma_zilnm_cov)
    
    time_bind <- log10(cbind(difftime_pr_y,difftime_pcoa_y,difftime_tsne_y,difftime_zifa,difftime_nbva,difftime_zip,difftime_zinb))
    
    colnames(m_bind) <- colnames(m_bind2) <-colnames(time_bind) <- colnames(beta_bind) <- colnames(beta_bind2) <-
      c("PCA","PCoA","tsne","ZIFA","GLLVM-NB","ZIPPCA-Poi","ZIPPCA-NB")
    
    # mse <- rep(paste("MSE1",expression(f_i),sep = ":"),sum(nrow(melt(m_bind))))
    # m_melt <- cbind(melt(m_bind),mse)
    # mse <- rep(paste("MSE2",expression(f_i),sep = ":"),sum(nrow(melt(m_bind2))))
    # m_melt2 <- cbind(melt(m_bind2),mse)
    # mse <- rep(paste("MSE1",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind))))
    # beta_melt <- cbind(melt(beta_bind),mse)
    # mse <- rep(paste("MSE2",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind2))))
    # beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("F_MSE1 ",sum(nrow(melt(m_bind))))
    m_melt <- cbind(melt(m_bind),mse)
    mse <- rep("F_MSE2",sum(nrow(melt(m_bind2))))
    m_melt2 <- cbind(melt(m_bind2),mse)
    mse <- rep("B_MSE1",sum(nrow(melt(beta_bind))))
    beta_melt <- cbind(melt(beta_bind),mse)
    mse <- rep("B_MSE2",sum(nrow(melt(beta_bind2))))
    beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("G_MSE1 ",sum(nrow(melt(gamma_bind))))
    g_melt <- cbind(melt(gamma_bind),mse)
    mse <- rep("time",sum(nrow(melt(time_bind))))
    time_melt <- cbind(melt(time_bind),mse)
    
    all_melt7 <- rbind(m_melt,m_melt2,beta_melt,beta_melt2,time_melt)
    setting <- rep("7",sum(nrow(all_melt7)))
    all_melt7 <- cbind(all_melt7,setting)
    colnames(all_melt7)[2] <- "method"
  }
  load("~/mbDenoise/reproducibility/simulation/sim_gllvm1.RData")
  {
    m_bind <- cbind(m_pr_y,m_pcoa_y,m_tsne_y,m_zifa,m_nbva,m_zip,m_zinb)
    m_bind2 <- cbind(m_pr_y2,m_pcoa_y2,m_tsne_y2,m_zifa2,m_nbva2,m_zip2,m_zinb2)
    
    beta_bind <- cbind(beta_pr_y,NaN,NaN,beta_zifa,beta_nbva,beta_zip,
                       beta_zinb)
    beta_bind2 <- cbind(beta_pr_y2,NaN,NaN,beta_zifa2,beta_nbva2,beta_zip2,
                        beta_zinb2)
    gamma_bind <- cbind(NaN,NaN,NaN,NaN,NaN,gamma_nbva_cov,NaN,gamma_zip_cov,NaN,
                        gamma_zinb_cov,NaN,gamma_zilnm_cov)
    
    time_bind <- log10(cbind(difftime_pr_y,difftime_pcoa_y,difftime_tsne_y,difftime_zifa,difftime_nbva,difftime_zip,difftime_zinb))
    
    colnames(m_bind) <- colnames(m_bind2) <-colnames(time_bind) <- colnames(beta_bind) <- colnames(beta_bind2) <-
      c("PCA","PCoA","tsne","ZIFA","GLLVM-NB","ZIPPCA-Poi","ZIPPCA-NB")
    
    
    # mse <- rep(paste("MSE1",expression(f_i),sep = ":"),sum(nrow(melt(m_bind))))
    # m_melt <- cbind(melt(m_bind),mse)
    # mse <- rep(paste("MSE2",expression(f_i),sep = ":"),sum(nrow(melt(m_bind2))))
    # m_melt2 <- cbind(melt(m_bind2),mse)
    # mse <- rep(paste("MSE1",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind))))
    # beta_melt <- cbind(melt(beta_bind),mse)
    # mse <- rep(paste("MSE2",expression(beta_j),sep = ":"),sum(nrow(melt(beta_bind2))))
    # beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("F_MSE1 ",sum(nrow(melt(m_bind))))
    m_melt <- cbind(melt(m_bind),mse)
    mse <- rep("F_MSE2",sum(nrow(melt(m_bind2))))
    m_melt2 <- cbind(melt(m_bind2),mse)
    mse <- rep("B_MSE1",sum(nrow(melt(beta_bind))))
    beta_melt <- cbind(melt(beta_bind),mse)
    mse <- rep("B_MSE2",sum(nrow(melt(beta_bind2))))
    beta_melt2 <- cbind(melt(beta_bind2),mse)
    mse <- rep("G_MSE1 ",sum(nrow(melt(gamma_bind))))
    g_melt <- cbind(melt(gamma_bind),mse)
    mse <- rep("time",sum(nrow(melt(time_bind))))
    time_melt <- cbind(melt(time_bind),mse)
    
    all_melt8 <- rbind(m_melt,m_melt2,beta_melt,beta_melt2,time_melt)
    setting <- rep("8",sum(nrow(all_melt8)))
    all_melt8 <- cbind(all_melt8,setting)
    colnames(all_melt8)[2] <- "method"
  }
  
}

df <- rbind(all_melt,all_melt2,all_melt3,all_melt7,all_melt8)
df$mse <- factor(df$mse,labels = c("F1","F2","B1","B2","time"))
df$setting <- factor(df$setting,labels = c("M1","M2","M3","M4","M5"))
df$method <- factor(df$method,labels = c("PPCA-NB","PCA","PCoA","t-SNE","ZIFA","mbDenoise-zinb",
                                         "mbDenoise-zip"))
colourCount = length(unique(df$method))

ggplot(data=df, aes(x=factor(method,
                             level=c("PCA","PCoA","t-SNE","ZIFA","PPCA-NB","mbDenoise-zip","mbDenoise-zinb")),
                    y=value,fill=factor(method,level=c("PCA","PCoA","t-SNE","ZIFA","PPCA-NB","mbDenoise-zip","mbDenoise-zinb")),
                    colour=factor(method,level=c("PCA","PCoA","t-SNE","ZIFA","PPCA-NB","mbDenoise-zip","mbDenoise-zinb"))))+
  geom_boxplot(outlier.alpha = 0.18,outlier.size = 0.8,lwd=0.2)+
  stat_boxplot(geom = "errorbar",width=0.35)+ 
  scale_fill_manual(values =  c(colorRampPalette(brewer.pal(7, "Paired"))(colourCount)[-6],colorRampPalette(brewer.pal(7, "Paired"))(colourCount)[6]))+
  scale_colour_manual(values =  c(colorRampPalette(brewer.pal(7, "Paired"))(colourCount)[-6],colorRampPalette(brewer.pal(7, "Paired"))(colourCount)[6]))+
  guides(fill = guide_legend(title = "Method"),color=guide_legend(title = "Method"))+
  labs(x = "",y = "", title = "")+
  theme_bw()+theme(strip.text =element_text(size = 10),axis.text=element_text(size = 6),legend.title=element_text(size = 8),
                   legend.text = element_text(size = 7),axis.ticks.x =element_blank()) +
  theme(axis.text.x =element_blank()) +
  facet_grid(mse ~ setting,scales="free")


