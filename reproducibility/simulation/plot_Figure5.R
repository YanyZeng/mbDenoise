library(reshape)
library(ggplot2)
library(ggpubr)
library(grid)
library(cowplot)

ref_p <-100

##M7
{
  load("~/mbDenoise/reproducibility/simulation/sim_nb1.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_1 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_1 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_1 <- 2*(power_poi_1*(1-FDR_poi_1))/(power_poi_1+(1-FDR_poi_1))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_1 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_1 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_1 <- 2*(power_ppca_1*(1-FDR_ppca_1))/(power_ppca_1+(1-FDR_ppca_1))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_1 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_1 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_1 <- 2*(power_nb_1*(1-FDR_nb_1))/(power_nb_1+(1-FDR_nb_1))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_1 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_1 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_1 <- 2*(power_ln_1*(1-FDR_ln_1))/(power_ln_1+(1-FDR_ln_1))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_1 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_1 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_1 <- 2*(power_poi_cov_1*(1-FDR_poi_cov_1))/(power_poi_cov_1+(1-FDR_poi_cov_1))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_1 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_1 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_1 <- 2*(power_nb_cov_1*(1-FDR_nb_cov_1))/(power_nb_cov_1+(1-FDR_nb_cov_1))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_1 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_1 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_1 <- 2*(power_ln_cov_1*(1-FDR_ln_cov_1))/(power_ln_cov_1+(1-FDR_ln_cov_1))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_1 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_1 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_1 <- 2*(power_edger_1*(1-FDR_edger_1))/(power_edger_1+(1-FDR_edger_1))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_1 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_1 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_1 <- 2*(power_deseq_1*(1-FDR_deseq_1))/(power_deseq_1+(1-FDR_deseq_1))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_1 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_1 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_1 <- 2*(power_meta_1*(1-FDR_meta_1))/(power_meta_1+(1-FDR_meta_1))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_12 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_12 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_12 <- 2*(power_meta_12*(1-FDR_meta_12))/(power_meta_12+(1-FDR_meta_12))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_1 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_1 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_1 <- 2*(power_real_1*(1-FDR_real_1))/(power_real_1+(1-FDR_real_1))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_12 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_12 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_12 <- 2*(power_real_12*(1-FDR_real_12))/(power_real_12+(1-FDR_real_12))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_1 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_1 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_1 <- 2*(power_mbi_1*(1-FDR_mbi_1))/(power_mbi_1+(1-FDR_mbi_1))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_1 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_1 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_1 <- 2*(power_sav_1*(1-FDR_sav_1))/(power_sav_1+(1-FDR_sav_1))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_nb1_1.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_2 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_2 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_2 <- 2*(power_poi_2*(1-FDR_poi_2))/(power_poi_2+(1-FDR_poi_2))
    
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_2 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_2 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_2 <- 2*(power_ppca_2*(1-FDR_ppca_2))/(power_ppca_2+(1-FDR_ppca_2))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_2 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_2 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_2 <- 2*(power_nb_2*(1-FDR_nb_2))/(power_nb_2+(1-FDR_nb_2))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_2 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_2 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_2 <- 2*(power_ln_2*(1-FDR_ln_2))/(power_ln_2+(1-FDR_ln_2))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_2 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_2 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_2 <- 2*(power_poi_cov_2*(1-FDR_poi_cov_2))/(power_poi_cov_2+(1-FDR_poi_cov_2))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_2 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_2 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_2 <- 2*(power_nb_cov_2*(1-FDR_nb_cov_2))/(power_nb_cov_2+(1-FDR_nb_cov_2))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_2 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_2 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_2 <- 2*(power_ln_cov_2*(1-FDR_ln_cov_2))/(power_ln_cov_2+(1-FDR_ln_cov_2))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_2 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_2 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_2 <- 2*(power_edger_2*(1-FDR_edger_2))/(power_edger_2+(1-FDR_edger_2))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_2 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_2 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_2 <- 2*(power_deseq_2*(1-FDR_deseq_2))/(power_deseq_2+(1-FDR_deseq_2))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_2 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_2 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_2 <- 2*(power_meta_2*(1-FDR_meta_2))/(power_meta_2+(1-FDR_meta_2))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_22 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_22 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_22 <- 2*(power_meta_22*(1-FDR_meta_22))/(power_meta_22+(1-FDR_meta_22))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_2 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_2 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_2 <- 2*(power_real_2*(1-FDR_real_2))/(power_real_2+(1-FDR_real_2))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_22 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_22 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_22 <- 2*(power_real_22*(1-FDR_real_22))/(power_real_22+(1-FDR_real_22))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_2 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_2 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_2 <- 2*(power_mbi_2*(1-FDR_mbi_2))/(power_mbi_2+(1-FDR_mbi_2))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_2 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_2 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_2 <- 2*(power_sav_2*(1-FDR_sav_2))/(power_sav_2+(1-FDR_sav_2))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_nb1_2.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_3 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_3 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_3 <- 2*(power_poi_3*(1-FDR_poi_3))/(power_poi_3+(1-FDR_poi_3))
    
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_3 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_3 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_3 <- 2*(power_ppca_3*(1-FDR_ppca_3))/(power_ppca_3+(1-FDR_ppca_3))
    
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_3 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_3 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_3 <- 2*(power_nb_3*(1-FDR_nb_3))/(power_nb_3+(1-FDR_nb_3))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_3 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_3 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_3 <- 2*(power_ln_3*(1-FDR_ln_3))/(power_ln_3+(1-FDR_ln_3))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_3 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_3 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_3 <- 2*(power_poi_cov_3*(1-FDR_poi_cov_3))/(power_poi_cov_3+(1-FDR_poi_cov_3))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_3 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_3 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_3 <- 2*(power_nb_cov_3*(1-FDR_nb_cov_3))/(power_nb_cov_3+(1-FDR_nb_cov_3))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_3 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_3 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_3 <- 2*(power_ln_cov_3*(1-FDR_ln_cov_3))/(power_ln_cov_3+(1-FDR_ln_cov_3))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_3 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_3 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_3 <- 2*(power_edger_3*(1-FDR_edger_3))/(power_edger_3+(1-FDR_edger_3))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_3 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_3 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_3 <- 2*(power_deseq_3*(1-FDR_deseq_3))/(power_deseq_3+(1-FDR_deseq_3))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_3 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_3 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_3 <- 2*(power_meta_3*(1-FDR_meta_3))/(power_meta_3+(1-FDR_meta_3))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_32 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_32 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_32 <- 2*(power_meta_32*(1-FDR_meta_32))/(power_meta_32+(1-FDR_meta_32))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_3 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_3 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_3 <- 2*(power_real_3*(1-FDR_real_3))/(power_real_3+(1-FDR_real_3))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_32 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_32 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_32 <- 2*(power_real_32*(1-FDR_real_32))/(power_real_32+(1-FDR_real_32))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_3 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_3 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_3 <- 2*(power_mbi_3*(1-FDR_mbi_3))/(power_mbi_3+(1-FDR_mbi_3))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_3 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_3 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_3 <- 2*(power_sav_3*(1-FDR_sav_3))/(power_sav_3+(1-FDR_sav_3))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_nb1_3.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_4 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_4 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_4 <- 2*(power_poi_4*(1-FDR_poi_4))/(power_poi_4+(1-FDR_poi_4))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_4 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_4 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_4 <- 2*(power_ppca_4*(1-FDR_ppca_4))/(power_ppca_4+(1-FDR_ppca_4))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_4 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_4 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_4 <- 2*(power_nb_4*(1-FDR_nb_4))/(power_nb_4+(1-FDR_nb_4))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_4 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_4 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_4 <- 2*(power_ln_4*(1-FDR_ln_4))/(power_ln_4+(1-FDR_ln_4))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_4 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_4 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_4 <- 2*(power_poi_cov_4*(1-FDR_poi_cov_4))/(power_poi_cov_4+(1-FDR_poi_cov_4))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_4 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_4 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_4 <- 2*(power_nb_cov_4*(1-FDR_nb_cov_4))/(power_nb_cov_4+(1-FDR_nb_cov_4))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_4 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_4 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_4 <- 2*(power_ln_cov_4*(1-FDR_ln_cov_4))/(power_ln_cov_4+(1-FDR_ln_cov_4))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_4 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_4 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_4 <- 2*(power_edger_4*(1-FDR_edger_4))/(power_edger_4+(1-FDR_edger_4))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_4 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_4 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_4 <- 2*(power_deseq_4*(1-FDR_deseq_4))/(power_deseq_4+(1-FDR_deseq_4))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_4 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_4 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_4 <- 2*(power_meta_4*(1-FDR_meta_4))/(power_meta_4+(1-FDR_meta_4))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_42 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_42 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_42 <- 2*(power_meta_42*(1-FDR_meta_42))/(power_meta_42+(1-FDR_meta_42))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_4 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_4 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_4 <- 2*(power_real_4*(1-FDR_real_4))/(power_real_4+(1-FDR_real_4))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_42 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_42 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_42 <- 2*(power_real_42*(1-FDR_real_42))/(power_real_42+(1-FDR_real_42))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_4 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_4 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_4 <- 2*(power_mbi_4*(1-FDR_mbi_4))/(power_mbi_4+(1-FDR_mbi_4))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_4 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_4 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_4 <- 2*(power_sav_4*(1-FDR_sav_4))/(power_sav_4+(1-FDR_sav_4))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_nb1_4.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_5 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_5 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_5 <- 2*(power_poi_5*(1-FDR_poi_5))/(power_poi_5+(1-FDR_poi_5))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_5 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_5 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_5 <- 2*(power_ppca_5*(1-FDR_ppca_5))/(power_ppca_5+(1-FDR_ppca_5))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_5 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_5 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_5 <- 2*(power_nb_5*(1-FDR_nb_5))/(power_nb_5+(1-FDR_nb_5))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_5 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_5 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_5 <- 2*(power_ln_5*(1-FDR_ln_5))/(power_ln_5+(1-FDR_ln_5))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_5 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_5 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_5 <- 2*(power_poi_cov_5*(1-FDR_poi_cov_5))/(power_poi_cov_5+(1-FDR_poi_cov_5))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_5 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_5 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_5 <- 2*(power_nb_cov_5*(1-FDR_nb_cov_5))/(power_nb_cov_5+(1-FDR_nb_cov_5))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_5 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_5 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_5 <- 2*(power_ln_cov_5*(1-FDR_ln_cov_5))/(power_ln_cov_5+(1-FDR_ln_cov_5))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_5 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_5 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_5 <- 2*(power_edger_5*(1-FDR_edger_5))/(power_edger_5+(1-FDR_edger_5))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_5 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_5 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_5 <- 2*(power_deseq_5*(1-FDR_deseq_5))/(power_deseq_5+(1-FDR_deseq_5))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_5 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_5 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_5 <- 2*(power_meta_5*(1-FDR_meta_5))/(power_meta_5+(1-FDR_meta_5))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_52 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_52 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_52 <- 2*(power_meta_52*(1-FDR_meta_52))/(power_meta_52+(1-FDR_meta_52))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_5 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_5 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_5 <- 2*(power_real_5*(1-FDR_real_5))/(power_real_5+(1-FDR_real_5))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_52 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_52 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_52 <- 2*(power_real_52*(1-FDR_real_52))/(power_real_52+(1-FDR_real_52))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_5 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_5 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_5 <- 2*(power_mbi_5*(1-FDR_mbi_5))/(power_mbi_5+(1-FDR_mbi_5))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_5 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_5 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_5 <- 2*(power_sav_5*(1-FDR_sav_5))/(power_sav_5+(1-FDR_sav_5))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_nb1_5.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_6 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_6 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_6 <- 2*(power_poi_6*(1-FDR_poi_6))/(power_poi_6+(1-FDR_poi_6))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_6 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_6 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_6 <- 2*(power_ppca_6*(1-FDR_ppca_6))/(power_ppca_6+(1-FDR_ppca_6))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_6 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_6 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_6 <- 2*(power_nb_6*(1-FDR_nb_6))/(power_nb_6+(1-FDR_nb_6))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_6 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_6 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_6 <- 2*(power_ln_6*(1-FDR_ln_6))/(power_ln_6+(1-FDR_ln_6))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_6 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_6 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_6 <- 2*(power_poi_cov_6*(1-FDR_poi_cov_6))/(power_poi_cov_6+(1-FDR_poi_cov_6))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_6 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_6 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_6 <- 2*(power_nb_cov_6*(1-FDR_nb_cov_6))/(power_nb_cov_6+(1-FDR_nb_cov_6))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_6 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_6 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_6 <- 2*(power_ln_cov_6*(1-FDR_ln_cov_6))/(power_ln_cov_6+(1-FDR_ln_cov_6))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_6 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_6 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_6 <- 2*(power_edger_6*(1-FDR_edger_6))/(power_edger_6+(1-FDR_edger_6))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_6 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_6 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_6 <- 2*(power_deseq_6*(1-FDR_deseq_6))/(power_deseq_6+(1-FDR_deseq_6))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_6 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_6 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_6 <- 2*(power_meta_6*(1-FDR_meta_6))/(power_meta_6+(1-FDR_meta_6))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_62 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_62 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_62 <- 2*(power_meta_62*(1-FDR_meta_62))/(power_meta_62+(1-FDR_meta_62))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_6 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_6 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_6 <- 2*(power_real_6*(1-FDR_real_6))/(power_real_6+(1-FDR_real_6))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_62 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_62 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_62 <- 2*(power_real_62*(1-FDR_real_62))/(power_real_62+(1-FDR_real_62))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_6 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_6 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_6 <- 2*(power_mbi_6*(1-FDR_mbi_6))/(power_mbi_6+(1-FDR_mbi_6))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_6 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_6 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_6 <- 2*(power_sav_6*(1-FDR_sav_6))/(power_sav_6+(1-FDR_sav_6))
    
    
  }
  
  {
    Recall_A <- c(power_real_1,power_real_2,power_real_3,power_real_4,power_real_5,power_real_6)
    Recall_B <- c(power_real_12,power_real_22,power_real_32,power_real_42,power_real_52,power_real_62)
    Recall_C <- c(power_nb_1,power_nb_2,power_nb_3,power_nb_4,power_nb_5,power_nb_6)
    Recall_D <- c(power_poi_1,power_poi_2,power_poi_3,power_poi_4,power_poi_5,power_poi_6)
    Recall_E <- c(power_ppca_1,power_ppca_2,power_ppca_3,power_ppca_4,power_ppca_5,power_ppca_6)
    Recall_F <- c(power_nb_cov_1,power_nb_cov_2,power_nb_cov_3,power_nb_cov_4,power_nb_cov_5,power_nb_cov_6)
    Recall_G <- c(power_poi_cov_1,power_poi_cov_2,power_poi_cov_3,power_poi_cov_4,power_poi_cov_5,power_poi_cov_6)
    Recall_H <- c(power_ln_cov_1,power_ln_cov_2,power_ln_cov_3,power_ln_cov_4,power_ln_cov_5,power_ln_cov_6)
    Recall_I <- c(power_deseq_1,power_deseq_2,power_deseq_3,power_deseq_4,power_deseq_5,power_deseq_6)
    Recall_J <- c(power_edger_1,power_edger_2,power_edger_3,power_edger_4,power_edger_5,power_edger_6)
    Recall_K <- c(power_meta_1,power_meta_2,power_meta_3,power_meta_4,power_meta_5,power_meta_6)
    Recall_L <- c(power_meta_12,power_meta_22,power_meta_32,power_meta_42,power_meta_52,power_meta_62)
    Recall_M <- c(power_sav_1,power_sav_2,power_sav_3,power_sav_4,power_sav_5,power_sav_6)
    Recall_N <- c(power_mbi_1,power_mbi_2,power_mbi_3,power_mbi_4,power_mbi_5,power_mbi_6)
  }
  {
    Precision_A <- c(1-FDR_real_1,1-FDR_real_2,1-FDR_real_3,1-FDR_real_4,1-FDR_real_5,1-FDR_real_6)
    Precision_B <- c(1-FDR_real_12,1-FDR_real_22,1-FDR_real_32,1-FDR_real_42,1-FDR_real_52,1-FDR_real_62)
    Precision_C <- c(c(1-FDR_nb_1,1-FDR_nb_2,1-FDR_nb_3,1-FDR_nb_4,1-FDR_nb_5,1-FDR_nb_6))
    Precision_D <- c(1-FDR_poi_1,1-FDR_poi_2,1-FDR_poi_3,1-FDR_poi_4,1-FDR_poi_5,1-FDR_poi_6)
    Precision_E <- c(1-FDR_ppca_1,1-FDR_ppca_2,1-FDR_ppca_3,1-FDR_ppca_4,1-FDR_ppca_5,1-FDR_ppca_6)
    Precision_F <- c(c(1-FDR_nb_cov_1,1-FDR_nb_cov_2,1-FDR_nb_cov_3,1-FDR_nb_cov_4,1-FDR_nb_cov_5,1-FDR_nb_cov_6))
    Precision_G <- c(1-FDR_poi_cov_1,1-FDR_poi_cov_2,1-FDR_poi_cov_3,1-FDR_poi_cov_4,1-FDR_poi_cov_5,1-FDR_poi_cov_6)
    Precision_H <- c(1-FDR_ln_cov_1,1-FDR_ln_cov_2,1-FDR_ln_cov_3,1-FDR_ln_cov_4,1-FDR_ln_cov_5,1-FDR_ln_cov_6)
    Precision_I <- c(1-FDR_deseq_1,1-FDR_deseq_2,1-FDR_deseq_3,1-FDR_deseq_4,1-FDR_deseq_5,1-FDR_deseq_6)
    Precision_J <- c(1-FDR_edger_1,1-FDR_edger_2,1-FDR_edger_3,1-FDR_edger_4,1-FDR_edger_5,1-FDR_edger_6)
    Precision_K <- c(1-FDR_meta_1,1-FDR_meta_2,1-FDR_meta_3,1-FDR_meta_4,1-FDR_meta_5,1-FDR_meta_6)
    Precision_L <- c(1-FDR_meta_12,1-FDR_meta_22,1-FDR_meta_32,1-FDR_meta_42,1-FDR_meta_52,1-FDR_meta_62)
    Precision_M <- c(1-FDR_sav_1,1-FDR_sav_2,1-FDR_sav_3,1-FDR_sav_4,1-FDR_sav_5,1-FDR_sav_6)
    Precision_N <- c(1-FDR_mbi_1,1-FDR_mbi_2,1-FDR_mbi_3,1-FDR_mbi_4,1-FDR_mbi_5,1-FDR_mbi_6)
    
  }
  {
    F1_A <- c(score_real_1,score_real_2,score_real_3,score_real_4,score_real_5,score_real_6)
    F1_B <- c(score_real_12,score_real_22,score_real_32,score_real_42,score_real_52,score_real_62)
    F1_C <- c(score_nb_1,score_nb_2,score_nb_3,score_nb_4,score_nb_5,score_nb_6)
    F1_D <- c(score_poi_1,score_poi_2,score_poi_3,score_poi_4,score_poi_5,score_poi_6)
    F1_E <- c(score_ppca_1,score_ppca_2,score_ppca_3,score_ppca_4,score_ppca_5,score_ppca_6)
    F1_F <- c(score_nb_cov_1,score_nb_cov_2,score_nb_cov_3,score_nb_cov_4,score_nb_cov_5,score_nb_cov_6)
    F1_G <- c(score_poi_cov_1,score_poi_cov_2,score_poi_cov_3,score_poi_cov_4,score_poi_cov_5,score_poi_cov_6)
    F1_H <- c(score_ln_cov_1,score_ln_cov_2,score_ln_cov_3,score_ln_cov_4,score_ln_cov_5,score_ln_cov_6)
    F1_I <- c(score_deseq_1,score_deseq_2,score_deseq_3,score_deseq_4,score_deseq_5,score_deseq_6)
    F1_J <- c(score_edger_1,score_edger_2,score_edger_3,score_edger_4,score_edger_5,score_edger_6)
    F1_K <- c(score_meta_1,score_meta_2,score_meta_3,score_meta_4,score_meta_5,score_meta_6)
    F1_L <- c(score_meta_12,score_meta_22,score_meta_32,score_meta_42,score_meta_52,score_meta_62)
    F1_M <- c(score_sav_1,score_sav_2,score_sav_3,score_sav_4,score_sav_5,score_sav_6)
    F1_N <- c(score_mbi_1,score_mbi_2,score_mbi_3,score_mbi_4,score_mbi_5,score_mbi_6)
  }
  
  df <- cbind(melt(c(Recall_B,Recall_E,Recall_F,Recall_G,Recall_I,Recall_J,Recall_K,
                     Recall_M,Recall_N)),
              melt(c(Precision_B,Precision_E,Precision_F,Precision_G,Precision_I,Precision_J,Precision_K,
                     Precision_M,Precision_N)),
              melt(c(F1_B,F1_E,F1_F,F1_G,F1_I,F1_J,F1_K,
                     F1_M,F1_N)),
              melt(c(rep("No imputation_t",6),rep("PPCA-NB-cov",6),
                     rep("mbDenoise-ZINB-cov_t",6),rep("mbDenoise-ZIP-cov_t",6),rep("DESeq2",6),rep("edgeR",6),
                     rep("metagenomeSeq:fitFeaturemodel",6),rep("SAVER_t",6),rep("mbImpute_t",6))))
  colnames(df)[4] <- "method";colnames(df)[1] <- "Recall";colnames(df)[2] <- "Precision";colnames(df)[3] <- "F1 score";
  x <- rep(c(0:5),9*3)
  df1 <- cbind(x,melt(df))
  setting <- rep("1",sum(nrow(df1)))
  df1 <- cbind(df1,setting)
}
##M8
{
  load("~/mbDenoise/reproducibility/simulation/sim_poi1.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_1 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_1 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_1 <- 2*(power_poi_1*(1-FDR_poi_1))/(power_poi_1+(1-FDR_poi_1))
    
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_1 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_1 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_1 <- 2*(power_ppca_1*(1-FDR_ppca_1))/(power_ppca_1+(1-FDR_ppca_1))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_1 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_1 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_1 <- 2*(power_nb_1*(1-FDR_nb_1))/(power_nb_1+(1-FDR_nb_1))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_1 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_1 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_1 <- 2*(power_ln_1*(1-FDR_ln_1))/(power_ln_1+(1-FDR_ln_1))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_1 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_1 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_1 <- 2*(power_poi_cov_1*(1-FDR_poi_cov_1))/(power_poi_cov_1+(1-FDR_poi_cov_1))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_1 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_1 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_1 <- 2*(power_nb_cov_1*(1-FDR_nb_cov_1))/(power_nb_cov_1+(1-FDR_nb_cov_1))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_1 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_1 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_1 <- 2*(power_ln_cov_1*(1-FDR_ln_cov_1))/(power_ln_cov_1+(1-FDR_ln_cov_1))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_1 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_1 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_1 <- 2*(power_edger_1*(1-FDR_edger_1))/(power_edger_1+(1-FDR_edger_1))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_1 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_1 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_1 <- 2*(power_deseq_1*(1-FDR_deseq_1))/(power_deseq_1+(1-FDR_deseq_1))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_1 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_1 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_1 <- 2*(power_meta_1*(1-FDR_meta_1))/(power_meta_1+(1-FDR_meta_1))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_12 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_12 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_12 <- 2*(power_meta_12*(1-FDR_meta_12))/(power_meta_12+(1-FDR_meta_12))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_1 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_1 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_1 <- 2*(power_real_1*(1-FDR_real_1))/(power_real_1+(1-FDR_real_1))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_12 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_12 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_12 <- 2*(power_real_12*(1-FDR_real_12))/(power_real_12+(1-FDR_real_12))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_1 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_1 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_1 <- 2*(power_mbi_1*(1-FDR_mbi_1))/(power_mbi_1+(1-FDR_mbi_1))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_1 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_1 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_1 <- 2*(power_sav_1*(1-FDR_sav_1))/(power_sav_1+(1-FDR_sav_1))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_poi1_1.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_2 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_2 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_2 <- 2*(power_poi_2*(1-FDR_poi_2))/(power_poi_2+(1-FDR_poi_2))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_2 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_2 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_2 <- 2*(power_ppca_2*(1-FDR_ppca_2))/(power_ppca_2+(1-FDR_ppca_2))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_2 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_2 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_2 <- 2*(power_nb_2*(1-FDR_nb_2))/(power_nb_2+(1-FDR_nb_2))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_2 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_2 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_2 <- 2*(power_ln_2*(1-FDR_ln_2))/(power_ln_2+(1-FDR_ln_2))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_2 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_2 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_2 <- 2*(power_poi_cov_2*(1-FDR_poi_cov_2))/(power_poi_cov_2+(1-FDR_poi_cov_2))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_2 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_2 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_2 <- 2*(power_nb_cov_2*(1-FDR_nb_cov_2))/(power_nb_cov_2+(1-FDR_nb_cov_2))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_2 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_2 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_2 <- 2*(power_ln_cov_2*(1-FDR_ln_cov_2))/(power_ln_cov_2+(1-FDR_ln_cov_2))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_2 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_2 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_2 <- 2*(power_edger_2*(1-FDR_edger_2))/(power_edger_2+(1-FDR_edger_2))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_2 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_2 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_2 <- 2*(power_deseq_2*(1-FDR_deseq_2))/(power_deseq_2+(1-FDR_deseq_2))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_2 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_2 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_2 <- 2*(power_meta_2*(1-FDR_meta_2))/(power_meta_2+(1-FDR_meta_2))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_22 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_22 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_22 <- 2*(power_meta_22*(1-FDR_meta_22))/(power_meta_22+(1-FDR_meta_22))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_2 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_2 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_2 <- 2*(power_real_2*(1-FDR_real_2))/(power_real_2+(1-FDR_real_2))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_22 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_22 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_22 <- 2*(power_real_22*(1-FDR_real_22))/(power_real_22+(1-FDR_real_22))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_2 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_2 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_2 <- 2*(power_mbi_2*(1-FDR_mbi_2))/(power_mbi_2+(1-FDR_mbi_2))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_2 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_2 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_2 <- 2*(power_sav_2*(1-FDR_sav_2))/(power_sav_2+(1-FDR_sav_2))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_poi1_2.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_3 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_3 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_3 <- 2*(power_poi_3*(1-FDR_poi_3))/(power_poi_3+(1-FDR_poi_3))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_3 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_3 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_3 <- 2*(power_ppca_3*(1-FDR_ppca_3))/(power_ppca_3+(1-FDR_ppca_3))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_3 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_3 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_3 <- 2*(power_nb_3*(1-FDR_nb_3))/(power_nb_3+(1-FDR_nb_3))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_3 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_3 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_3 <- 2*(power_ln_3*(1-FDR_ln_3))/(power_ln_3+(1-FDR_ln_3))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_3 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_3 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_3 <- 2*(power_poi_cov_3*(1-FDR_poi_cov_3))/(power_poi_cov_3+(1-FDR_poi_cov_3))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_3 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_3 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_3 <- 2*(power_nb_cov_3*(1-FDR_nb_cov_3))/(power_nb_cov_3+(1-FDR_nb_cov_3))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_3 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_3 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_3 <- 2*(power_ln_cov_3*(1-FDR_ln_cov_3))/(power_ln_cov_3+(1-FDR_ln_cov_3))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_3 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_3 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_3 <- 2*(power_edger_3*(1-FDR_edger_3))/(power_edger_3+(1-FDR_edger_3))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_3 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_3 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_3 <- 2*(power_deseq_3*(1-FDR_deseq_3))/(power_deseq_3+(1-FDR_deseq_3))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_3 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_3 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_3 <- 2*(power_meta_3*(1-FDR_meta_3))/(power_meta_3+(1-FDR_meta_3))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_32 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_32 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_32 <- 2*(power_meta_32*(1-FDR_meta_32))/(power_meta_32+(1-FDR_meta_32))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_3 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_3 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_3 <- 2*(power_real_3*(1-FDR_real_3))/(power_real_3+(1-FDR_real_3))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_32 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_32 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_32 <- 2*(power_real_32*(1-FDR_real_32))/(power_real_32+(1-FDR_real_32))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_3 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_3 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_3 <- 2*(power_mbi_3*(1-FDR_mbi_3))/(power_mbi_3+(1-FDR_mbi_3))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_3 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_3 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_3 <- 2*(power_sav_3*(1-FDR_sav_3))/(power_sav_3+(1-FDR_sav_3))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_poi1_3.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_4 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_4 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_4 <- 2*(power_poi_4*(1-FDR_poi_4))/(power_poi_4+(1-FDR_poi_4))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_4 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_4 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_4 <- 2*(power_ppca_4*(1-FDR_ppca_4))/(power_ppca_4+(1-FDR_ppca_4))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_4 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_4 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_4 <- 2*(power_nb_4*(1-FDR_nb_4))/(power_nb_4+(1-FDR_nb_4))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_4 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_4 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_4 <- 2*(power_ln_4*(1-FDR_ln_4))/(power_ln_4+(1-FDR_ln_4))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_4 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_4 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_4 <- 2*(power_poi_cov_4*(1-FDR_poi_cov_4))/(power_poi_cov_4+(1-FDR_poi_cov_4))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_4 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_4 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_4 <- 2*(power_nb_cov_4*(1-FDR_nb_cov_4))/(power_nb_cov_4+(1-FDR_nb_cov_4))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_4 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_4 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_4 <- 2*(power_ln_cov_4*(1-FDR_ln_cov_4))/(power_ln_cov_4+(1-FDR_ln_cov_4))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_4 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_4 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_4 <- 2*(power_edger_4*(1-FDR_edger_4))/(power_edger_4+(1-FDR_edger_4))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_4 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_4 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_4 <- 2*(power_deseq_4*(1-FDR_deseq_4))/(power_deseq_4+(1-FDR_deseq_4))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_4 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_4 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_4 <- 2*(power_meta_4*(1-FDR_meta_4))/(power_meta_4+(1-FDR_meta_4))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_42 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_42 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_42 <- 2*(power_meta_42*(1-FDR_meta_42))/(power_meta_42+(1-FDR_meta_42))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_4 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_4 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_4 <- 2*(power_real_4*(1-FDR_real_4))/(power_real_4+(1-FDR_real_4))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_42 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_42 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_42 <- 2*(power_real_42*(1-FDR_real_42))/(power_real_42+(1-FDR_real_42))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_4 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_4 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_4 <- 2*(power_mbi_4*(1-FDR_mbi_4))/(power_mbi_4+(1-FDR_mbi_4))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_4 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_4 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_4 <- 2*(power_sav_4*(1-FDR_sav_4))/(power_sav_4+(1-FDR_sav_4))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_poi1_4.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_5 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_5 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_5 <- 2*(power_poi_5*(1-FDR_poi_5))/(power_poi_5+(1-FDR_poi_5))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_5 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_5 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_5 <- 2*(power_ppca_5*(1-FDR_ppca_5))/(power_ppca_5+(1-FDR_ppca_5))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_5 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_5 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_5 <- 2*(power_nb_5*(1-FDR_nb_5))/(power_nb_5+(1-FDR_nb_5))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_5 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_5 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_5 <- 2*(power_ln_5*(1-FDR_ln_5))/(power_ln_5+(1-FDR_ln_5))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_5 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_5 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_5 <- 2*(power_poi_cov_5*(1-FDR_poi_cov_5))/(power_poi_cov_5+(1-FDR_poi_cov_5))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_5 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_5 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_5 <- 2*(power_nb_cov_5*(1-FDR_nb_cov_5))/(power_nb_cov_5+(1-FDR_nb_cov_5))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_5 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_5 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_5 <- 2*(power_ln_cov_5*(1-FDR_ln_cov_5))/(power_ln_cov_5+(1-FDR_ln_cov_5))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_5 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_5 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_5 <- 2*(power_edger_5*(1-FDR_edger_5))/(power_edger_5+(1-FDR_edger_5))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_5 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_5 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_5 <- 2*(power_deseq_5*(1-FDR_deseq_5))/(power_deseq_5+(1-FDR_deseq_5))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_5 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_5 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_5 <- 2*(power_meta_5*(1-FDR_meta_5))/(power_meta_5+(1-FDR_meta_5))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_52 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_52 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_52 <- 2*(power_meta_52*(1-FDR_meta_52))/(power_meta_52+(1-FDR_meta_52))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_5 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_5 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_5 <- 2*(power_real_5*(1-FDR_real_5))/(power_real_5+(1-FDR_real_5))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_52 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_52 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_52 <- 2*(power_real_52*(1-FDR_real_52))/(power_real_52+(1-FDR_real_52))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_5 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_5 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_5 <- 2*(power_mbi_5*(1-FDR_mbi_5))/(power_mbi_5+(1-FDR_mbi_5))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_5 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_5 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_5 <- 2*(power_sav_5*(1-FDR_sav_5))/(power_sav_5+(1-FDR_sav_5))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_poi1_5.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_6 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_6 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_6 <- 2*(power_poi_6*(1-FDR_poi_6))/(power_poi_6+(1-FDR_poi_6))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_6 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_6 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_6 <- 2*(power_ppca_6*(1-FDR_ppca_6))/(power_ppca_6+(1-FDR_ppca_6))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_6 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_6 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_6 <- 2*(power_nb_6*(1-FDR_nb_6))/(power_nb_6+(1-FDR_nb_6))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_6 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_6 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_6 <- 2*(power_ln_6*(1-FDR_ln_6))/(power_ln_6+(1-FDR_ln_6))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_6 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_6 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_6 <- 2*(power_poi_cov_6*(1-FDR_poi_cov_6))/(power_poi_cov_6+(1-FDR_poi_cov_6))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_6 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_6 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_6 <- 2*(power_nb_cov_6*(1-FDR_nb_cov_6))/(power_nb_cov_6+(1-FDR_nb_cov_6))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_6 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_6 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_6 <- 2*(power_ln_cov_6*(1-FDR_ln_cov_6))/(power_ln_cov_6+(1-FDR_ln_cov_6))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_6 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_6 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_6 <- 2*(power_edger_6*(1-FDR_edger_6))/(power_edger_6+(1-FDR_edger_6))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_6 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_6 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_6 <- 2*(power_deseq_6*(1-FDR_deseq_6))/(power_deseq_6+(1-FDR_deseq_6))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_6 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_6 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_6 <- 2*(power_meta_6*(1-FDR_meta_6))/(power_meta_6+(1-FDR_meta_6))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_62 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_62 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_62 <- 2*(power_meta_62*(1-FDR_meta_62))/(power_meta_62+(1-FDR_meta_62))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_6 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_6 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_6 <- 2*(power_real_6*(1-FDR_real_6))/(power_real_6+(1-FDR_real_6))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_62 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_62 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_62 <- 2*(power_real_62*(1-FDR_real_62))/(power_real_62+(1-FDR_real_62))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_6 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_6 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_6 <- 2*(power_mbi_6*(1-FDR_mbi_6))/(power_mbi_6+(1-FDR_mbi_6))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_6 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_6 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_6 <- 2*(power_sav_6*(1-FDR_sav_6))/(power_sav_6+(1-FDR_sav_6))
    
    
  }
  
  {
    Recall_A <- c(power_real_1,power_real_2,power_real_3,power_real_4,power_real_5,power_real_6)
    Recall_B <- c(power_real_12,power_real_22,power_real_32,power_real_42,power_real_52,power_real_62)
    Recall_C <- c(power_nb_1,power_nb_2,power_nb_3,power_nb_4,power_nb_5,power_nb_6)
    Recall_D <- c(power_poi_1,power_poi_2,power_poi_3,power_poi_4,power_poi_5,power_poi_6)
    Recall_E <- c(power_ppca_1,power_ppca_2,power_ppca_3,power_ppca_4,power_ppca_5,power_ppca_6)
    Recall_F <- c(power_nb_cov_1,power_nb_cov_2,power_nb_cov_3,power_nb_cov_4,power_nb_cov_5,power_nb_cov_6)
    Recall_G <- c(power_poi_cov_1,power_poi_cov_2,power_poi_cov_3,power_poi_cov_4,power_poi_cov_5,power_poi_cov_6)
    Recall_H <- c(power_ln_cov_1,power_ln_cov_2,power_ln_cov_3,power_ln_cov_4,power_ln_cov_5,power_ln_cov_6)
    Recall_I <- c(power_deseq_1,power_deseq_2,power_deseq_3,power_deseq_4,power_deseq_5,power_deseq_6)
    Recall_J <- c(power_edger_1,power_edger_2,power_edger_3,power_edger_4,power_edger_5,power_edger_6)
    Recall_K <- c(power_meta_1,power_meta_2,power_meta_3,power_meta_4,power_meta_5,power_meta_6)
    Recall_L <- c(power_meta_12,power_meta_22,power_meta_32,power_meta_42,power_meta_52,power_meta_62)
    Recall_M <- c(power_sav_1,power_sav_2,power_sav_3,power_sav_4,power_sav_5,power_sav_6)
    Recall_N <- c(power_mbi_1,power_mbi_2,power_mbi_3,power_mbi_4,power_mbi_5,power_mbi_6)
  }
  {
    Precision_A <- c(1-FDR_real_1,1-FDR_real_2,1-FDR_real_3,1-FDR_real_4,1-FDR_real_5,1-FDR_real_6)
    Precision_B <- c(1-FDR_real_12,1-FDR_real_22,1-FDR_real_32,1-FDR_real_42,1-FDR_real_52,1-FDR_real_62)
    Precision_C <- c(c(1-FDR_nb_1,1-FDR_nb_2,1-FDR_nb_3,1-FDR_nb_4,1-FDR_nb_5,1-FDR_nb_6))
    Precision_D <- c(1-FDR_poi_1,1-FDR_poi_2,1-FDR_poi_3,1-FDR_poi_4,1-FDR_poi_5,1-FDR_poi_6)
    Precision_E <- c(1-FDR_ppca_1,1-FDR_ppca_2,1-FDR_ppca_3,1-FDR_ppca_4,1-FDR_ppca_5,1-FDR_ppca_6)
    Precision_F <- c(c(1-FDR_nb_cov_1,1-FDR_nb_cov_2,1-FDR_nb_cov_3,1-FDR_nb_cov_4,1-FDR_nb_cov_5,1-FDR_nb_cov_6))
    Precision_G <- c(1-FDR_poi_cov_1,1-FDR_poi_cov_2,1-FDR_poi_cov_3,1-FDR_poi_cov_4,1-FDR_poi_cov_5,1-FDR_poi_cov_6)
    Precision_H <- c(1-FDR_ln_cov_1,1-FDR_ln_cov_2,1-FDR_ln_cov_3,1-FDR_ln_cov_4,1-FDR_ln_cov_5,1-FDR_ln_cov_6)
    Precision_I <- c(1-FDR_deseq_1,1-FDR_deseq_2,1-FDR_deseq_3,1-FDR_deseq_4,1-FDR_deseq_5,1-FDR_deseq_6)
    Precision_J <- c(1-FDR_edger_1,1-FDR_edger_2,1-FDR_edger_3,1-FDR_edger_4,1-FDR_edger_5,1-FDR_edger_6)
    Precision_K <- c(1-FDR_meta_1,1-FDR_meta_2,1-FDR_meta_3,1-FDR_meta_4,1-FDR_meta_5,1-FDR_meta_6)
    Precision_L <- c(1-FDR_meta_12,1-FDR_meta_22,1-FDR_meta_32,1-FDR_meta_42,1-FDR_meta_52,1-FDR_meta_62)
    Precision_M <- c(1-FDR_sav_1,1-FDR_sav_2,1-FDR_sav_3,1-FDR_sav_4,1-FDR_sav_5,1-FDR_sav_6)
    Precision_N <- c(1-FDR_mbi_1,1-FDR_mbi_2,1-FDR_mbi_3,1-FDR_mbi_4,1-FDR_mbi_5,1-FDR_mbi_6)
    
  }
  {
    F1_A <- c(score_real_1,score_real_2,score_real_3,score_real_4,score_real_5,score_real_6)
    F1_B <- c(score_real_12,score_real_22,score_real_32,score_real_42,score_real_52,score_real_62)
    F1_C <- c(score_nb_1,score_nb_2,score_nb_3,score_nb_4,score_nb_5,score_nb_6)
    F1_D <- c(score_poi_1,score_poi_2,score_poi_3,score_poi_4,score_poi_5,score_poi_6)
    F1_E <- c(score_ppca_1,score_ppca_2,score_ppca_3,score_ppca_4,score_ppca_5,score_ppca_6)
    F1_F <- c(score_nb_cov_1,score_nb_cov_2,score_nb_cov_3,score_nb_cov_4,score_nb_cov_5,score_nb_cov_6)
    F1_G <- c(score_poi_cov_1,score_poi_cov_2,score_poi_cov_3,score_poi_cov_4,score_poi_cov_5,score_poi_cov_6)
    F1_H <- c(score_ln_cov_1,score_ln_cov_2,score_ln_cov_3,score_ln_cov_4,score_ln_cov_5,score_ln_cov_6)
    F1_I <- c(score_deseq_1,score_deseq_2,score_deseq_3,score_deseq_4,score_deseq_5,score_deseq_6)
    F1_J <- c(score_edger_1,score_edger_2,score_edger_3,score_edger_4,score_edger_5,score_edger_6)
    F1_K <- c(score_meta_1,score_meta_2,score_meta_3,score_meta_4,score_meta_5,score_meta_6)
    F1_L <- c(score_meta_12,score_meta_22,score_meta_32,score_meta_42,score_meta_52,score_meta_62)
    F1_M <- c(score_sav_1,score_sav_2,score_sav_3,score_sav_4,score_sav_5,score_sav_6)
    F1_N <- c(score_mbi_1,score_mbi_2,score_mbi_3,score_mbi_4,score_mbi_5,score_mbi_6)
  }
  
  df <- cbind(melt(c(Recall_B,Recall_E,Recall_F,Recall_G,Recall_I,Recall_J,Recall_K,
                     Recall_M,Recall_N)),
              melt(c(Precision_B,Precision_E,Precision_F,Precision_G,Precision_I,Precision_J,Precision_K,
                     Precision_M,Precision_N)),
              melt(c(F1_B,F1_E,F1_F,F1_G,F1_I,F1_J,F1_K,
                     F1_M,F1_N)),
              melt(c(rep("No imputation_t",6),rep("PPCA-NB-cov",6),
                     rep("mbDenoise-ZINB-cov_t",6),rep("mbDenoise-ZIP-cov_t",6),rep("DESeq2",6),rep("edgeR",6),
                     rep("metagenomeSeq:fitFeaturemodel",6),rep("SAVER_t",6),rep("mbImpute_t",6))))
  colnames(df)[4] <- "method";colnames(df)[1] <- "Recall";colnames(df)[2] <- "Precision";colnames(df)[3] <- "F1 score";
  x <- rep(c(0:5),9*3)
  df2 <- cbind(x,melt(df))
  setting <- rep("2",sum(nrow(df2)))
  df2 <- cbind(df2,setting)
}
##M9
{
  load("~/mbDenoise/reproducibility/simulation/sim_lnm1.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_1 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_1 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_1 <- 2*(power_poi_1*(1-FDR_poi_1))/(power_poi_1+(1-FDR_poi_1))
    
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_1 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_1 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_1 <- 2*(power_ppca_1*(1-FDR_ppca_1))/(power_ppca_1+(1-FDR_ppca_1))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_1 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_1 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_1 <- 2*(power_nb_1*(1-FDR_nb_1))/(power_nb_1+(1-FDR_nb_1))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_1 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_1 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_1 <- 2*(power_ln_1*(1-FDR_ln_1))/(power_ln_1+(1-FDR_ln_1))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_1 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_1 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_1 <- 2*(power_poi_cov_1*(1-FDR_poi_cov_1))/(power_poi_cov_1+(1-FDR_poi_cov_1))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_1 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_1 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_1 <- 2*(power_nb_cov_1*(1-FDR_nb_cov_1))/(power_nb_cov_1+(1-FDR_nb_cov_1))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_1 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_1 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_1 <- 2*(power_ln_cov_1*(1-FDR_ln_cov_1))/(power_ln_cov_1+(1-FDR_ln_cov_1))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_1 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_1 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_1 <- 2*(power_edger_1*(1-FDR_edger_1))/(power_edger_1+(1-FDR_edger_1))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_1 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_1 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_1 <- 2*(power_deseq_1*(1-FDR_deseq_1))/(power_deseq_1+(1-FDR_deseq_1))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_1 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_1 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_1 <- 2*(power_meta_1*(1-FDR_meta_1))/(power_meta_1+(1-FDR_meta_1))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_12 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_12 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_12 <- 2*(power_meta_12*(1-FDR_meta_12))/(power_meta_12+(1-FDR_meta_12))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_1 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_1 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_1 <- 2*(power_real_1*(1-FDR_real_1))/(power_real_1+(1-FDR_real_1))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_12 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_12 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_12 <- 2*(power_real_12*(1-FDR_real_12))/(power_real_12+(1-FDR_real_12))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_1 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_1 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_1 <- 2*(power_mbi_1*(1-FDR_mbi_1))/(power_mbi_1+(1-FDR_mbi_1))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_1 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_1 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_1 <- 2*(power_sav_1*(1-FDR_sav_1))/(power_sav_1+(1-FDR_sav_1))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_lnm1_0.3.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_2 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_2 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_2 <- 2*(power_poi_2*(1-FDR_poi_2))/(power_poi_2+(1-FDR_poi_2))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_2 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_2 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_2 <- 2*(power_ppca_2*(1-FDR_ppca_2))/(power_ppca_2+(1-FDR_ppca_2))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_2 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_2 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_2 <- 2*(power_nb_2*(1-FDR_nb_2))/(power_nb_2+(1-FDR_nb_2))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_2 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_2 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_2 <- 2*(power_ln_2*(1-FDR_ln_2))/(power_ln_2+(1-FDR_ln_2))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_2 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_2 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_2 <- 2*(power_poi_cov_2*(1-FDR_poi_cov_2))/(power_poi_cov_2+(1-FDR_poi_cov_2))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_2 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_2 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_2 <- 2*(power_nb_cov_2*(1-FDR_nb_cov_2))/(power_nb_cov_2+(1-FDR_nb_cov_2))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_2 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_2 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_2 <- 2*(power_ln_cov_2*(1-FDR_ln_cov_2))/(power_ln_cov_2+(1-FDR_ln_cov_2))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_2 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_2 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_2 <- 2*(power_edger_2*(1-FDR_edger_2))/(power_edger_2+(1-FDR_edger_2))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_2 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_2 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_2 <- 2*(power_deseq_2*(1-FDR_deseq_2))/(power_deseq_2+(1-FDR_deseq_2))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_2 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_2 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_2 <- 2*(power_meta_2*(1-FDR_meta_2))/(power_meta_2+(1-FDR_meta_2))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_22 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_22 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_22 <- 2*(power_meta_22*(1-FDR_meta_22))/(power_meta_22+(1-FDR_meta_22))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_2 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_2 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_2 <- 2*(power_real_2*(1-FDR_real_2))/(power_real_2+(1-FDR_real_2))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_22 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_22 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_22 <- 2*(power_real_22*(1-FDR_real_22))/(power_real_22+(1-FDR_real_22))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_2 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_2 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_2 <- 2*(power_mbi_2*(1-FDR_mbi_2))/(power_mbi_2+(1-FDR_mbi_2))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_2 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_2 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_2 <- 2*(power_sav_2*(1-FDR_sav_2))/(power_sav_2+(1-FDR_sav_2))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_lnm1_0.6.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_3 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_3 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_3 <- 2*(power_poi_3*(1-FDR_poi_3))/(power_poi_3+(1-FDR_poi_3))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_3 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_3 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_3 <- 2*(power_ppca_3*(1-FDR_ppca_3))/(power_ppca_3+(1-FDR_ppca_3))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_3 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_3 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_3 <- 2*(power_nb_3*(1-FDR_nb_3))/(power_nb_3+(1-FDR_nb_3))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_3 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_3 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_3 <- 2*(power_ln_3*(1-FDR_ln_3))/(power_ln_3+(1-FDR_ln_3))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_3 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_3 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_3 <- 2*(power_poi_cov_3*(1-FDR_poi_cov_3))/(power_poi_cov_3+(1-FDR_poi_cov_3))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_3 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_3 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_3 <- 2*(power_nb_cov_3*(1-FDR_nb_cov_3))/(power_nb_cov_3+(1-FDR_nb_cov_3))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_3 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_3 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_3 <- 2*(power_ln_cov_3*(1-FDR_ln_cov_3))/(power_ln_cov_3+(1-FDR_ln_cov_3))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_3 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_3 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_3 <- 2*(power_edger_3*(1-FDR_edger_3))/(power_edger_3+(1-FDR_edger_3))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_3 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_3 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_3 <- 2*(power_deseq_3*(1-FDR_deseq_3))/(power_deseq_3+(1-FDR_deseq_3))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_3 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_3 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_3 <- 2*(power_meta_3*(1-FDR_meta_3))/(power_meta_3+(1-FDR_meta_3))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_32 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_32 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_32 <- 2*(power_meta_32*(1-FDR_meta_32))/(power_meta_32+(1-FDR_meta_32))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_3 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_3 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_3 <- 2*(power_real_3*(1-FDR_real_3))/(power_real_3+(1-FDR_real_3))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_32 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_32 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_32 <- 2*(power_real_32*(1-FDR_real_32))/(power_real_32+(1-FDR_real_32))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_3 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_3 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_3 <- 2*(power_mbi_3*(1-FDR_mbi_3))/(power_mbi_3+(1-FDR_mbi_3))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_3 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_3 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_3 <- 2*(power_sav_3*(1-FDR_sav_3))/(power_sav_3+(1-FDR_sav_3))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_lnm1_0.9.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_4 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_4 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_4 <- 2*(power_poi_4*(1-FDR_poi_4))/(power_poi_4+(1-FDR_poi_4))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_4 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_4 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_4 <- 2*(power_ppca_4*(1-FDR_ppca_4))/(power_ppca_4+(1-FDR_ppca_4))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_4 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_4 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_4 <- 2*(power_nb_4*(1-FDR_nb_4))/(power_nb_4+(1-FDR_nb_4))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_4 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_4 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_4 <- 2*(power_ln_4*(1-FDR_ln_4))/(power_ln_4+(1-FDR_ln_4))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_4 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_4 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_4 <- 2*(power_poi_cov_4*(1-FDR_poi_cov_4))/(power_poi_cov_4+(1-FDR_poi_cov_4))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_4 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_4 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_4 <- 2*(power_nb_cov_4*(1-FDR_nb_cov_4))/(power_nb_cov_4+(1-FDR_nb_cov_4))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_4 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_4 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_4 <- 2*(power_ln_cov_4*(1-FDR_ln_cov_4))/(power_ln_cov_4+(1-FDR_ln_cov_4))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_4 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_4 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_4 <- 2*(power_edger_4*(1-FDR_edger_4))/(power_edger_4+(1-FDR_edger_4))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_4 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_4 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_4 <- 2*(power_deseq_4*(1-FDR_deseq_4))/(power_deseq_4+(1-FDR_deseq_4))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_4 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_4 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_4 <- 2*(power_meta_4*(1-FDR_meta_4))/(power_meta_4+(1-FDR_meta_4))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_42 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_42 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_42 <- 2*(power_meta_42*(1-FDR_meta_42))/(power_meta_42+(1-FDR_meta_42))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_4 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_4 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_4 <- 2*(power_real_4*(1-FDR_real_4))/(power_real_4+(1-FDR_real_4))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_42 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_42 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_42 <- 2*(power_real_42*(1-FDR_real_42))/(power_real_42+(1-FDR_real_42))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_4 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_4 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_4 <- 2*(power_mbi_4*(1-FDR_mbi_4))/(power_mbi_4+(1-FDR_mbi_4))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_4 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_4 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_4 <- 2*(power_sav_4*(1-FDR_sav_4))/(power_sav_4+(1-FDR_sav_4))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_lnm1_1.2.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_5 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_5 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_5 <- 2*(power_poi_5*(1-FDR_poi_5))/(power_poi_5+(1-FDR_poi_5))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_5 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_5 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_5 <- 2*(power_ppca_5*(1-FDR_ppca_5))/(power_ppca_5+(1-FDR_ppca_5))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_5 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_5 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_5 <- 2*(power_nb_5*(1-FDR_nb_5))/(power_nb_5+(1-FDR_nb_5))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_5 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_5 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_5 <- 2*(power_ln_5*(1-FDR_ln_5))/(power_ln_5+(1-FDR_ln_5))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_5 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_5 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_5 <- 2*(power_poi_cov_5*(1-FDR_poi_cov_5))/(power_poi_cov_5+(1-FDR_poi_cov_5))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_5 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_5 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_5 <- 2*(power_nb_cov_5*(1-FDR_nb_cov_5))/(power_nb_cov_5+(1-FDR_nb_cov_5))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_5 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_5 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_5 <- 2*(power_ln_cov_5*(1-FDR_ln_cov_5))/(power_ln_cov_5+(1-FDR_ln_cov_5))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_5 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_5 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_5 <- 2*(power_edger_5*(1-FDR_edger_5))/(power_edger_5+(1-FDR_edger_5))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_5 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_5 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_5 <- 2*(power_deseq_5*(1-FDR_deseq_5))/(power_deseq_5+(1-FDR_deseq_5))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_5 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_5 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_5 <- 2*(power_meta_5*(1-FDR_meta_5))/(power_meta_5+(1-FDR_meta_5))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_52 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_52 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_52 <- 2*(power_meta_52*(1-FDR_meta_52))/(power_meta_52+(1-FDR_meta_52))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_5 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_5 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_5 <- 2*(power_real_5*(1-FDR_real_5))/(power_real_5+(1-FDR_real_5))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_52 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_52 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_52 <- 2*(power_real_52*(1-FDR_real_52))/(power_real_52+(1-FDR_real_52))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_5 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_5 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_5 <- 2*(power_mbi_5*(1-FDR_mbi_5))/(power_mbi_5+(1-FDR_mbi_5))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_5 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_5 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_5 <- 2*(power_sav_5*(1-FDR_sav_5))/(power_sav_5+(1-FDR_sav_5))
    
    
  }
  load("~/mbDenoise/reproducibility/simulation/sim_lnm1_1.5.RData")
  {
    All_P_poi <- apply(zip_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_6 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_6 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_6 <- 2*(power_poi_6*(1-FDR_poi_6))/(power_poi_6+(1-FDR_poi_6))
    
    All_P_ppca <- apply(ppcanb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ppca <- apply(ppcanb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ppca <- All_P_ppca-TP_ppca
    power_ppca_6 <- mean(na.omit(TP_ppca/(n.w*0.5)))
    FDR_ppca_6 <- mean(na.omit(ifelse(All_P_ppca==0,0,FP_ppca/All_P_ppca)))
    score_ppca_6 <- 2*(power_ppca_6*(1-FDR_ppca_6))/(power_ppca_6+(1-FDR_ppca_6))
    
    #nb
    All_P_nb <- apply(zinb_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_6 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_6 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_6 <- 2*(power_nb_6*(1-FDR_nb_6))/(power_nb_6+(1-FDR_nb_6))
    
    #ln
    All_P_ln <- apply(zilnm_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_6 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_6 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_6 <- 2*(power_ln_6*(1-FDR_ln_6))/(power_ln_6+(1-FDR_ln_6))
    
    All_P_poi <- apply(zip_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_poi <- apply(zip_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_poi <- All_P_poi-TP_poi
    power_poi_cov_6 <- mean(na.omit(TP_poi/(n.w*0.5)))
    FDR_poi_cov_6 <- mean(na.omit(ifelse(All_P_poi==0,0,FP_poi/All_P_poi)))
    score_poi_cov_6 <- 2*(power_poi_cov_6*(1-FDR_poi_cov_6))/(power_poi_cov_6+(1-FDR_poi_cov_6))
    
    #nb
    All_P_nb <- apply(zinb_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_nb <- apply(zinb_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_nb <- All_P_nb-TP_nb
    power_nb_cov_6 <- mean(na.omit(TP_nb/(n.w*0.5)))
    FDR_nb_cov_6 <- mean(na.omit(ifelse(All_P_nb==0,0,FP_nb/All_P_nb)))
    score_nb_cov_6 <- 2*(power_nb_cov_6*(1-FDR_nb_cov_6))/(power_nb_cov_6+(1-FDR_nb_cov_6))
    
    #ln
    All_P_ln <- apply(zilnm_cov_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_ln <- apply(zilnm_cov_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_ln <- All_P_ln-TP_ln
    power_ln_cov_6 <- mean(na.omit(TP_ln/(n.w*0.5)))
    FDR_ln_cov_6 <- mean(na.omit(ifelse(All_P_ln==0,0,FP_ln/All_P_ln)))
    score_ln_cov_6 <- 2*(power_ln_cov_6*(1-FDR_ln_cov_6))/(power_ln_cov_6+(1-FDR_ln_cov_6))
    
    #edger
    All_P_edger <- apply(edger_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_edger <- apply(edger_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_edger <- All_P_edger-TP_edger
    power_edger_6 <- mean(na.omit(TP_edger/(n.w*0.5)))
    FDR_edger_6 <- mean(na.omit(ifelse(All_P_edger==0,0,FP_edger/All_P_edger)))
    score_edger_6 <- 2*(power_edger_6*(1-FDR_edger_6))/(power_edger_6+(1-FDR_edger_6))
    
    #deseq
    All_P_deseq <- apply(deseq_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_deseq <- apply(deseq_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_deseq <- All_P_deseq-TP_deseq
    power_deseq_6 <- mean(na.omit(TP_deseq/(n.w*0.5)))
    FDR_deseq_6 <- mean(na.omit(ifelse(All_P_deseq==0,0,FP_deseq/All_P_deseq)))
    score_deseq_6 <- 2*(power_deseq_6*(1-FDR_deseq_6))/(power_deseq_6+(1-FDR_deseq_6))
    
    #metagenome
    All_P_meta<- apply(meta_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_6 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_6 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_6 <- 2*(power_meta_6*(1-FDR_meta_6))/(power_meta_6+(1-FDR_meta_6))
    
    All_P_meta<- apply(meta_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_meta<- apply(meta_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_meta<- All_P_meta-TP_meta
    power_meta_62 <- mean(na.omit(TP_meta/(n.w*0.5)))
    FDR_meta_62 <- mean(na.omit(ifelse(All_P_meta==0,0,FP_meta/All_P_meta)))
    score_meta_62 <- 2*(power_meta_62*(1-FDR_meta_62))/(power_meta_62+(1-FDR_meta_62))
    
    #real
    All_P_real <- apply(real_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_6 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_6 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_6 <- 2*(power_real_6*(1-FDR_real_6))/(power_real_6+(1-FDR_real_6))
    
    All_P_real <- apply(real_p_adj2[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_real <- apply(real_p_adj2[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_real <- All_P_real-TP_real
    power_real_62 <- mean(na.omit(TP_real/(n.w*0.5)))
    FDR_real_62 <- mean(na.omit(ifelse(All_P_real==0,0,FP_real/All_P_real)))
    score_real_62 <- 2*(power_real_62*(1-FDR_real_62))/(power_real_62+(1-FDR_real_62))
    
    All_P_mbi <- apply(mbimpute_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_mbi <- apply(mbimpute_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_mbi <- All_P_mbi-TP_mbi
    power_mbi_6 <- mean(na.omit(TP_mbi/(n.w*0.5)))
    FDR_mbi_6 <- mean(na.omit(ifelse(All_P_mbi==0,0,FP_mbi/All_P_mbi)))
    score_mbi_6 <- 2*(power_mbi_6*(1-FDR_mbi_6))/(power_mbi_6+(1-FDR_mbi_6))
    
    All_P_sav <- apply(saver_p_adj[,1:ref_p],2,function(x){sum(x<0.05)})
    TP_sav <- apply(saver_p_adj[(1:(n.w*0.5)),1:ref_p],2,function(x){sum(x<0.05)})
    FP_sav <- All_P_sav-TP_sav
    power_sav_6 <- mean(na.omit(TP_sav/(n.w*0.5)))
    FDR_sav_6 <- mean(na.omit(ifelse(All_P_sav==0,0,FP_sav/All_P_sav)))
    score_sav_6 <- 2*(power_sav_6*(1-FDR_sav_6))/(power_sav_6+(1-FDR_sav_6))
    
    
  }
  
  {
    Recall_A <- c(power_real_1,power_real_2,power_real_3,power_real_4,power_real_5,power_real_6)
    Recall_B <- c(power_real_12,power_real_22,power_real_32,power_real_42,power_real_52,power_real_62)
    Recall_C <- c(power_nb_1,power_nb_2,power_nb_3,power_nb_4,power_nb_5,power_nb_6)
    Recall_D <- c(power_poi_1,power_poi_2,power_poi_3,power_poi_4,power_poi_5,power_poi_6)
    Recall_E <- c(power_ppca_1,power_ppca_2,power_ppca_3,power_ppca_4,power_ppca_5,power_ppca_6)
    Recall_F <- c(power_nb_cov_1,power_nb_cov_2,power_nb_cov_3,power_nb_cov_4,power_nb_cov_5,power_nb_cov_6)
    Recall_G <- c(power_poi_cov_1,power_poi_cov_2,power_poi_cov_3,power_poi_cov_4,power_poi_cov_5,power_poi_cov_6)
    Recall_H <- c(power_ln_cov_1,power_ln_cov_2,power_ln_cov_3,power_ln_cov_4,power_ln_cov_5,power_ln_cov_6)
    Recall_I <- c(power_deseq_1,power_deseq_2,power_deseq_3,power_deseq_4,power_deseq_5,power_deseq_6)
    Recall_J <- c(power_edger_1,power_edger_2,power_edger_3,power_edger_4,power_edger_5,power_edger_6)
    Recall_K <- c(power_meta_1,power_meta_2,power_meta_3,power_meta_4,power_meta_5,power_meta_6)
    Recall_L <- c(power_meta_12,power_meta_22,power_meta_32,power_meta_42,power_meta_52,power_meta_62)
    Recall_M <- c(power_sav_1,power_sav_2,power_sav_3,power_sav_4,power_sav_5,power_sav_6)
    Recall_N <- c(power_mbi_1,power_mbi_2,power_mbi_3,power_mbi_4,power_mbi_5,power_mbi_6)
  }
  {
    Precision_A <- c(1-FDR_real_1,1-FDR_real_2,1-FDR_real_3,1-FDR_real_4,1-FDR_real_5,1-FDR_real_6)
    Precision_B <- c(1-FDR_real_12,1-FDR_real_22,1-FDR_real_32,1-FDR_real_42,1-FDR_real_52,1-FDR_real_62)
    Precision_C <- c(c(1-FDR_nb_1,1-FDR_nb_2,1-FDR_nb_3,1-FDR_nb_4,1-FDR_nb_5,1-FDR_nb_6))
    Precision_D <- c(1-FDR_poi_1,1-FDR_poi_2,1-FDR_poi_3,1-FDR_poi_4,1-FDR_poi_5,1-FDR_poi_6)
    Precision_E <- c(1-FDR_ppca_1,1-FDR_ppca_2,1-FDR_ppca_3,1-FDR_ppca_4,1-FDR_ppca_5,1-FDR_ppca_6)
    Precision_F <- c(c(1-FDR_nb_cov_1,1-FDR_nb_cov_2,1-FDR_nb_cov_3,1-FDR_nb_cov_4,1-FDR_nb_cov_5,1-FDR_nb_cov_6))
    Precision_G <- c(1-FDR_poi_cov_1,1-FDR_poi_cov_2,1-FDR_poi_cov_3,1-FDR_poi_cov_4,1-FDR_poi_cov_5,1-FDR_poi_cov_6)
    Precision_H <- c(1-FDR_ln_cov_1,1-FDR_ln_cov_2,1-FDR_ln_cov_3,1-FDR_ln_cov_4,1-FDR_ln_cov_5,1-FDR_ln_cov_6)
    Precision_I <- c(1-FDR_deseq_1,1-FDR_deseq_2,1-FDR_deseq_3,1-FDR_deseq_4,1-FDR_deseq_5,1-FDR_deseq_6)
    Precision_J <- c(1-FDR_edger_1,1-FDR_edger_2,1-FDR_edger_3,1-FDR_edger_4,1-FDR_edger_5,1-FDR_edger_6)
    Precision_K <- c(1-FDR_meta_1,1-FDR_meta_2,1-FDR_meta_3,1-FDR_meta_4,1-FDR_meta_5,1-FDR_meta_6)
    Precision_L <- c(1-FDR_meta_12,1-FDR_meta_22,1-FDR_meta_32,1-FDR_meta_42,1-FDR_meta_52,1-FDR_meta_62)
    Precision_M <- c(1-FDR_sav_1,1-FDR_sav_2,1-FDR_sav_3,1-FDR_sav_4,1-FDR_sav_5,1-FDR_sav_6)
    Precision_N <- c(1-FDR_mbi_1,1-FDR_mbi_2,1-FDR_mbi_3,1-FDR_mbi_4,1-FDR_mbi_5,1-FDR_mbi_6)
    
  }
  {
    F1_A <- c(score_real_1,score_real_2,score_real_3,score_real_4,score_real_5,score_real_6)
    F1_B <- c(score_real_12,score_real_22,score_real_32,score_real_42,score_real_52,score_real_62)
    F1_C <- c(score_nb_1,score_nb_2,score_nb_3,score_nb_4,score_nb_5,score_nb_6)
    F1_D <- c(score_poi_1,score_poi_2,score_poi_3,score_poi_4,score_poi_5,score_poi_6)
    F1_E <- c(score_ppca_1,score_ppca_2,score_ppca_3,score_ppca_4,score_ppca_5,score_ppca_6)
    F1_F <- c(score_nb_cov_1,score_nb_cov_2,score_nb_cov_3,score_nb_cov_4,score_nb_cov_5,score_nb_cov_6)
    F1_G <- c(score_poi_cov_1,score_poi_cov_2,score_poi_cov_3,score_poi_cov_4,score_poi_cov_5,score_poi_cov_6)
    F1_H <- c(score_ln_cov_1,score_ln_cov_2,score_ln_cov_3,score_ln_cov_4,score_ln_cov_5,score_ln_cov_6)
    F1_I <- c(score_deseq_1,score_deseq_2,score_deseq_3,score_deseq_4,score_deseq_5,score_deseq_6)
    F1_J <- c(score_edger_1,score_edger_2,score_edger_3,score_edger_4,score_edger_5,score_edger_6)
    F1_K <- c(score_meta_1,score_meta_2,score_meta_3,score_meta_4,score_meta_5,score_meta_6)
    F1_L <- c(score_meta_12,score_meta_22,score_meta_32,score_meta_42,score_meta_52,score_meta_62)
    F1_M <- c(score_sav_1,score_sav_2,score_sav_3,score_sav_4,score_sav_5,score_sav_6)
    F1_N <- c(score_mbi_1,score_mbi_2,score_mbi_3,score_mbi_4,score_mbi_5,score_mbi_6)
  }
  
  df <- cbind(melt(c(Recall_B,Recall_E,Recall_F,Recall_G,Recall_I,Recall_J,Recall_K,
                     Recall_M,Recall_N)),
              melt(c(Precision_B,Precision_E,Precision_F,Precision_G,Precision_I,Precision_J,Precision_K,
                     Precision_M,Precision_N)),
              melt(c(F1_B,F1_E,F1_F,F1_G,F1_I,F1_J,F1_K,
                     F1_M,F1_N)),
              melt(c(rep("No imputation_t",6),rep("PPCA-NB-cov",6),
                     rep("mbDenoise-ZINB-cov_t",6),rep("mbDenoise-ZIP-cov_t",6),rep("DESeq2",6),rep("edgeR",6),
                     rep("metagenomeSeq:fitFeaturemodel",6),rep("SAVER_t",6),rep("mbImpute_t",6))))
  colnames(df)[4] <- "method";colnames(df)[1] <- "Recall";colnames(df)[2] <- "Precision";colnames(df)[3] <- "F1 score";
  x <- rep(c(0:5),9*3)
  #x <- rep(c(0,0.3,0.6,0.9,1.2,1.5),9*3)
  df3 <- cbind(x,melt(df))
  setting <- rep("3",sum(nrow(df3)))
  df3 <- cbind(df3,setting)
}

df_all <- rbind(df1,df2,df3)
df_all$method <- factor(df_all$method,labels = c("DESeq2" ,"edgeR" ,"mbDenoise-zinb-cov",
                                                 "mbDenoise-zip-cov", "mbImpute",                   
                                                 "metagenomeSeq","t test","PPCA-NB-cov","SAVER"))
df_all$setting <- factor(df_all$setting,labels = c("M7","M8","M9"))
colourCount = length(unique(df_all$method))

p <- ggplot(data=df_all, aes(x=x, y=value, shape=factor(method,
                                                        level=c("t test","DESeq2","edgeR","metagenomeSeq","PPCA-NB-cov","mbImpute","SAVER","mbDenoise-zip-cov",
                                                                "mbDenoise-zinb-cov")),
                             colour=factor(method,level=c("t test","DESeq2","edgeR","metagenomeSeq","PPCA-NB-cov",
                                                          "mbImpute","SAVER","mbDenoise-zip-cov",
                                                          "mbDenoise-zinb-cov")),
                             linetype=factor(method,level=c("t test","DESeq2","edgeR","metagenomeSeq","PPCA-NB-cov",
                                                            "mbImpute","SAVER","mbDenoise-zip-cov",
                                                            "mbDenoise-zinb-cov")))) + 
  geom_line(size=0.65) + geom_point(size=2)+
  scale_shape_manual(name  ="Method",labels=c(expression(paste(italic(t), "-test")),"DESeq2","edgeR","metagenomeSeq","PPCA-NB-cov",
                                              "mbImpute","SAVER","mbDenoise-zip-cov",
                                              "mbDenoise-zinb-cov"),values=c(0:8))+
  scale_linetype_manual(name  ="Method",labels=c(expression(paste(italic(t), "-test")),"DESeq2","edgeR","metagenomeSeq","PPCA-NB-cov",
                                                 "mbImpute","SAVER","mbDenoise-zip-cov",
                                                 "mbDenoise-zinb-cov"),
                        values=c('dashed', 'dashed', 'dashed',"dashed","dashed",
                                 'dashed', 'dashed', 'dashed',"dashed","dashed"))+
  scale_color_manual(name  ="Method",labels=c(expression(paste(italic(t), "-test")),"DESeq2","edgeR","metagenomeSeq","PPCA-NB-cov",
                                              "mbImpute","SAVER","mbDenoise-zip-cov",
                                              "mbDenoise-zinb-cov"),
                     values = c(colorRampPalette(brewer.pal(9, "Paired"))(colourCount)[-6],colorRampPalette(brewer.pal(9, "Paired"))(colourCount)[6]))+
  theme_bw()+
  theme(strip.text =element_text(size = 13),axis.title.x =element_text(size=14)) +
  labs(x = "Effect size",y = "", title = "")+ 
  guides(fill = guide_legend(title = "Method"))+theme(legend.text.align = 0)+
  facet_grid(variable ~ setting,scales="free_x")

p


