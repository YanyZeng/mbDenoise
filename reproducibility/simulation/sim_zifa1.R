setwd("~/mbDenoise/reproducibility")

n.n = 50
n.w = 150
n.factors = 2
iter = 100
library(vegan)
library(MASS)
library(reticulate)
library(transport)
library(edgeR)
library(DESeq2)
library(mbImpute)
library(SAVER)
library(stats)
library(coin)
library(metagenomeSeq)
library(glmnet)
library(DirichletMultinomial)
library('Matrix')
library('parallel')
library('doParallel')
library('foreach')
library('optimx')
library('trustOptim')
library('ZIPFA')
library(mbDenoise)

source("other_software/proximal_gradient.R")
source("other_software/basic_functions.R")
source("other_software/tune_proximal_gradient.R")
source("other_software/naive_methods.R")
source_python('other_software/ZIFA_copy.py')

{
  zero.prob <- rep(0,iter)
  zeroinfl.prob <- rep(0,iter)
  depth <- rep(0,iter)

  m_pr_y <- rep(0,iter)
  beta_pr_y <- rep(0,iter)
  beta_pr_y2 <- rep(0,iter)
  m_pr_y2 <- rep(0,iter)
  difftime_pr_y <- rep(0,iter)
  
  m_pcoa_y <- rep(0,iter)
  m_pcoa_y2 <- rep(0,iter)
  difftime_pcoa_y <- rep(0,iter)

  m_tsne_y <- rep(0,iter)
  m_tsne_y2 <- rep(0,iter)
  difftime_tsne_y <- rep(0,iter)
  
  beta_zifa2 <- rep(0,iter)
  beta_zifa <- rep(0,iter)
  m_zifa <- rep(0,iter)
  m_zifa2 <- rep(0,iter)
  difftime_zifa <- rep(0,iter)
  
  m_nbva <- rep(0,iter)
  beta_nbva2 <- rep(0,iter)
  beta_nbva <- rep(0,iter)
  m_nbva2 <- rep(0,iter)
  difftime_nbva <- rep(0,iter)
  
  m_nbva_cov <- rep(0,iter)
  beta_nbva_cov2 <- rep(0,iter)
  beta_nbva_cov <- rep(0,iter)
  m_nbva_cov2 <- rep(0,iter)
  gamma_nbva_cov <- rep(0,iter)
  gamma_nbva_cov2 <- rep(0,iter)
  difftime_nbva_cov <- rep(0,iter)
  
  m_zip <- rep(0,iter)
  beta_zip2 <- rep(0,iter)
  beta_zip <- rep(0,iter)
  m_zip2 <- rep(0,iter)
  difftime_zip <- rep(0,iter)
  
  m_zinb <- rep(0,iter)
  beta_zinb2 <- rep(0,iter)
  beta_zinb <- rep(0,iter)
  m_zinb2 <- rep(0,iter)
  difftime_zinb <- rep(0,iter)
  
  m_zip_cov <- rep(0,iter)
  beta_zip_cov2 <- rep(0,iter)
  beta_zip_cov <- rep(0,iter)
  m_zip_cov2 <- rep(0,iter)
  difftime_zip_cov <- rep(0,iter)
  gamma_zip_cov <- rep(0,iter)
  gamma_zip_cov2 <- rep(0,iter)
  
  m_zinb_cov <- rep(0,iter)
  beta_zinb_cov2 <- rep(0,iter)
  beta_zinb_cov <- rep(0,iter)
  m_zinb_cov2 <- rep(0,iter)
  difftime_zinb_cov <- rep(0,iter)
  gamma_zinb_cov <- rep(0,iter)
  gamma_zinb_cov2 <- rep(0,iter)
  
  difftime_pg <- rep(0,iter)
  difftime_zr <- rep(0,iter)
  difftime_svt <- rep(0,iter)
  
  a_0.5 <- rep(0,iter)
  a_nb <- rep(0,iter)
  a_poi <- rep(0,iter)
  a_svt <- rep(0,iter)
  a_multi <- rep(0,iter)
  a_nb_cov <- rep(0,iter)
  a_poi_cov <- rep(0,iter)
 
  b_0.5 <- rep(0,iter)
  b_nb <- rep(0,iter)
  b_poi <- rep(0,iter)
  b_svt <- rep(0,iter)
  b_multi <- rep(0,iter)
  b_nb_cov <- rep(0,iter)
  b_poi_cov <- rep(0,iter)

  c_0.5 <- rep(0,iter)
  c_nb <- rep(0,iter)
  c_poi <- rep(0,iter)
  c_svt <- rep(0,iter)
  c_multi <- rep(0,iter)
  c_nb_cov <- rep(0,iter)
  c_poi_cov <- rep(0,iter)

  d_0.5 <- rep(0,iter)
  d_nb <- rep(0,iter)
  d_poi <- rep(0,iter)
  d_svt <- rep(0,iter)
  d_multi <- rep(0,iter)
  d_nb_cov <- rep(0,iter)
  d_poi_cov <- rep(0,iter)

  a_dmn <- rep(0,iter)
  b_dmn <- rep(0,iter)
  c_dmn <- rep(0,iter)
  d_dmn <- rep(0,iter)
  
  a_mbi <- rep(0,iter)
  a_sav <- rep(0,iter)
  b_mbi <- rep(0,iter)
  b_sav <- rep(0,iter)
  c_mbi <- rep(0,iter)
  c_sav <- rep(0,iter)
  d_mbi <- rep(0,iter)
  d_sav <- rep(0,iter)

  a_ppcanb <- rep(0,iter)
  a_ppcanb_cov <- rep(0,iter)
  b_ppcanb <- rep(0,iter)
  b_ppcanb_cov <- rep(0,iter)
  c_ppcanb <- rep(0,iter)
  c_ppcanb_cov <- rep(0,iter)
  d_ppcanb <- rep(0,iter)
  d_ppcanb_cov <- rep(0,iter)
  
  ppcanb_cor <-  ppcanb_cov_cor <- vector("list", iter)
  ppcanb_mse2 <-  ppcanb_cov_mse2 <- vector("list", iter)
  ppcanb_wasserstein1d  <-  ppcanb_cov_wasserstein1d <- vector("list", iter)
  ppcanb_mean <- ppcanb_cov_mean <- vector("list", iter)
  ppcanb_sd <- ppcanb_cov_sd <- vector("list", iter)
  ppcanb_p_adj <-ppcanb_cov_p_adj <-matrix(0,n.w,iter)
  ppcanb_p <-ppcanb_cov_p <-matrix(0,n.w,iter)
  
  noimpute_cor <- zinb_cor <- zip_cor <- sav_cor <- mbi_cor <- vector("list", iter)
  zinb_cov_cor <- zip_cov_cor <- vector("list", iter)

  noimpute_wasserstein1d <- zinb_wasserstein1d <- zip_wasserstein1d <- mbi_wasserstein1d <- sav_wasserstein1d <- rep(0,iter)
  zinb_cov_wasserstein1d <- zip_cov_wasserstein1d <-  rep(0,iter)
  complete_mean <- complete_sd <- noimpute_mean <- noimpute_sd <- vector("list", iter)
  zinb_mean <- zinb_sd <- zip_mean<- zip_sd <-  vector("list", iter)
  mbi_mean <- mbi_sd <- sav_mean <- sav_sd <-  vector("list", iter)
  zinb_cov_mean <- zinb_cov_sd <- zip_cov_mean<- zip_cov_sd <-  vector("list", iter)
  
  noimpute_mse2 <- rep(0,iter)
  zinb_mse2<- rep(0,iter)
  zip_mse2<- rep(0,iter)
  zinb_cov_mse2<- rep(0,iter)
  zip_cov_mse2<- rep(0,iter)
  sav_mse2<- rep(0,iter)
  mbi_mse2<- rep(0,iter)

  zip_p_adj <-zip_cov_p_adj <-matrix(0,n.w,iter)
  zinb_p_adj <- zinb_cov_p_adj <-matrix(0,n.w,iter)
  real_p_adj <-real_p_adj2 <- matrix(0,n.w,iter)
  zip_p <- zinb_p <- real_p <- real_p2 <- matrix(0,n.w,iter)
  zip_cov_p <- zinb_cov_p <-   matrix(0,n.w,iter)
  edger_p_adj <-  matrix(0,n.w,iter);
  deseq_p_adj <- matrix(0,n.w,iter)
  meta_p_adj <- matrix(0,n.w,iter);meta_p_adj2 <- matrix(0,n.w,iter);
  meta_p <- matrix(0,n.w,iter);meta_p2 <- matrix(0,n.w,iter);
  mbimpute_p_adj <- matrix(0,n.w,iter)
  saver_p_adj <-  matrix(0,n.w,iter)
  mbimpute_p <-  matrix(0,n.w,iter)
  saver_p <- matrix(0,n.w,iter)
  difftime_mbi <- difftime_sav <- rep(0,iter)
  difftime_dmn <- rep(0,iter)
 
}

for(k in 1:iter){
  print(k)
  set.seed(k)
  
  f <- matrix(0,nrow = n.n, ncol = n.factors)
  for(i in 1:n.n){
    f[i,] <- rnorm(n.factors, mean = 0, sd = 1)
  }
  betaj <- matrix(0,nrow = n.w, ncol = n.factors)
  for(j in 1:n.w){
    betaj[j,] <- runif(n.factors, -0.5, 0.5)
    
  }
  sigma <- diag(runif(n.w,0.9,1.1)*0.3,n.w,n.w)
  beta0 <- runif(n.w,2.7,3.3)
  
  u <- matrix(beta0,n.n,n.w,byrow=TRUE) + f %*% t(betaj)
  X <- matrix(0,n.n,n.w,byrow = TRUE)
  for(i in 1:n.n){
    X[i,] <- mvrnorm(n=1, mu =u[i,], Sigma = sigma)
  }
  eta <- exp(-0.1*(X^2))
  X <- round(exp(X),0)
  
  z <- matrix(0,n.n,n.w)
  for(i in 1:n.n){
    z[i,] <- rbinom(n.w, size=1, prob=eta[i,])
  }
  X[z==1]=0
  
  X_cov<- c(rep(1,n.n/2),rep(0,n.n/2))
  X_ori <- (1-z)*exp(u)
  gamma <- rep(0,n.w)
  
  Y <- log2(1+X)
  zerorow <- which(rowSums(X)==0)
  if(length(zerorow) >0 ){
    X <- X[-zerorow,];X_ori <- X_ori[-zerorow,];Y <- Y[-zerorow,];X_cov<-X_cov[-zerorow];
    f <- f[-zerorow,];
  }
  zerocol <- which(colSums(X)==0)
  if(length(zerocol) >0 ){
    X <- X[,-zerocol];X_ori <- X_ori[,-zerocol];Y <- Y[,-zerocol];
    betaj <- t(t(betaj)[,-zerocol])};
  
  n.nn <- nrow(X)
  n.mm <- ncol(X)
  
  depth[k] <- summary(rowSums(X))[6]/summary(rowSums(X))[1]
  zero.prob[k]  <- sum(X==0)/(dim(X)[1]*dim(X)[2])
  
  ##mbDenoise-zip-cov
  t7 <- Sys.time()
  re_zip_cov <- tryCatch({ZIPPCApn(X,X_cov,family="poisson")},error=function(e){NULL})
  t8 <- Sys.time()
  difftime_zip_cov[k] <- difftime(t8,t7,units = "secs")
  if(is.null(re_zip_cov$params$factor_coefs_j) || is.null(re_zip_cov$lvs$factor_scores )){
    m_zip_cov[k] <- NaN;m_zip_cov2[k] <- NaN;beta_zip_cov[k] <- NaN;beta_zip_cov2[k] <- NaN
    gamma_zip_cov[k] <- NaN;gamma_zip_cov2[k] <- NaN
  }else{
    qa <-qr.Q(qr(betaj))
    qb <- qr.Q(qr(re_zip_cov$params$factor_coefs_j))
    qa2 <-qr.Q(qr(f))
    qb2 <- qr.Q(qr(re_zip_cov$lvs$factor_scores ))
    qc <-qr.Q(qr(gamma))
    qc2 <- qr.Q(qr(re_zip_cov$params$gamma))
    m_zip_cov[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2)
    beta_zip_cov[k] <- sum((qa %*% t(qa) - qb %*% t(qb))^2 )
    gamma_zip_cov[k] <- sum((qc %*% t(qc) - qc2 %*% t(qc2))^2 )
    m_zip_cov2[k] <- (procrustes(f,re_zip_cov$lvs$factor_scores,symmetric = TRUE)$ss)
    beta_zip_cov2[k] <- (procrustes(betaj,re_zip_cov$params$factor_coefs_j,symmetric = TRUE)$ss)
    gamma_zip_cov2[k] <- tryCatch({(procrustes(gamma,re_zip_cov$params$gamma,symmetric = TRUE)$ss)},error=function(e){NaN})
  }
  
  ##mbDenoise-zinb-cov
  t7 <- Sys.time()
  re_zinb_cov <- tryCatch({ZIPPCApn(X,X_cov)},error=function(e){NULL})
  t8 <- Sys.time()
  difftime_zinb_cov[k] <- difftime(t8,t7,units = "secs")
  if(is.null(re_zinb_cov$params$factor_coefs_j) || is.null(re_zinb_cov$lvs$factor_scores )){
    m_zinb_cov[k] <- NaN;m_zinb_cov2[k] <- NaN;beta_zinb_cov[k] <- NaN;beta_zinb_cov2[k] <- NaN
    gamma_zinb_cov[k] <- NaN;gamma_zinb_cov2[k] <- NaN
  }else{
    qa <-qr.Q(qr(betaj))
    qb <- qr.Q(qr(re_zinb_cov$params$factor_coefs_j))
    qa2 <-qr.Q(qr(f))
    qb2 <- qr.Q(qr(re_zinb_cov$lvs$factor_scores))
    qc <-qr.Q(qr(gamma))
    qc2 <- qr.Q(qr(re_zinb_cov$params$gamma))
    m_zinb_cov[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2)
    beta_zinb_cov[k] <- sum((qa %*% t(qa) - qb %*% t(qb))^2 )
    gamma_zinb_cov[k] <- sum((qc %*% t(qc) - qc2 %*% t(qc2))^2 )
    m_zinb_cov2[k] <- (procrustes(f,re_zinb_cov$lvs$factor_scores,symmetric = TRUE)$ss)
    beta_zinb_cov2[k] <- (procrustes(betaj,re_zinb_cov$params$factor_coefs_j,symmetric = TRUE)$ss)
    gamma_zinb_cov2[k] <- tryCatch({(procrustes(gamma,re_zinb_cov$params$gamma,symmetric = TRUE)$ss)},error=function(e){NaN})
  }
  
  ##PPCA-NB-cov
  t1 <- Sys.time()
  re_nbva_cov <- tryCatch({gllvm::gllvm(X,as.data.frame(X_cov),family = "negative.binomial",method = "VA",Lambda.struc = "diagonal",row.eff="fixed")},
                          error=function(e){NULL})
  t2 <- Sys.time()
  difftime_nbva_cov[k] <- difftime(t2,t1,units = "secs")
  if(is.null(re_nbva_cov$params$theta) || is.null(re_nbva_cov$lvs)){
    beta_nbva[k] <- NaN;m_nbva[k] <- NaN;beta_nbva2[k] <- NaN;m_nbva2[k] <- NaN
  }else{
    qa <-qr.Q(qr(betaj))
    qb <- qr.Q(qr(re_nbva_cov$params$theta))
    qa2 <-qr.Q(qr(f))
    qb2 <- qr.Q(qr(re_nbva_cov$lvs))
    qc <-qr.Q(qr(gamma))
    qc2 <- qr.Q(qr(re_nbva_cov$params$Xcoef))
    m_nbva_cov[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2)
    beta_nbva_cov[k] <- sum((qa %*% t(qa) - qb %*% t(qb))^2 )
    gamma_nbva_cov[k] <- sum((qc %*% t(qc) - qc2 %*% t(qc2))^2 )
    m_nbva_cov2[k] <- (procrustes(f,re_nbva_cov$lvs,symmetric = TRUE)$ss)
    beta_nbva_cov2[k] <- (procrustes(betaj,re_nbva_cov$params$theta,symmetric = TRUE)$ss)
    gamma_nbva_cov2[k] <- tryCatch({(procrustes(gamma,re_nbva_cov$params$Xcoef,symmetric = TRUE)$ss)},error=function(e){NaN})
  }
  
  ##PPCA-NB
  t1 <- Sys.time()
  re_nbva <- tryCatch({gllvm::gllvm(X,family = "negative.binomial",method = "VA",Lambda.struc = "diagonal",row.eff= "fixed")},
                      error=function(e){NULL})
  t2 <- Sys.time()
  difftime_nbva[k] <- difftime(t2,t1,units = "secs")
  if(is.null(re_nbva$params$theta) || is.null(re_nbva$lvs)){
    beta_nbva[k] <- NaN;m_nbva[k] <- NaN;beta_nbva2[k] <- NaN;m_nbva2[k] <- NaN
  }else{
    qa <-qr.Q(qr(betaj))
    qb <- qr.Q(qr(re_nbva$params$theta))
    qa2 <-qr.Q(qr(f))
    qb2 <- qr.Q(qr(re_nbva$lvs))
    m_nbva[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2)
    beta_nbva[k] <- sum((qa %*% t(qa) - qb %*% t(qb))^2 )
    m_nbva2[k] <- (procrustes(f,re_nbva$lvs,symmetric = TRUE)$ss)
    beta_nbva2[k] <- (procrustes(betaj,re_nbva$params$theta,symmetric = TRUE)$ss)
  }
  
  ##ZIFA
  t3 <- Sys.time()
  re_zifa <- tryCatch({fitModel(Y,n.factors)},
                      error=function(e){NULL})
  t4 <- Sys.time()
  difftime_zifa[k] <- difftime(t4,t3,units = "secs")
  if(is.null(re_zifa)){
    m_zifa2[k] <- NaN;beta_zifa2[k] <- NaN;m_zifa[k] <- NaN;beta_zifa[k] <- NaN
  }else{
    qa <-qr.Q(qr(betaj))
    qb <- qr.Q(qr(re_zifa[2][[1]]$A))
    qa2 <-qr.Q(qr(f))
    qb2 <- qr.Q(qr(re_zifa[[1]]))
    m_zifa[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2 )
    beta_zifa[k] <- sum((qa %*% t(qa) - qb %*% t(qb))^2 )
    m_zifa2[k] <- (procrustes(f,re_zifa[[1]],symmetric = TRUE)$ss)
    beta_zifa2[k] <- (procrustes(betaj,re_zifa[2][[1]]$A,symmetric = TRUE)$ss)
  }
  
  ##t-SNE
  t7 <- Sys.time()
  tsne_out <- tryCatch({Rtsne::Rtsne(Y,perplexity=round((n.nn-2)/3), initial_dims = n.mm)},
                       error=function(e){ NULL})
  t8 <- Sys.time()
  difftime_tsne_y[k] <- difftime(t8,t7,units = "secs")
  if(is.null(tsne_out)){
    m_tsne_y[k] <- NaN;m_tsne_y2[k] <- NaN;
  }else{
    qa2 <-qr.Q(qr(f))
    qb2 <- qr.Q(qr(tsne_out$Y))
    m_tsne_y[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2 )
    m_tsne_y2[k] <- (procrustes(f,tsne_out$Y,symmetric = TRUE)$ss)
  }
  
  ##PCA
  t7 <- Sys.time()
  pr <- prcomp(Y)
  t8 <- Sys.time()
  difftime_pr_y[k] <- difftime(t8,t7,units = "secs")
  qa <-qr.Q(qr(betaj))
  qb <- qr.Q(qr(pr$rotation[,1:2]))
  qa2 <-qr.Q(qr(f))
  qb2 <- qr.Q(qr(pr$x[,1:2]))
  m_pr_y[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2 )
  beta_pr_y[k] <- sum((qa %*% t(qa) - qb %*% t(qb))^2 )
  m_pr_y2[k] <- (procrustes(f,pr$x[,1:2],symmetric = TRUE)$ss)
  beta_pr_y2[k] <- (procrustes(betaj,pr$rotation[,1:2],symmetric = TRUE)$ss)
  
  #### PCoA (Bray-Curtis)
  t7 <- Sys.time()
  dist_bray <- vegan::vegdist(Y, method="bray")
  pcoa <- ape::pcoa(dist_bray)
  t8 <- Sys.time()
  difftime_pcoa_y[k] <- difftime(t8,t7,units = "secs")
  qa2 <-qr.Q(qr(f))
  qb2 <- qr.Q(qr(pcoa$vectors[,1:2]))
  m_pcoa_y[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2 )
  m_pcoa_y2[k] <- (procrustes(f,pcoa$vectors[,1:2],symmetric = TRUE)$ss)

  ##mbDenoise-zip
  t7 <- Sys.time()
  re_zip <- tryCatch({ZIPPCApn(X,family="poisson")},error=function(e){NULL})
  t8 <- Sys.time()
  difftime_zip[k] <- difftime(t8,t7,units = "secs")
  if(is.null(re_zip$params$factor_coefs_j) || is.null(re_zip$lvs$factor_scores )){
    m_zip[k] <- NaN;m_zip2[k] <- NaN;beta_zip[k] <- NaN;beta_zip2[k] <- NaN
  }else{
    qa <-qr.Q(qr(betaj))
    qb <- qr.Q(qr(re_zip$params$factor_coefs_j))
    qa2 <-qr.Q(qr(f))
    qb2 <- qr.Q(qr(re_zip$lvs$factor_scores ))
    beta_zip[k] <- sum((qa %*% t(qa) - qb %*% t(qb))^2 )
    m_zip[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2 )
    m_zip2[k] <- (procrustes(f,re_zip$lvs$factor_scores,symmetric = TRUE)$ss)
    beta_zip2[k] <- (procrustes(betaj,re_zip$params$factor_coefs_j,symmetric = TRUE)$ss)
  }
  
  ##mbDenoise-zinb
  t7 <- Sys.time()
  re_zinb <- tryCatch({ZIPPCApn(X)},error=function(e){NULL})
  t8 <- Sys.time()
  difftime_zinb[k] <- difftime(t8,t7,units = "secs")
  if(is.null(re_zinb$params$factor_coefs_j) || is.null(re_zinb$lvs$factor_scores )){
    m_zinb[k] <- NaN;m_zinb2[k] <- NaN;beta_zinb[k] <- NaN;beta_zinb2[k] <- NaN
  }else{
    qa <-qr.Q(qr(betaj))
    qb <- qr.Q(qr(re_zinb$params$factor_coefs_j))
    qa2 <-qr.Q(qr(f))
    qb2 <- qr.Q(qr(re_zinb$lvs$factor_scores))
    m_zinb[k] <- sum((qa2 %*% t(qa2) - qb2 %*% t(qb2))^2 )
    beta_zinb[k] <- sum((qa %*% t(qa) - qb %*% t(qb))^2 )
    m_zinb2[k] <- (procrustes(f,re_zinb$lvs$factor_scores,symmetric = TRUE)$ss)
    beta_zinb2[k] <- (procrustes(betaj,re_zinb$params$factor_coefs_j,symmetric = TRUE)$ss)
  }
  
  ##mbImpute
  t7 <- Sys.time()
  mbi <- tryCatch({mbImpute( condition = X_cov,otu_tab = X,parallel = TRUE, ncores = 4)},
                  error=function(e){ NULL})
  t8 <- Sys.time()
  difftime_mbi[k] <- difftime(t8,t7,units = "secs")
  
  ##SAVER
  t7 <- Sys.time()
  sav <- tryCatch({saver(t(X), ncores = 12)},
                  error=function(e){ NULL})
  t8 <- Sys.time()
  difftime_sav[k] <- difftime(t8,t7,units = "secs")
  

  ppcanb_cov <- tryCatch({exp(re_nbva_cov$lvs %*% t(re_nbva_cov$params$theta) + matrix(re_nbva_cov$params$Xcoef,n.nn,n.mm,byrow=TRUE)*X_cov+matrix(re_nbva_cov$params$beta0,n.nn,n.mm,byrow=TRUE))},
                         error=function(e){NaN})
  ppcanb <- tryCatch({exp(re_nbva$lvs %*% t(re_nbva$params$theta) + matrix(re_nbva$params$beta0,n.nn,n.mm,byrow=TRUE))},
                     error=function(e){NaN})
  ##mse
    noimpute_mse2[k] <-tryCatch({(sum((log2(X_ori+1)-log2(X+1))^2,na.rm=T)/(n.nn*n.mm))},
                                error=function(e){NaN})
    zinb_mse2[k] <-tryCatch({(sum((log2(X_ori+1)-log2(re_zinb$muz+1))^2,na.rm=T)/(n.nn*n.mm))},
                            error=function(e){NaN})
    zip_mse2[k] <-tryCatch({(sum((log2(X_ori+1)-log2(re_zip$muz+1))^2,na.rm=T)/(n.nn*n.mm))},
                           error=function(e){NaN})
    zinb_cov_mse2[k] <-tryCatch({(sum((log2(X_ori+1)-log2(re_zinb_cov$muz+1))^2,na.rm=T)/(n.nn*n.mm))},
                                error=function(e){NaN})
    zip_cov_mse2[k] <-tryCatch({(sum((log2(X_ori+1)-log2(re_zip_cov$muz+1))^2,na.rm=T)/(n.nn*n.mm))},
                               error=function(e){NaN})
    sav_mse2[k] <-tryCatch({(sum((log2(X_ori+1)-log2(t(sav$estimate)+1))^2,na.rm=T)/(n.nn*n.mm))},
                           error=function(e){NaN})
    mbi_mse2[k] <-tryCatch({(sum((log2(X_ori+1)-log2(mbi$imp_count_mat_norm+1))^2,na.rm=T)/(n.nn*n.mm))},
                           error=function(e){NaN})
    ppcanb_mse2[k] <-tryCatch({(sum((log2(X_ori+1)-log2(ppcanb+1))^2,na.rm=T)/(n.nn*n.mm))},
                              error=function(e){NaN})
    ppcanb_cov_mse2[k] <-tryCatch({(sum((log2(X_ori+1)-log2(ppcanb_cov+1))^2,na.rm=T)/(n.nn*n.mm))},
                                  error=function(e){NaN})
  ##wasserstein1d
    complete_mean[[k]]<- colMeans(X_ori,na.rm = T)
    complete_sd[[k]] <-  apply(na.omit(X_ori),2,sd)
    noimpute_mean[[k]] <- colMeans(X,na.rm = T)
    noimpute_sd[[k]] <-  apply(na.omit(X),2,sd)
    noimpute_wasserstein1d[k] <- wasserstein1d(colMeans(X_ori,na.rm = T)/apply(na.omit(X_ori),2,sd), colMeans(X,na.rm = T)/apply(na.omit(X),2,sd))
    
    zinb_mean[[k]] <- tryCatch({colMeans(re_zinb$muz,na.rm = T)},error=function(e){NaN})
    zinb_sd[[k]] <-  tryCatch({apply(na.omit(re_zinb$muz),2,sd)},error=function(e){NaN})
    zinb_wasserstein1d[k] <- tryCatch({wasserstein1d(colMeans(X_ori,na.rm = T)/apply(na.omit(X_ori),2,sd), colMeans(re_zinb$muz,na.rm = T)/apply(na.omit(re_zinb$muz),2,sd))},error=function(e){NaN})
    
    zip_mean[[k]] <- tryCatch({colMeans(re_zip$muz,na.rm = T)},error=function(e){NaN})
    zip_sd[[k]] <-  tryCatch({apply(na.omit(re_zip$muz),2,sd)},error=function(e){NaN})
    zip_wasserstein1d[k] <- tryCatch({wasserstein1d(colMeans(X_ori,na.rm = T)/apply(na.omit(X_ori),2,sd), colMeans(re_zip$muz,na.rm = T)/apply(na.omit(re_zip$muz),2,sd))},error=function(e){NaN})
    
    zinb_cov_mean[[k]] <- tryCatch({colMeans(re_zinb_cov$muz,na.rm = T)},error=function(e){NaN})
    zinb_cov_sd[[k]] <-  tryCatch({apply(na.omit(re_zinb_cov$muz),2,sd)},error=function(e){NaN})
    zinb_cov_wasserstein1d[k] <- tryCatch({wasserstein1d(colMeans(X_ori,na.rm = T)/apply(na.omit(X_ori),2,sd), colMeans(re_zinb_cov$muz,na.rm = T)/apply(na.omit(re_zinb_cov$muz),2,sd))},error=function(e){NaN})
    
    zip_cov_mean[[k]] <- tryCatch({colMeans(re_zip_cov$muz,na.rm = T)},error=function(e){NaN})
    zip_cov_sd[[k]] <-  tryCatch({apply(na.omit(re_zip_cov$muz),2,sd)},error=function(e){NaN})
    zip_cov_wasserstein1d[k] <- tryCatch({wasserstein1d(colMeans(X_ori,na.rm = T)/apply(na.omit(X_ori),2,sd), colMeans(re_zip_cov$muz,na.rm = T)/apply(na.omit(re_zip_cov$muz),2,sd))},error=function(e){NaN})
    
    ppcanb_mean[[k]] <- tryCatch({colMeans(ppcanb,na.rm = T)},error=function(e){NaN})
    ppcanb_sd[[k]] <-  tryCatch({apply(na.omit(ppcanb),2,sd)},error=function(e){NaN})
    ppcanb_wasserstein1d[k] <- tryCatch({wasserstein1d(colMeans(X_ori,na.rm = T)/apply(na.omit(X_ori),2,sd), colMeans(ppcanb,na.rm = T)/apply(na.omit(ppcanb),2,sd))},error=function(e){NaN})
    
    ppcanb_cov_mean[[k]] <- tryCatch({colMeans(ppcanb_cov,na.rm = T)},error=function(e){NaN})
    ppcanb_cov_sd[[k]] <-  tryCatch({apply(na.omit(ppcanb_cov),2,sd)},error=function(e){NaN})
    ppcanb_cov_wasserstein1d[k] <- tryCatch({wasserstein1d(colMeans(X_ori,na.rm = T)/apply(na.omit(X_ori),2,sd), colMeans(ppcanb_cov,na.rm = T)/apply(na.omit(ppcanb_cov),2,sd))},error=function(e){NaN})
    
    if(is.null(mbi)){
      mbi_mean[[k]] <- NaN;mbi_sd[[k]] <- NaN;mbi_wasserstein1d[k] <- NaN;
    }else{
      mbi_mean[[k]] <- colMeans(mbi$imp_count_mat_norm,na.rm = T)
      mbi_sd[[k]] <-  apply(na.omit(mbi$imp_count_mat_norm),2,sd)
      mbi_wasserstein1d[k] <- wasserstein1d(colMeans(X_ori,na.rm = T)/apply(na.omit(X_ori),2,sd),colMeans(mbi$imp_count_mat_norm,na.rm = T)/apply(na.omit(mbi$imp_count_mat_norm),2,sd))
    }
    
    if(is.null(sav)){
      sav_mean[[k]] <- NaN;sav_sd[[k]] <- NaN;sav_wasserstein1d[k] <- NaN;
    }else{
      sav_mean[[k]] <- colMeans(t(sav$estimate),na.rm = T)
      sav_sd[[k]] <-  apply(na.omit(t(sav$estimate)),2,sd)
      sav_wasserstein1d[k] <- wasserstein1d(colMeans(X_ori,na.rm = T)/apply(na.omit(X_ori),2,sd), colMeans(t(sav$estimate),na.rm = T)/apply(na.omit(t(sav$estimate)),2,sd))
    }
  
  ##cor
  for(j in 1:ncol(X)){
    noimpute_cor[[k]][j] <- tryCatch({cor(X_ori[,j],X[,j],use = "na.or.complete")},
                                     error=function(e){NaN})
    zinb_cor[[k]][j] <-tryCatch({cor(X_ori[,j],re_zinb$muz[,j],use = "na.or.complete")},
                                error=function(e){NaN})
    zip_cor[[k]][j] <-tryCatch({cor(X_ori[,j],re_zip$muz[,j],use = "na.or.complete")},
                               error=function(e){NaN})
    zinb_cov_cor[[k]][j] <-tryCatch({cor(X_ori[,j],re_zinb_cov$muz[,j],use = "na.or.complete")},
                                    error=function(e){NaN})
    zip_cov_cor[[k]][j] <-tryCatch({cor(X_ori[,j],re_zip_cov$muz[,j],use = "na.or.complete")},
                                   error=function(e){NaN})
    sav_cor[[k]][j] <-tryCatch({cor(X_ori[,j],(t(sav$estimate))[,j],use = "na.or.complete")},
                               error=function(e){NaN})
    mbi_cor[[k]][j] <-tryCatch({cor(X_ori[,j],mbi$imp_count_mat_norm[,j],use = "na.or.complete")},
                               error=function(e){NaN})
    ppcanb_cor[[k]][j] <-tryCatch({cor(X_ori[,j],ppcanb[,j],use = "na.or.complete")},
                                  error=function(e){NaN})
    ppcanb_cov_cor[[k]][j] <-tryCatch({cor(X_ori[,j],ppcanb_cov[,j],use = "na.or.complete")},
                                      error=function(e){NaN})
    
  }
  
 
}

save.image("simulation/sim_zifa1.RData")
