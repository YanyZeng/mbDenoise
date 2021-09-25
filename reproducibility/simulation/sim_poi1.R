setwd("~/mbDenoise/reproducibility")

n.n = 60
n.w = 100
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
  Q.iter <-1
  n.mm <- 0
  X <- matrix(0,n.n,n.w,byrow = TRUE)
  while(n.mm!=n.w | sum(is.na(X)) >0){
    if(Q.iter>1){
      s <- sample(100:2000,1)
      set.seed(k+s)
      print(k+s)
    }else{
      set.seed(k)
      print(k)
    }
    f <- matrix(0,nrow = n.n, ncol = n.factors)
    for(i in 1:n.n){
      f[i,] <- rnorm(n.factors, mean = 0, sd = 1)
    }
    betaj <- matrix(0,nrow = n.w, ncol = n.factors)
    for(j in 1:n.w){
      betaj[j,] <- runif(n.factors,-3,3)
    }
    alpha <- runif(n.n,-5,5)
    beta0 <- rep(0,n.w)
    g <- rep(0,n.w*0.5*0.5)
    ## in DA test, 0= 1,2,3,4,5, such as g <- rep(1,n.w*0.5*0.5)
    gamma <- c(g,-g,rep(0,(n.w-n.w*0.5)))
    X_cov<- c(rep(1,n.n/2),rep(0,n.n/2))
    
    ll <- f %*% t(betaj) + matrix(gamma,n.n,n.w,byrow=TRUE)*X_cov +matrix(alpha,n.n,n.w)+matrix(beta0,n.n,n.w,byrow=TRUE)
    exp_mat <- exp(ll)
    eta_mat <- matrix(0.25,n.n,n.w,byrow=TRUE)
    
    z <- matrix(0,n.n,n.w,byrow = TRUE)
    for(i in 1:n.n){
      z[i,] <- rbinom(n.w, size=1, prob=eta_mat[i,])
    }
    sum <- rowSums((1-z)*exp_mat)
    Qn_z <- (1-z)*exp_mat/sum
    
    sum <- rowSums(exp_mat)
    Qn <- exp_mat/sum
    
    X <- matrix(0,n.n,n.w,byrow = TRUE)
    for(i in 1:n.n){
      for(j in 1:n.w){
        X[i,j] <- rpois(n=1,lambda = exp_mat[i,j])
      }
    }
    X[z==1]=0
    colnames(X) <- c(1:ncol(X))
    
    Y <- log2(1+X)
    X_ori <- (1-z)*exp(f %*% t(betaj)+ matrix(gamma,n.n,n.w,byrow=TRUE)*X_cov +matrix(beta0,n.n,n.w,byrow=TRUE))
    
    zerorow <- which(rowSums(X)==0)
    if(length(zerorow) >0 ){
      X <- X[-zerorow,];X_ori <- X_ori[-zerorow,];Y <- Y[-zerorow,];X_cov<-X_cov[-zerorow];
      f <- f[-zerorow,];Qn <- Qn[-zerorow,];z <- z[-zerorow,]
    }
    zerocol <- which(colSums(X)==0)
    if(length(zerocol) >0 ){
      X <- X[,-zerocol];X_ori <- X_ori[,-zerocol];Y <- Y[,-zerocol];
      betaj <- t(t(betaj)[,-zerocol]);Qn <- Qn[,-zerocol];z <- z[,-zerocol]
    }
    
    n.nn <- nrow(X)
    n.mm <- ncol(X)
    
    depth[k] <- summary(rowSums(X))[6]/summary(rowSums(X))[1]
    zero.prob[k]  <- sum(X==0)/(dim(X)[1]*dim(X)[2])
    zeroinfl.prob[k] <- sum(z==1)/(sum(X==0))
    
    Q.iter <- Q.iter+1
  }
  
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
  
  # pmr
  t7 <- Sys.time()
  Xhat_pg <- tryCatch({autoTuneProxGradient(W = X)$X_hat},
                      error=function(e){NULL})
  t8 <- Sys.time()
  difftime_pg[k] <- difftime(t8,t7,units = "secs")
  
  # zr
  t7 <- Sys.time()
  Xhat_zr05 <- zr(W = X, alpha = 0.5) # replaced by 0.5
  t8 <- Sys.time()
  difftime_zr[k] <- difftime(t8,t7,units = "secs")
  
  # svt
  t7 <- Sys.time()
  Xhat_svt05 <- svt(W = X, r = 20, alpha = 0.5) # replaced by 0.5
  t8 <- Sys.time()
  difftime_svt[k] <- difftime(t8,t7,units = "secs")
  
  ##dmm
  t7 <- Sys.time()
  dmn <- tryCatch({DirichletMultinomial::dmn(X,k=10)},
                  error=function(e){NULL})
  t8 <- Sys.time()
  difftime_dmn[k] <- difftime(t8,t7,units = "secs")

  alpha3 <- dmn@fit$Estimate
  comesti3_s <- list()
  for (i in 1:ncol(alpha3)) {
    comesti3_s[[i]] <- t(apply(X,1,function(x){ (x+alpha3[,i])/(sum(x+alpha3[,i]))}))*dmn@group[,i]
  }
  QQQ_dmn <- Reduce('+', comesti3_s)
  
  QQQ_mbi <- tryCatch({(zr(W =mbi$imp_count_mat_lognorm, alpha = 0.5) )},
                      error=function(e){ NaN})
  QQQ_sav<-tryCatch({zr(W =t(sav$estimate), alpha = 0.5) },
                    error=function(e){ NaN})
  ppcanb_cov_q <- tryCatch({exp(re_nbva_cov$lvs %*% t(re_nbva_cov$params$theta) + matrix(re_nbva_cov$params$Xcoef,n.nn,n.mm,byrow=TRUE)*X_cov+matrix(re_nbva_cov$params$beta0,n.nn,n.mm,byrow=TRUE))/rowSums(exp(re_nbva_cov$lvs %*% t(re_nbva_cov$params$theta) + matrix(re_nbva_cov$params$Xcoef,n.nn,n.mm,byrow=TRUE)*X_cov+matrix(re_nbva_cov$params$beta0,n.nn,n.mm,byrow=TRUE)))},
                           error=function(e){NaN})
  ppcanb_q <- tryCatch({exp(re_nbva$lvs %*% t(re_nbva$params$theta) + matrix(re_nbva$params$beta0,n.nn,n.mm,byrow=TRUE))/rowSums(exp(re_nbva$lvs %*% t(re_nbva$params$theta) + matrix(re_nbva$params$beta0,n.nn,n.mm,byrow=TRUE)))},
                       error=function(e){NaN})
  
  ppcanb_cov <- tryCatch({exp(re_nbva_cov$lvs %*% t(re_nbva_cov$params$theta) + matrix(re_nbva_cov$params$Xcoef,n.nn,n.mm,byrow=TRUE)*X_cov+matrix(re_nbva_cov$params$beta0,n.nn,n.mm,byrow=TRUE))},
                         error=function(e){NaN})
  ppcanb <- tryCatch({exp(re_nbva$lvs %*% t(re_nbva$params$theta) + matrix(re_nbva$params$beta0,n.nn,n.mm,byrow=TRUE))},
                     error=function(e){NaN})
  
  # ###1. frobenius norm error
  if(is.null(re_zinb)){a_nb[k] <- NaN}else{a_nb[k] <- sqrt((sum((re_zinb$Q-Qn)^2)))}
  if(is.null(re_zip)){a_poi[k] <- NaN}else{a_poi[k] <- sqrt((sum((re_zip$Q-Qn)^2)))}
  if(is.null(Xhat_pg)){a_multi[k] <- NaN}else{a_multi[k] <- sqrt((sum((Xhat_pg-Qn)^2)))}
  if(is.null(mbi)){a_mbi[k] <- NaN}else{a_mbi[k] <- sqrt((sum((QQQ_mbi-Qn)^2)))}
  if(is.null(sav)){a_sav[k] <- NaN}else{a_sav[k] <- sqrt((sum((QQQ_sav-Qn)^2)))}
  if(is.null(re_zinb_cov)){a_nb_cov[k] <- NaN}else{a_nb_cov[k] <- sqrt((sum((re_zinb_cov$Q-Qn)^2)))}
  if(is.null(re_zip_cov)){a_poi_cov[k] <- NaN}else{a_poi_cov[k] <- sqrt((sum((re_zip_cov$Q-Qn)^2)))}
  if(is.null(dmn)){a_dmn[k] <- NaN}else{a_dmn[k] <- sqrt((sum((QQQ_dmn-Qn)^2)))}
  if(is.null(Xhat_zr05)){a_0.5[k] <- NaN}else{a_0.5[k] <- sqrt((sum((Xhat_zr05-Qn)^2)))}
  if(is.null(Xhat_svt05)){a_svt[k] <- NaN}else{a_svt[k] <- sqrt((sum((Xhat_svt05-Qn)^2)))}
  if(is.null(re_nbva$params$theta) || is.null(re_nbva$lvs)){a_ppcanb[k] <- NaN}else{a_ppcanb[k] <- sqrt((sum((ppcanb_q-Qn)^2)))}
  if(is.null(re_nbva_cov$params$theta) || is.null(re_nbva_cov$lvs)){a_ppcanb_cov[k] <- NaN}else{a_ppcanb_cov[k] <- sqrt((sum((ppcanb_cov_q-Qn)^2)))}
  
  # ### 2. KL
  if(is.null(re_zinb)){b_nb[k] <- NaN}else{b_nb[k] <-  sum(Qn*log(Qn/re_zinb$Q))/n.nn}
  if(is.null(re_zip)){b_poi[k] <- NaN}else{b_poi[k] <- sum(Qn*log(Qn/re_zip$Q))/n.nn}
  if(is.null(Xhat_pg)){b_multi[k] <- NaN}else{b_multi[k] <- sum(Qn*log(Qn/Xhat_pg))/n.nn}
  if(is.null(mbi)){b_mbi[k] <- NaN}else{b_mbi[k] <-  sum(Qn*log(Qn/QQQ_mbi))/n.nn}
  if(is.null(sav)){b_sav[k] <- NaN}else{b_sav[k] <-  sum(Qn*log(Qn/(QQQ_sav)))/n.nn}
  if(is.null(re_zinb_cov)){b_nb_cov[k] <- NaN}else{b_nb_cov[k] <-  sum(Qn*log(Qn/re_zinb_cov$Q))/n.nn}
  if(is.null(re_zip_cov)){b_poi_cov[k] <- NaN}else{b_poi_cov[k] <- sum(Qn*log(Qn/re_zip_cov$Q))/n.nn}
  if(is.null(dmn)){b_dmn[k] <- NaN}else{ b_dmn[k] <-  sum(Qn*log(Qn/QQQ_dmn))/n.n}
  if(is.null(Xhat_zr05)){b_0.5[k] <- NaN}else{b_0.5[k] <- sum(Qn*log(Qn/Xhat_zr05))/n.nn}
  if(is.null(Xhat_svt05)){b_svt[k] <- NaN}else{b_svt[k] <- sum(Qn*log(Qn/Xhat_svt05))/n.nn}
  if(is.null(re_nbva$params$theta) || is.null(re_nbva$lvs)){b_ppcanb[k] <- NaN}else{b_ppcanb[k] <-  sum(Qn*log(Qn/ppcanb_q))/n.nn}
  if(is.null(re_nbva_cov$params$theta) || is.null(re_nbva_cov$lvs)){b_ppcanb_cov[k] <- NaN}else{b_ppcanb_cov[k] <-  sum(Qn*log(Qn/ppcanb_cov_q))/n.nn}
  
  # ## 3. shannon index mse
  if(is.null(re_zinb)){c_nb[k] <- NaN}else{c_nb[k] <- (sum((-rowSums(Qn*log(Qn))+rowSums(re_zinb$Q*log(re_zinb$Q)))^2))/n.nn}
  if(is.null(re_zip)){c_poi[k] <- NaN}else{c_poi[k] <- (sum((-rowSums(Qn*log(Qn))+rowSums(re_zip$Q*log(re_zip$Q)))^2))/n.nn}
  if(is.null(Xhat_pg)){c_multi[k] <- NaN}else{c_multi[k] <- (sum((-rowSums(Qn*log(Qn))+rowSums(Xhat_pg*log(Xhat_pg)))^2))/n.nn}
  if(is.null(mbi)){c_mbi[k] <- NaN}else{c_mbi[k] <-  (sum((-rowSums(Qn*log(Qn))+rowSums(QQQ_mbi*log((QQQ_mbi))))^2))/n.nn}
  if(is.null(sav)){c_sav[k] <- NaN}else{c_sav[k] <-  (sum((-rowSums(Qn*log(Qn))+rowSums((QQQ_sav)*log((QQQ_sav))))^2))/n.nn}
  if(is.null(re_zinb_cov)){c_nb_cov[k] <- NaN}else{c_nb_cov[k] <- (sum((-rowSums(Qn*log(Qn))+rowSums(re_zinb_cov$Q*log(re_zinb_cov$Q)))^2))/n.nn}
  if(is.null(re_zip_cov)){c_poi_cov[k] <- NaN}else{c_poi_cov[k] <- (sum((-rowSums(Qn*log(Qn))+rowSums(re_zip_cov$Q*log(re_zip_cov$Q)))^2))/n.nn}
  if(is.null(dmn)){c_dmn[k] <- NaN}else{ c_dmn[k] <-(sum((-rowSums(Qn*log(Qn))+rowSums(QQQ_dmn*log(QQQ_dmn)))^2))/n.n}
  if(is.null(Xhat_zr05)){c_0.5[k] <- NaN}else{c_0.5[k] <- (sum((-rowSums(Qn*log(Qn))+rowSums(Xhat_zr05*log(Xhat_zr05)))^2))/n.nn}
  if(is.null(Xhat_svt05)){c_svt[k] <- NaN}else{c_svt[k] <- (sum((-rowSums(Qn*log(Qn))+rowSums(Xhat_svt05*log(Xhat_svt05)))^2))/n.nn}
  if(is.null(re_nbva$params$theta) || is.null(re_nbva$lvs)){c_ppcanb[k] <- NaN}else{c_ppcanb[k] <- (sum((-rowSums(Qn*log(Qn))+rowSums(ppcanb_q*log(ppcanb_q)))^2))/n.nn}
  if(is.null(re_nbva_cov$params$theta) || is.null(re_nbva_cov$lvs)){c_ppcanb_cov[k] <- NaN}else{c_ppcanb_cov[k] <- (sum((-rowSums(Qn*log(Qn))+rowSums(ppcanb_cov_q*log(ppcanb_cov_q)))^2))/n.nn}
  
  # ###4.simpson index mse
  if(is.null(re_zinb)){d_nb[k] <- NaN}else{d_nb[k] <- (sum((rowSums(Qn^2)-rowSums((re_zinb$Q^2)))^2))/n.nn}
  if(is.null(re_zip)){d_poi[k] <- NaN}else{d_poi[k] <- (sum((rowSums(Qn^2)-rowSums((re_zip$Q^2)))^2))/n.nn}
  if(is.null(Xhat_pg)){d_multi[k] <- NaN}else{d_multi[k] <- (sum((rowSums(Qn^2)-rowSums((Xhat_pg^2)))^2))/n.nn}
  if(is.null(mbi)){d_mbi[k] <- NaN}else{d_mbi[k] <- (sum((rowSums(Qn^2)-rowSums((QQQ_mbi)))^2))/n.nn}
  if(is.null(sav)){d_sav[k] <- NaN}else{d_sav[k] <- (sum((rowSums(Qn^2)-rowSums(((QQQ_sav)^2)))^2))/n.nn}
  if(is.null(re_zinb_cov)){d_nb_cov[k] <- NaN}else{d_nb_cov[k] <- (sum((rowSums(Qn^2)-rowSums((re_zinb_cov$Q^2)))^2))/n.nn}
  if(is.null(re_zip_cov)){d_poi_cov[k] <- NaN}else{d_poi_cov[k] <- (sum((rowSums(Qn^2)-rowSums((re_zip_cov$Q^2)))^2))/n.nn}
  if(is.null(dmn)){d_dmn[k] <- NaN}else{ d_dmn[k] <-(sum((rowSums(Qn^2)-rowSums((QQQ_dmn^2)))^2))/n.n}
  if(is.null(Xhat_zr05)){d_0.5[k] <- NaN}else{d_0.5[k] <- (sum((rowSums(Qn^2)-rowSums((Xhat_zr05^2)))^2))/n.nn}
  if(is.null(Xhat_svt05)){d_svt[k] <- NaN}else{d_svt[k] <- (sum((rowSums(Qn^2)-rowSums((Xhat_svt05^2)))^2))/n.nn}
  if(is.null(re_nbva$params$theta) || is.null(re_nbva$lvs)){d_ppcanb[k] <- NaN}else{d_ppcanb[k] <- (sum((rowSums(Qn^2)-rowSums((ppcanb_q^2)))^2))/n.nn}
  if(is.null(re_nbva_cov$params$theta) || is.null(re_nbva_cov$lvs)){d_ppcanb_cov[k] <- NaN}else{d_ppcanb_cov[k] <- (sum((rowSums(Qn^2)-rowSums((ppcanb_cov_q^2)))^2))/n.nn}
  
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
  
  ##DESeq2
  group_gile <- data.frame(cbind(rownames(c(1:nrow(X))),X_cov))
  colnames(group_gile) <- c("X_cov")
  dds <- tryCatch({DESeqDataSetFromMatrix(countData = t(X+1), colData = group_gile, design = ~X_cov)},
                  error=function(e){ NULL})
  if(is.null(dds)){
    deseq_p_adj[,k] <- NaN;
  }else{
    dds <- DESeq(dds,quiet=T)
    res <- results(dds)#提取分析结果
    deseq_p_adj[,k] <-res$padj}
  
  ##edgeR
  X_tr <- data.frame(t(X))
  y<-DGEList(counts=X_tr,genes=rownames(X_tr))
  y<-edgeR::calcNormFactors(y)
  design <- model.matrix(~X_cov)
  Disp <- estimateDisp(y, design, robust = TRUE)
  fit <- glmFit(Disp, design)
  qlf <- glmLRT(fit, coef=2)
  edger_p_adj[,k] <- p.adjust(qlf$table$PValue, method = "BH")
  
  ##metagenomeSeq:fitFeatureModel
  data_c <- t(X)
  rownames(data_c) <- NULL
  samp <- data.frame(X_cov)
  phenotypeData = AnnotatedDataFrame(samp)
  X_MRexperiment <- metagenomeSeq::newMRexperiment(data_c,phenoData=phenotypeData,featureData=NULL)
  XX <- cumNorm(X_MRexperiment, p = 0.5)
  pd <- pData(XX)
  mod <- model.matrix(~1+X_cov, data = pd)
  fit <- tryCatch({fitFeatureModel(XX, mod)},
                  error=function(e){ NULL})
  if(is.null(fit)){
    meta_p_adj[,k] <- NaN;
  }else{
    meta <- fit@pvalues
    meta_p_adj[,k] <- p.adjust(meta,method = "fdr")
  }
  
  for(l in 1:ncol(X)){
    zip_p[l,k] <- tryCatch({t.test(log2(re_zip$muz+1)[,l]~X_cov)$p.value},error=function(e){ NaN})
    zinb_p[l,k] <- tryCatch({t.test(log2(re_zinb$muz+1)[,l]~X_cov)$p.value},error=function(e){ NaN})
    zip_cov_p[l,k] <- tryCatch({t.test(log2(re_zip_cov$muz+1)[,l]~X_cov)$p.value},error=function(e){ NaN})
    zinb_cov_p[l,k] <- tryCatch({t.test(log2(re_zinb_cov$muz+1)[,l]~X_cov)$p.value},error=function(e){ NaN})
    saver_p[l,k] <- tryCatch({t.test((t(log2(sav$estimate+1)))[,l]~X_cov)$p.value},error=function(e){ NaN})
    mbimpute_p[l,k] <- tryCatch({t.test((log2(mbi$imp_count_mat_norm+1))[,l]~X_cov)$p.value},error=function(e){ NaN})
    real_p2[l,k]<- tryCatch({t.test(log2(X[,l]+1)~X_cov)$p.value},error=function(e){ NaN})
    ppcanb_p[l,k] <- tryCatch({t.test(log2(ppcanb+1)[,l]~X_cov)$p.value},error=function(e){ NaN})
    ppcanb_cov_p[l,k] <- tryCatch({t.test(log2(ppcanb_cov+1)[,l]~X_cov)$p.value},error=function(e){ NaN})
    
  }
  ppcanb_p_adj[,k] <-  p.adjust(ppcanb_p[,k],method = "BH")
  ppcanb_cov_p_adj[,k] <-  p.adjust(ppcanb_cov_p[,k],method = "BH")
  real_p_adj2[,k] <-  p.adjust(real_p2[,k],method = "BH")
  zinb_p_adj[,k] <- p.adjust(zinb_p[,k],method = "BH")
  zip_p_adj[,k] <- p.adjust(zip_p[,k],method = "BH")
  zinb_cov_p_adj[,k] <- p.adjust(zinb_cov_p[,k],method = "BH")
  zip_cov_p_adj[,k] <- p.adjust(zip_cov_p[,k],method = "BH")
  saver_p_adj[,k] <- p.adjust(saver_p[,k],method = "BH")
  mbimpute_p_adj[,k] <- p.adjust(mbimpute_p[,k],method = "BH")
  
}

save.image("simulation/sim_poi1.RData")
