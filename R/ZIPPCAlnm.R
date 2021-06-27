#' @title ZIPPCAlnm
#' @description Microbiome data denoising framework (mbDenoise) with zero-inflated probabilistic PCA with logistical normal multinomial model (ZIPPCA-LNM),
#'              which can be used for downstream statistical analysis including ordination, compositional
#'              normalization, differential abundance analysis, etc.
#' @param X matrix of observations.
#' @param V vector of the sample covariate.
#' @param n.factors the rank or number of factors, after dimensional reduction. Defaults to 2.
#' @param rank logical, if TRUE, the rank or number of factors, is chosen from 1 to 5 by BIC (Bayesian information criterion). Defaults to FALSE.
#' @param trace logical, defaults to \code{FALSE}. if \code{TRUE} each current iteration step information will be printed.
#' @param maxit maximum number of iterations within \code{optim} and \code{constrOptim} function, defaults to 100.
#' @param parallel logical, if TRUE, use parallel toolbox to accelerate.

#' @return
#'
#'
#'  \item{VLB }{ variational lower bound of log likelihood}
#'  \item{lvs}{list of latent variables
#'  \itemize{
#'    \item{pi }{ the probabilities of excess zeros}
#'    \item{factor_scores }{ coordinates or factor scores in low-dimensional subspace}
#'    \item{factor_scores2 }{ coordinates or factor scores in low-dimensional subspace with defalt rank 2, which is suitable for visualization.}
#'    }}
#'  \item{params}{list of model parameters
#'  \itemize{
#'    \item{factor_coefs_j }{ coefficients of latent variables fator scores or factor loadings}
#'    \item{factor_coefs_0 }{ taxon-specific intercepts}
#'    \item{gamma }{ coeffcients of sample covariate}
#'    \item{tuo }{ taxon-specific parameter of zero-inflation probability}
#'    \item{c }{ sample-specific parameter of zero-inflation probability}
#'    }}
#'  \item{Q }{ the underlying composition of microbiome data}
#'  \item{muz }{ the denoised counts of microbiome data}
#'  \item{bic}{ the number of the rank selection, chosen by BIC type information criterion}

#' @examples
#' n.n = 60
#' n.w = 100
#' n.factors = 2
#' set.seed(1)
#' si <- diag(n.factors)
#' me <- c(0,0)
#' f <- matrix(0,nrow = n.n, ncol = n.factors)
#' for(i in 1:n.n){
#'  f[i,] <- rnorm(n.factors, mean = 0, sd = 1)
#' }
#' betaj <- matrix(0,nrow = n.w, ncol = n.factors)
#' for(j in 1:n.w){
#'   betaj[j,] <- runif(n.factors,-1,1)
#' }
#' alpha <- rep(0,n.n)
#' beta0 <- rep(0,n.w)
#' g <- rep(0,n.w*0.5*0.5)
#' gamma <- c(g,-g,rep(0,(n.w-n.w*0.5)))
#' X_cov<- c(rep(1,n.n/2),rep(0,n.n/2))
#' ll <- f %*% t(betaj) +matrix(alpha,n.n,n.w)+matrix(beta0,n.n,n.w,byrow=TRUE)
#' exp_mat <- exp(ll)
#' eta_mat <- matrix(0.25,n.n,n.w,byrow=TRUE)
#' z <- matrix(0,n.n,n.w,byrow = TRUE)
#' for(i in 1:n.n){
#'   z[i,] <- rbinom(n.w, size=1, prob=eta_mat[i,])
#' }
#' sum <- rowSums((1-z)*exp_mat)
#' Qn_z <- (1-z)*exp_mat/sum
#' sum <- rowSums(exp_mat)
#' Qn <- exp_mat/sum
#' X <- matrix(0,n.n,n.w,byrow = TRUE)
#' for(i in 1:n.n){
#' X[i,] <- rmultinom(1, size = runif(1,80,800), prob = Qn[i,])
#' }
#' X[z==1]=0

#' zerorow <- which(rowSums(X)==0)
#' if(length(zerorow) >0 ){
#'    X <- X[-zerorow,];X_cov<-X_cov[-zerorow];f <- f[-zerorow,];
#'    Qn <- Qn[-zerorow,];Qn_z <- Qn_z[-zerorow,];
#' }
#' zerocol <- which(colSums(X)==0)
#' if(length(zerocol) >0 ){
#'   X <- X[,-zerocol];betaj <- t(t(betaj)[,-zerocol]);
#'   Qn <- Qn[,-zerocol];Qn_z <- Qn_z[,-zerocol];
#' }
#' re_zilnm_cov <- ZIPPCAlnm(X,X_cov)
#' re_zilnm <- ZIPPCAlnm(X)
#'
#' @export

ZIPPCAlnm <- function(X, V=NULL, n.factors=2, rank=FALSE,
                      trace = FALSE, maxit = 100, parallel=TRUE){

  ZILNMVA <- function(X,V,n.factors,trace,maxit) {

    n.s<-nrow(X); n.f<-ncol(X);
    M <- rowSums(X)
    if(is.null(V)){Y <- 0
    }else if(is.numeric(V)){Y <- V
    }else{
      Y <- as.numeric(as.factor(V))-1}

    out.list <- list()

    ### Initialization 1
    pzero.col <- apply(X, 2, function(x) {sum(x==0)/n.s})
    z.hat <- pi <- new.pi <- t(ifelse(t(X)==0, pzero.col, 0))
    tuo_j <- car::logit(colMeans(pi))
    c_i <- car::logit(rowMeans(pi))
    #eta <- new.eta <- round(apply(pi, 2, mean), 6)
    sigma <- new.sigma <- matrix(1,n.s,n.factors)
    factor_coefs_0 <-  new.factor_coefs_0 <-  rep(1,n.f)
    gamma <-  new.gamma <-  rep(1,n.f)

    X.rc <- scale(log(X+0.05),scale = T,center = T)
    re <- svd(X.rc,n.factors,n.factors)
    factor_coefs_j <- new.factor_coefs_j <- re$v
    if(n.factors==1){
      factor_scores <- new.factor_scores <- re$u * (re$d[1])
    }else{factor_scores <- new.factor_scores <- re$u %*% diag(re$d[1:n.factors])}


    ### Initialization 2
    cur.VLB <- -1e6; iter <- 1; ratio <- 10; diff=1e5;eps = 1e-4;max.iter = 100;
    b.cur.logfunc <- -1e6;b0.cur.logfunc <- -1e6;f.cur.logfunc <- -1e6;ga.cur.logfunc <- -1e6

    while((diff> eps*(abs(cur.VLB)+eps)) && iter <= max.iter) {
      if(trace) cat("Iteration:", iter, "\n")
      ## VLB
      VLB_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*ll
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))
        fun2 <- function(i) { 0.5 * (sum(log(t(new.sigma)[,i])) - sum(t(new.sigma)[,i]) - sum(new.factor_scores[i,]^2))}

        y <- sum(y1)+sum(y2)+sum(sapply(1:n.s,fun2))
        return(y)
      }

      ###optim f
      f_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        lf <- new.factor_scores %*% t(new.factor_coefs_j)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*lf
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))
        fun2 <- function(i) { 0.5 * (- sum(new.factor_scores[i,]^2))}
        y <- sum(y1)+sum(y2)+sum(sapply(1:n.s,fun2))

        return(y)
      }
      f_grad_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        sum <- exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))
        fun2 <- function(i) { -new.factor_scores[i,]+(X[i,])%*%new.factor_coefs_j-(M[i]*sum[i,])%*%new.factor_coefs_j}
        f_grad <- t(sapply(1:n.s,fun2))
        return(c(f_grad))
      }

      q <- try(optim(c(factor_scores), x=new.factor_coefs_0,b=new.factor_coefs_j,s=new.sigma,g=new.gamma, method = "BFGS", fn = f_f_ini, gr = f_grad_f_ini, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_scores <- factor_scores;
      }else{
        if(iter > 1 && f.cur.logfunc > q$value){if(trace)
          cat("Optimization of mu did not improve on iteration step ",iter,"\n");
          new.factor_scores <- factor_scores;
        }else{
          if(trace) cat("Variational parameters mu updated","\n")
          new.factor_scores <- matrix(q$par,n.s,n.factors);
          new.factor_scores <- scale(new.factor_scores)
          if(q$convergence != 0) { if(trace) cat("Optimization of mu did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update sigma
      beta2 <- new.factor_coefs_j^2
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
      sum <- M*exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

      for(i in 1:n.s){
        if(n.factors==1){
          new.sigma[i] <- 1/(1+(apply((sum)[i,]*beta2,2,sum)))
        }else{
          new.sigma[i,] <- 1/(1+(diag(diag(apply((sum)[i,]*beta2,2,sum)))))
        }}

      ###update beta
      beta_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        lf <- new.factor_scores %*% t(new.factor_coefs_j)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)

        y1 <- X*lf
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))
        y <- sum(y1)+sum(y2)

        return(y)
      }
      beta_grad_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        b3 <- NULL
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        sum <- M*exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

        for(p in 1:n.factors){
          b3 <- c(b3,(sweep(X,1,new.factor_scores[,p],"*") -sweep((new.sigma[,p]%*%t(new.factor_coefs_j[,p])),1,new.factor_scores[,p],"+") * sum))
        }

        b3 <- matrix(b3,n.s,n.f*n.factors)
        return(c(colSums(b3)))
      }

      q <- try(optim(c(factor_coefs_j), x=new.factor_coefs_0,f=new.factor_scores,s=new.sigma, g=new.gamma,method = "BFGS", fn = beta_f_ini, gr = beta_grad_f_ini, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_coefs_j <- factor_coefs_j;
      }else{
        if(iter > 1 && b.cur.logfunc > q$value){if(trace)
          cat("Optimization of beta did not improve on iteration step ",iter,"\n");
          new.factor_coefs_j <- factor_coefs_j;
        }else{
          if(trace) cat("Model parameters beta updated","\n")
          new.factor_coefs_j <- matrix(q$par,n.f,n.factors);
          #new.factor_coefs_j <- scale(new.factor_coefs_j)
          if(q$convergence != 0) { if(trace) cat("Optimization of beta did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update beta0
      b0_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        lab <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)

        y1 <- X*lab
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))
        y <- sum(y1)+sum(y2)

        return(y)
      }
      b0_grad_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        grad <- X-M*exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

        return(c(colSums(grad)))
      }

      q <- try(optim(c(factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma,g=new.gamma,method = "BFGS", fn = b0_f_ini, gr = b0_grad_f_ini, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_coefs_0 <- factor_coefs_0
      }else{
        if(iter > 1 && b0.cur.logfunc > q$value){if(trace)
          cat("Optimization of beta did not improve on iteration step ",iter,"\n");
          new.factor_coefs_0 <- factor_coefs_0
        }else{
          if(trace) cat("Model parameters beta updated","\n")
          new.factor_coefs_0 <- q$par;
          if(q$convergence != 0) { if(trace) cat("Optimization of beta did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update gamma
      gam_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        lab <- matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)

        y1 <- X*lab
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))
        y <- sum(y1)+sum(y2)

        return(y)
      }
      gam_grad_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        grad <- X*Y-M*exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*Y

        return(c(colSums(grad)))
      }

      q <- try(optim(c(gamma),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma,x=new.factor_coefs_0,method = "BFGS", fn = gam_f_ini, gr = gam_grad_f_ini, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.gamma <- gamma
      }else{
        if(iter > 1 && ga.cur.logfunc > q$value){if(trace)
          cat("Optimization of gamma did not improve on iteration step ",iter,"\n");
          new.gamma <- gamma
        }else{
          if(trace) cat("Model parameters gamma updated","\n")
          new.gamma <- q$par;
          if(q$convergence != 0) { if(trace) cat("Optimization of gamma did not converge on iteration step ", iter,"\n") }
        }
      }

      q1 <- list(value = beta_f_ini(c(new.factor_coefs_j),x=new.factor_coefs_0,f=new.factor_scores,s=new.sigma,g=new.gamma))
      b.new.logfunc <- q1$value
      b.cur.logfunc <- b.new.logfunc

      q2 <- list(value = f_f_ini(c(new.factor_scores),x=new.factor_coefs_0, b=new.factor_coefs_j,s=new.sigma,g=new.gamma))
      new.f.cur.logfunc <- q2$value
      f.cur.logfunc <- new.f.cur.logfunc

      q3 <- list(value = b0_f_ini(c(new.factor_coefs_0), f=new.factor_scores,b=new.factor_coefs_j,s=new.sigma,g=new.gamma))
      new.b0.cur.logfunc <- q3$value
      b0.cur.logfunc <- new.b0.cur.logfunc

      q4 <- list(value = gam_f_ini(c(new.gamma), f=new.factor_scores,b=new.factor_coefs_j,s=new.sigma,x=new.factor_coefs_0))
      new.ga.cur.logfunc <- q4$value
      ga.cur.logfunc <- new.ga.cur.logfunc

      ## Take values of VLB to define stopping rule
      q <- list(value = VLB_ini(c(new.factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma,g=new.gamma))
      new.VLB <-  tryCatch({q$value},error=function(e){NaN})
      diff=abs(new.VLB-cur.VLB)
      ratio <- abs(new.VLB/cur.VLB);
      if(trace) cat("New VLB:", new.VLB,"cur VLB:", cur.VLB, "Ratio of VLB", ratio, ". Difference in VLB:",diff,"\n")
      cur.VLB <- new.VLB

      if(is.na(cur.VLB)){
        sigma <- new.sigma <- matrix(1,n.s,n.factors)
        factor_coefs_0 <-  new.factor_coefs_0 <-  rep(1,n.f)
        gamma <-  new.gamma <-  rep(1,n.f)
        X.rc <- scale(log(X+0.05),scale = T,center = T)
        re <- svd(X.rc,n.factors,n.factors)
        factor_coefs_j <- new.factor_coefs_j <- re$v
        if(n.factors==1){
          factor_scores <- new.factor_scores <- re$u * (re$d[1])
        }else{factor_scores <- new.factor_scores <- re$u %*% diag(re$d[1:n.factors])}
        iter <- 101

      }else{
        factor_coefs_0 <-new.factor_coefs_0
        factor_coefs_j <- new.factor_coefs_j
        factor_scores <- new.factor_scores
        sigma <- new.sigma
        gamma <- new.gamma
        iter <- iter + 1
      }
    }

    ###VA iteration
    cur.VLB <- -1e6; iter <- 1; ratio <- 10; diff=1e5;eps = 1e-4;max.iter = 100;
    b.cur.logfunc <- -1e6;b0.cur.logfunc <- -1e6;f.cur.logfunc <- -1e6;ga.cur.logfunc <- -1e6
    e.cur.logfunc <- -1e6;

    while((diff> eps*(abs(cur.VLB)+eps)) && iter <= max.iter) {
      if(trace) cat("Iteration:", iter, "\n")
      ## VLB
      VLB <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,e=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];
        new.tuo <- e[1:n.f]; e <- e[-(1:n.f)]
        new.ci <- e[1:n.s]

        new.pi[X==0] <- new.pi[X==0]-1e-8
        new.eta <- matrix(new.tuo,n.s, n.f,byrow = T)+ matrix(new.ci,n.s, n.f)

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*ll
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))
        fun2 <- function(i) { 0.5 * (sum(log(t(new.sigma)[,i])) - sum(t(new.sigma)[,i]) - sum(new.factor_scores[i,]^2))}
        y3 <- new.pi*new.eta-log(1+exp(new.eta))-(1-new.pi)*log(1-new.pi)-I(X==0)*new.pi*log(new.pi)

        y <- sum(y1)+sum(y2)+sum(sapply(1:n.s,fun2))+sum(y3,na.rm=T)
        return(y)
      }

      ###optim pi
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
      #alp <- log(M/(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))))))
      alp <- log(M/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))))))
      new.eta <- matrix(tuo_j,n.s, n.f,byrow = T)+ matrix(c_i,n.s, n.f)
      e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)+matrix(alp,n.s,n.f)
      sum1 <- exp( - exp(e.mat))
      new.pi <- 1/(1+exp(-new.eta)*sum1)
      new.pi[X!=0]=0

      ###optim f
      # # ###optim f
      f_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        lf <- new.factor_scores %*% t(new.factor_coefs_j)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*lf
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

        fun2 <- function(i) { 0.5 * (- sum(new.factor_scores[i,]^2))}
        y <- sum(y1)+sum(y2)+sum(sapply(1:n.s,fun2))

        return(y)
      }
      f_grad_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        sum <- I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))/(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

        # fun2 <- function(i) { -new.factor_scores[i,]+(X[i,])%*%new.factor_coefs_j-(X[i,]*sum[i,])%*%new.factor_coefs_j-
        #                       (new.re[i]*sum[i,])%*%new.factor_coefs_j}
        fun2 <- function(i) { -new.factor_scores[i,]+(X[i,])%*%new.factor_coefs_j-(M[i]*sum[i,])%*%new.factor_coefs_j}
        f_grad <- t(sapply(1:n.s,fun2))
        return(c(f_grad))
      }

      q <- try(optim(c(factor_scores), x=new.factor_coefs_0,b=new.factor_coefs_j,s=new.sigma,pi=new.pi,g=new.gamma, method = "BFGS", fn = f_f, gr = f_grad_f, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_scores <- factor_scores;
      }else{
        if(iter > 1 && f.cur.logfunc > q$value){if(trace)
          cat("Optimization of mu did not improve on iteration step ",iter,"\n");
          new.factor_scores <- factor_scores;
        }else{
          if(trace) cat("Variational parameters mu updated","\n")
          new.factor_scores <- matrix(q$par,n.s,n.factors);
          if(q$convergence != 0) { if(trace) cat("Optimization of mu did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update sigma
      beta2 <- new.factor_coefs_j^2
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
      sum <- M*I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))/(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

      for(i in 1:n.s){
        if(n.factors==1){
          new.sigma[i] <- 1/(1+(apply((sum)[i,]*beta2,2,sum)))
        }else{
          new.sigma[i,] <- 1/(1+(diag(diag(apply((sum)[i,]*beta2,2,sum)))))
        }}


      ###update beta
      beta_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        lf <- new.factor_scores %*% t(new.factor_coefs_j)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)

        y1 <- X*lf
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

        y <- sum(y1)+sum(y2)

        return(y)
      }
      beta_grad_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        b3 <- NULL
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        sum <- M*I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))/(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

        for(p in 1:n.factors){
          b3 <- c(b3,(sweep(X,1,new.factor_scores[,p],"*") -sweep((new.sigma[,p]%*%t(new.factor_coefs_j[,p])),1,new.factor_scores[,p],"+") * sum))
        }

        b3 <- matrix(b3,n.s,n.f*n.factors)
        return(c(colSums(b3)))
      }

      q <- try(optim(c(factor_coefs_j), x=new.factor_coefs_0,f=new.factor_scores,s=new.sigma,pi=new.pi, g=new.gamma, method = "BFGS", fn = beta_f, gr = beta_grad_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_coefs_j <- factor_coefs_j;
      }else{
        if(iter > 1 && b.cur.logfunc > q$value){if(trace)
          cat("Optimization of beta did not improve on iteration step ",iter,"\n");
          new.factor_coefs_j <- factor_coefs_j;
        }else{
          if(trace) cat("Model parameters beta updated","\n")
          new.factor_coefs_j <- matrix(q$par,n.f,n.factors);
          if(q$convergence != 0) { if(trace) cat("Optimization of beta did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update beta0
      b0_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        lab <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)

        y1 <- X*lab
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

        y <- sum(y1)+sum(y2)

        return(y)
      }
      b0_grad_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        grad <- X-M*I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))/(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

        return(c(colSums(grad)))
      }


      q <- try(optim(c(factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma, pi=new.pi,g=new.gamma, method = "BFGS", fn = b0_f, gr = b0_grad_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_coefs_0 <- factor_coefs_0
      }else{
        if(iter > 1 && b0.cur.logfunc > q$value){if(trace)
          cat("Optimization of beta did not improve on iteration step ",iter,"\n");
          new.factor_coefs_0 <- factor_coefs_0
        }else{
          if(trace) cat("Model parameters beta updated","\n")
          new.factor_coefs_0 <- q$par
          if(q$convergence != 0) { if(trace) cat("Optimization of beta did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update gamma
      ga_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        lab <- matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)

        y1 <- X*lab
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))

        y <- sum(y1)+sum(y2)

        return(y)
      }
      ga_grad_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        grad <- X*Y-M*I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))/(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*Y

        return(c(colSums(grad)))
      }


      q <- try(optim(c(gamma),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma, pi=new.pi,x=new.factor_coefs_0, method = "BFGS", fn = ga_f, gr = ga_grad_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.gamma <- gamma
      }else{
        if(iter > 1 && ga.cur.logfunc > q$value){if(trace)
          cat("Optimization of gamma did not improve on iteration step ",iter,"\n");
          new.gamma <- gamma
        }else{
          if(trace) cat("Model parameters gamma updated","\n")
          new.gamma <- q$par
          if(q$convergence != 0) { if(trace) cat("Optimization of gamma did not converge on iteration step ", iter,"\n") }
        }
      }

      ## update tuo,ci
      eta_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL,e=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];
        new.tuo <- e[1:n.f]; e <- e[-(1:n.f)]
        new.ci <- e[1:n.s]

        new.pi[X==0] <- new.pi[X==0]-1e-8
        new.eta <- matrix(new.tuo,n.s, n.f,byrow = T)+ matrix(new.ci,n.s, n.f)

        y1 <- new.pi*new.eta-log(1+exp(new.eta))
        y <- sum(y1)
        return(y)
      }
      eta_grad_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL,e=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];
        new.tuo <- e[1:n.f]; e <- e[-(1:n.f)]
        new.ci <- e[1:n.s]

        new.pi[X==0] <- new.pi[X==0]-1e-8
        new.eta <- matrix(new.tuo,n.s, n.f,byrow = T)+ matrix(new.ci,n.s, n.f)

        grad <- new.pi-exp(new.eta)/(1+exp(new.eta))

        return(c(colSums(grad), rowSums(grad)))
      }
      q <- try(optim(c(tuo_j,c_i),g=new.gamma,b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma, pi=new.pi,x=new.factor_coefs_0, method = "BFGS", fn = eta_f, gr = eta_grad_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.tuo <- tuo_j; new.ci <- c_i;
      }else{
        if(iter > 1 && e.cur.logfunc > q$value){if(trace)
          cat("Optimization of eta did not improve on iteration step ",iter,"\n");
          new.tuo <- tuo_j; new.ci <- c_i;
        }else{
          if(trace) cat("Model parameters eta updated","\n")
          new.tuo <- q$par[1:n.f]; x2 <- q$par[-(1:n.f)];
          new.ci <- x2[1:n.s]
          if(q$convergence != 0) { if(trace) cat("Optimization of eta did not converge on iteration step ", iter,"\n") }
        }
      }

      q1 <- list(value = beta_f(c(new.factor_coefs_j),x=new.factor_coefs_0,f=new.factor_scores,s=new.sigma,pi=new.pi,g=new.gamma))
      b.new.logfunc <- q1$value
      b.cur.logfunc <- b.new.logfunc

      q2 <- list(value = f_f(c(new.factor_scores),x=new.factor_coefs_0, b=new.factor_coefs_j,s=new.sigma,pi=new.pi,g=new.gamma))
      new.f.cur.logfunc <- q2$value
      f.cur.logfunc <- new.f.cur.logfunc

      q3 <- list(value = b0_f(c(new.factor_coefs_0), f=new.factor_scores,b=new.factor_coefs_j,s=new.sigma,pi=new.pi,g=new.gamma))
      new.b0.cur.logfunc <- q3$value
      b0.cur.logfunc <- new.b0.cur.logfunc

      q4 <- list(value = ga_f(c(new.gamma), f=new.factor_scores,b=new.factor_coefs_j,s=new.sigma,pi=new.pi,x=new.factor_coefs_0))
      new.ga.cur.logfunc <- q4$value
      ga.cur.logfunc <- new.ga.cur.logfunc

      q5 <- list(value = eta_f(c(new.tuo,new.ci),g=new.gamma, f=new.factor_scores,b=new.factor_coefs_j,s=new.sigma,pi=new.pi,x=new.factor_coefs_0))
      new.e.cur.logfunc <- q5$value
      e.cur.logfunc <- new.e.cur.logfunc

      ## Take values of VLB to define stopping rule
      q <- list(value = VLB(c(new.factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma,pi=new.pi,e=c(new.tuo,new.ci),g=new.gamma))
      new.VLB <- q$value
      diff=abs(new.VLB-cur.VLB)
      ratio <- abs(new.VLB/cur.VLB);
      if(trace) cat("New VLB:", new.VLB,"cur VLB:", cur.VLB, "Ratio of VLB", ratio, ". Difference in VLB:",diff,"\n")
      cur.VLB <- new.VLB

      factor_coefs_0 <-new.factor_coefs_0
      factor_coefs_j <- new.factor_coefs_j
      pi <- new.pi
      tuo_j <- new.tuo;
      c_i <- new.ci;
      factor_scores <- new.factor_scores
      sigma <- new.sigma
      gamma <- new.gamma
      iter <- iter + 1
    }

    ll <- matrix(factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(gamma,n.s,n.f,byrow=TRUE)*Y+factor_scores %*% t(factor_coefs_j)
    exp.mat <- exp(ll+0.5*(sigma) %*% t(factor_coefs_j^2))
    sum <- exp.mat/(rowSums(exp.mat))
    exp_z <- (1-new.pi)*exp(ll+0.5*(sigma) %*% t(factor_coefs_j^2))
    sum_z <- (1-new.pi)*exp.mat/(rowSums((1-new.pi)*exp.mat))


    if(iter > 99){
      print("ZILNMVA Not converging!")
    }

    ## print the output
    out.list$VLB <- cur.VLB
    out.list$iter=iter-1
    out.list$lvs$pi <- pi
    out.list$lvs$factor_scores <- factor_scores
    out.list$params$tuo <- tuo_j
    out.list$params$c <- c_i
    out.list$params$factor_coefs_j <- factor_coefs_j
    out.list$params$factor_coefs_0 <- factor_coefs_0
    out.list$params$gamma <- gamma
    out.list$lvs$sigma <- sigma
    out.list$Q <- sum
    out.list$mu <- exp.mat
    out.list$muz <- exp_z
    out.list$Qz <- sum_z

    return(out.list)
  }

  if(rank==FALSE){
    re <- tryCatch({ZILNMVA(X,V,n.factors,trace,maxit)},error=function(e){NaN})
    #re <- ZILNMVA(X,V,n.factors,trace,maxit)

  }else{
    if (parallel){
      #cl <- parallel::makeCluster(detectCores(logical = FALSE))
      cl <- parallel::makeCluster(getOption("cl.cores", 4))
      doParallel::registerDoParallel(cl)
    }
    out.list <- list()
    p <- ncol(X)
    n <- nrow(X)
    r <- 5
    fold <- 5
    beta <- list()
    beta0 <- list()
    f <- list()
    Q <- list()
    pi <- list()
    mu <- list()
    gamma <- list()
    sigma <- list()
    Qz <- list()
    muz <- list()
    tuo <- list()
    c <- list()
    iter <- rep(NaN,r)
    L <- rep(NaN,r)
    G_w <- rep(NaN,r);bic <- rep(NaN,r);

    if (parallel){
      Mres <- foreach::foreach(w=1:r) %dopar% {
        re <- tryCatch({ZILNMVA(X,V,n.factors=w,trace,maxit)},error=function(e){NaN})
        re
      }
      for(w in 1:r){
        if(!is.na(Mres[[w]]$VLB)){
          L[w] <- Mres[[w]]$VLB
          iter[w] <- Mres[[w]]$iter
          beta[[w]] <- Mres[[w]]$params$factor_coefs_j
          beta0[[w]] <- Mres[[w]]$params$factor_coefs_0
          gamma[[w]] <- Mres[[w]]$params$gamma
          f[[w]] <- Mres[[w]]$lvs$factor_scores
          sigma[[w]] <- Mres[[w]]$lvs$sigma
          pi[[w]] <- Mres[[w]]$lvs$pi
          mu[[w]] <- Mres[[w]]$mu
          Qz[[w]] <- Mres[[w]]$Qz
          Q[[w]] <- Mres[[w]]$Q
          muz[[w]] <- Mres[[w]]$muz
          tuo[[w]] <- Mres[[w]]$params$tuo
          c[[w]] <- Mres[[w]]$params$c
          G_w[w] <- w*p-w^2+2*w*n
          bic[w] <- -2*L[w]+(log(n)+log(p))*G_w[w]
        }else{
          L[w] <- NaN
          iter[w] <- NaN
          beta[[w]] <- NaN
          beta0[[w]] <- NaN
          tuo[[w]] <- NaN
          c[[w]] <- NaN
          gamma[[w]] <- NaN
          sigma[[w]] <- NaN
          f[[w]] <- NaN
          Q[[w]] <- NaN
          Qz[[w]] <- NaN
          pi[[w]] <- NaN
          G_w[w] <- NaN
          bic[w] <- NaN
          mu[[w]] <- NaN
          muz[[w]] <- NaN
        }
      }
    }else{
      for(w in 1:r){
        re <- tryCatch({ZILNMVA(X,V,n.factors=w,trace,maxit)},error=function(e){NaN})
        if(!is.na(re$VLB)){
          L[w] <- re$VLB
          iter[w] <- re$iter
          beta[[w]] <- re$params$factor_coefs_j
          beta0[[w]] <- re$params$factor_coefs_0
          gamma[[w]] <- re$params$gamma
          sigma[[w]] <- re$lvs$sigma
          f[[w]] <- re$lvs$factor_scores
          Q[[w]] <- re$Q
          mu[[w]] <- re$mu
          Qz[[w]] <- re$Qz
          muz[[w]] <- re$muz
          tuo[[w]] <- re$params$tuo
          c[[w]] <- re$params$c
          pi[[w]] <- re$lvs$pi
          G_w[w] <- w*p-w^2+2*w*n
          bic[w] <- -2*L[w]+(log(n)+log(p))*G_w[w]
        }else{
          L[w] <- NaN
          iter[w] <- NaN
          beta[[w]] <- NaN
          beta0[[w]] <- NaN
          tuo[[w]] <- NaN
          c[[w]] <- NaN
          gamma[[w]] <- NaN
          sigma[[w]] <- NaN
          f[[w]] <- NaN
          Q[[w]] <- NaN
          Qz[[w]] <- NaN
          pi[[w]] <- NaN
          G_w[w] <- NaN
          bic[w] <- NaN
          mu[[w]] <- NaN
          muz[[w]] <- NaN
        }
      }
    }
    out.list$bic <- which.min(bic)
    out.list$VLB <- L[[(which.min(bic))]]
    out.list$iter <- iter[[(which.min(bic))]]
    out.list$lvs$pi <- pi[[(which.min(bic))]]
    out.list$lvs$factor_scores <- f[[(which.min(bic))]]
    out.list$lvs$factor_scores2 <- f[[2]]
    out.list$lvs$sigma <- sigma[[(which.min(bic))]]
    out.list$params$factor_coefs_j <- beta[[(which.min(bic))]]
    out.list$params$factor_coefs_0 <- beta0[[(which.min(bic))]]
    out.list$params$gamma <- gamma[[(which.min(bic))]]
    out.list$params$tuo <- tuo[[(which.min(bic))]]
    out.list$params$c <- c[[(which.min(bic))]]
    out.list$Q <- Q[[(which.min(bic))]]
    out.list$mu <- mu[[(which.min(bic))]]
    out.list$muz <- muz[[(which.min(bic))]]
    out.list$Qz <- Qz[[(which.min(bic))]]

    if (parallel){
      parallel::stopCluster(cl = cl)
    }
    return(out.list)

  }
}
