require(vegan) # this package is used to calculate bray-curtis index
require(corpcor)
# source("basic_functions.R")

#----------------------------------------------------------------------------------------
#  zero replacement estimator
#  Input:
#                 W ------ n x p count matrix (row/column is sample/variable)
#             alpha ------ zero replacement parameter
#  Output:
#             X_hat ------ estimated zero-replacement composition
#-----------------------------------------------------------------------------------------
zr <- function(W, alpha = 0.5){
  W[W < alpha] <- alpha
  X_hat <- W/(rowSums(W)%*%matrix(1,1,ncol(W)))
  return(X_hat = X_hat)  
}

#----------------------------------------------------------------------------------------
#  SVT estimator
#  Input:
#                 W ------ n x p count matrix (row/column is sample/variable)
#                 r ------ selected rank
#             alpha ------ zero replacement parameter
#  Output:
#             X_hat ------ estimated SVT composition
#-----------------------------------------------------------------------------------------
svt <- function(W, r, alpha = 0.5){
  S <- fast.svd(W)
  u_r <- S$u[,1:r]
  v_r <- S$v[,1:r]
  d_r <- S$d[1:r]
  W_hat <- u_r %*% diag(d_r) %*% t(v_r)
  X_hat <- zr(W_hat, alpha)
  return(X_hat = X_hat)
}

#----------------------------------------------------------------------------------------
#  Tuning parameter for zero-replacement
#  Input:
#                 W ------ n x p count data matrix (row/column is sample/variable)
#       alpha_range ------ range of parameters
#  Output:
#         alpha_hat ------ selected parameter alpha
#             X_hat ------ estimated composition
#-----------------------------------------------------------------------------------------
tuneZr <- function(W, alpha_range){
  
  n <- nrow(W)
  p <- ncol(W)
  n_alpha <- length(alpha_range)  
  X_test <- W/(rowSums(W)%*%matrix(1,1,p))
  avg_x_test <- colSums(X_test) / (colSums(X_test > 0) + 1e-10)
  
  # save variables for debugging
  error_grid <- matrix(0, n_alpha, 1)
  
  for (i in 1:n_alpha){
    alpha <- alpha_range[i]
    X_hat_temp <- zr(W, alpha)
    rowSumComp <- rowSums(X_hat_temp)
    if (any(abs(rowSumComp-1)>1e-10)){
      stop("Bug! Row sum is not 1!")
    }
    X_hat_temp[X_test > 0] <- 0
    avg_x_hat <- colSums(X_hat_temp) / (colSums(X_hat_temp > 0) + 1e-10)
    error_grid[i] <- mean(abs(avg_x_test - avg_x_hat))
  }
  
  min_row_ind <- tail(which(error_grid == min(error_grid), arr.ind = TRUE), 1)[1]
  alpha_hat <- alpha_range[min_row_ind]
  X_hat <-  zr(W, alpha_hat)
  return(list(X_hat = X_hat, alpha_hat = alpha_hat))
}

#----------------------------------------------------------------------------------------
#  Zero-replacement automatically using cross validation
#  Input:
#                 W ------ n x p count data matrix (row/column is sample/variable)
#  Output:
#         alpha_hat ------ selected parameter alphaX
#             X_hat ------ estimated composition
#-----------------------------------------------------------------------------------------
autoTuneZr <- function(W){
  
  # set tuning parameter range
  p <- ncol(W)
  alpha_min <- min(1/p, 1e-3)
  alpha_max <- max(min(W[W>0]), 100)
  n_grid <- 1000
  index <- 1:n_grid
  alpha_range <- alpha_min * (alpha_max/alpha_min)^((index-1)/(n_grid-1))
  
  # results
  tune <- tuneZr(W, alpha_range)
  return(list(alpha_hat = tune$alpha_hat, X_hat = tune$X_hat))
}

#----------------------------------------------------------------------------------------
#  Tuning parameter for SVT
#  Input:
#                 W ------ n x p count data matrix (row/column is sample/variable)
#                 r ------ selected rank r
#       alpha_range ------ range of parameters
#  Output:
#        alphaX_hat ------ selected parameter alphaX
#             X_hat ------ estimated composition
#-----------------------------------------------------------------------------------------
tuneSvt <- function(W, r, alpha_range){
  
  n <- nrow(W)
  p <- ncol(W)
  n_alpha <- length(alpha_range)  
  X_test <- W/(rowSums(W)%*%matrix(1,1,p))
  avg_x_test <- colSums(X_test) / (colSums(X_test > 0) + 1e-10)
  
  # save variables for debugging
  error_grid <- matrix(0, n_alpha, 1)
  
  for (i in 1:n_alpha){
    alpha <- alpha_range[i]
    X_hat_temp <- svt(W, r, alpha)
    rowSumComp <- rowSums(X_hat_temp)
    if (any(abs(rowSumComp-1)>1e-10)){
      stop("Bug! Row sum is not 1!")
    }
    X_hat_temp[X_test > 0] <- 0
    avg_x_hat <- colSums(X_hat_temp) / (colSums(X_hat_temp > 0) + 1e-10)
    error_grid[i] <- mean(abs(avg_x_test - avg_x_hat))
  }
  
  min_row_ind <- tail(which(error_grid == min(error_grid), arr.ind = TRUE), 1)[1]
  alpha_hat <- alpha_range[min_row_ind]
  X_hat <- svt(W = W, r = r, alpha = alpha_hat)
  
  return(list(X_hat = X_hat, alpha_hat = alpha_hat))
}

#----------------------------------------------------------------------------------------
#  SVT automatically using cross validation
#  Input:
#                 W ------ n x p count data matrix (row/column is sample/variable)
#                 r ------ selected rank r
#  Output:
#         alpha_hat ------ selected parameter alphaX
#             X_hat ------ estimated composition
#-----------------------------------------------------------------------------------------
autoTuneSvt <- function(W, r){
  
  # set tuning parameter range
  p <- ncol(W)
  alpha_min <- min(1/p, 1e-3)
  alpha_max <- max(min(W[W>0]), 10)
  n_grid <- 1000
  index <- 1:n_grid
  alpha_range <- alpha_min * (alpha_max/alpha_min)^((index-1)/(n_grid-1))
  
  # results
  tune <- tuneSvt(W, r, alpha_range)
  return(list(alpha_hat = tune$alpha_hat, X_hat = tune$X_hat))
}