require(foreach)
require(doParallel)
source("proximal_gradient.R")
source("basic_functions.R")

#----------------------------------------------------------------------------------------
#  Tuning parameter for Proximal Gradient (nRep = 1)
#  Input:
#                 W ------ n x p count data matrix (row/column is sample/variable)
#            alphaX ------ lower bound of the constraint
#      lambda_range ------ the range of lambda
#          fracTest ------ the fraction of test set
#  Output:
#        lambda_hat ------ selected parameter lambda
#        alphaX_hat ------ selected parameter alphaX
#             X_hat ------ estimated composition
#-----------------------------------------------------------------------------------------
tuneProxGradient <- function(W, alphaX, lambda_range, fracTest = 1/5){
  n <- nrow(W)
  p <- ncol(W)
  n_lambda <- length(lambda_range)
  X_test <- W/(rowSums(W)%*%matrix(1,1,p))
  avg_x_test <- colSums(X_test) / (colSums(X_test > 0) + 1e-10)
  
  # save variables for debugging
  error_grid <- matrix(0, n_lambda, 1)
  
  for (i in 1:n_lambda){
    lambda <- lambda_range[i]
    X_pg_temp_res <- proxGradient(W = W, lambda = lambda, alphaX = alphaX, 
                                  betaX = p, L = 1e-4, gamma = 5, iterMax = 1e+4,
                                  epsilon = 1e-5)
    X_hat_temp <- X_pg_temp_res$XHat
    rowSumComp <- rowSums(X_hat_temp)
    if (any(abs(rowSumComp-1)>1e-10)){
      stop("Bug! Row sum is not 1!")
      }
    X_hat_temp[X_test > 0] <- 0
    avg_x_hat <- colSums(X_hat_temp) / (colSums(X_hat_temp > 0) + 1e-10)
    error_grid[i] <- mean(abs(avg_x_test - avg_x_hat))
    }
  
  min_row_ind <- tail(which(error_grid == min(error_grid), arr.ind = TRUE), 1)[1]
  lambda_hat <- lambda_range[min_row_ind]
  X_hat_res <- proxGradient(W = W, lambda = lambda_hat, alphaX = alphaX, 
                            betaX = p, L = 1e-4, gamma = 5, iterMax = 1e+4,
                            epsilon = 1e-5)
  
  return(list(lambda_hat = lambda_hat,X_hat = X_hat_res$XHat,S=X_hat_res$S,d = X_hat_res$d,d_vec=X_hat_res$d_vec))
}

#----------------------------------------------------------------------------------------
#  Proximal Gradient automatically using cross validation
#  Input:
#                 W ------ n x p count data matrix (row/column is sample/variable)
#            n_grid ------ number of grid for each parameter (lambda and alphaX)
#  Output:
#        lambda_hat ------ selected parameter lambda
#             X_hat ------ estimated composition
#-----------------------------------------------------------------------------------------
autoTuneProxGradient <- function(W, n_grid = 5){
  n <- nrow(W)
  p <- ncol(W)
  X <- W/(rowSums(W)%*%matrix(1,1,p))
  
  # set tuning parameter range
  alphaX <- min(min(X[X>0])*(1e-2)*p, 1)
  lambda_min <- 0.2
  lambda_max <- 0.4
  
  index <- 1:n_grid
  lambda_range <- lambda_min * (lambda_max/lambda_min)^((index-1)/(n_grid-1))
  
  # tuning parameter result
  tune <- tuneProxGradient(W, alphaX, lambda_range)
  return(list(lambda_hat = tune$lambda_hat, S=tune$S,X_hat = tune$X_hat,d=tune$d,d_vec=tune$d_vec))
}

