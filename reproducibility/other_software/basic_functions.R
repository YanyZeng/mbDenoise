#-----------------------------------------------------------------------------------------
#  percentage mean square error between matrix X (basis) and Y
#-----------------------------------------------------------------------------------------
diff_percentage_mse <- function(X, Y){
  non_zero_index <- which(X != 0, arr.ind = T)
  val <- sum(((X[non_zero_index]-Y[non_zero_index])^2)/(X[non_zero_index]^2))/nrow(non_zero_index)
  return(val = val)
}

#-----------------------------------------------------------------------------------------
#  percentage mean absolute error between matrix X (basis) and Y
#-----------------------------------------------------------------------------------------
diff_percentage_mae <- function(X, Y){
  non_zero_index <- which(X != 0, arr.ind = T)
  val <- sum((abs(X[non_zero_index]-Y[non_zero_index]))/abs(X[non_zero_index]))/nrow(non_zero_index)
  return(val = val)
}

#----------------------------------------------------------------------------------------
#  Parameter selection for lambda 
#  Input:
#                 W ------ n x p count data matrix (row/column is sample/variable)
#            alphaX ------ selected parameter alphaX
#             betaR ------ selected parameter betaR
#  Output:
#        lambda_hat ------ selected parameter lambda
#-----------------------------------------------------------------------------------------
tune_para <- function(w, alphaX, betaR){
  n <- nrow(w)
  p <- ncol(w)
  N <- sum(w)
  lambda_hat <- (4/(3*alphaX) + 4*sqrt(1/(9*alphaX^2)+betaR/alphaX)) * sqrt(p*max(n,p)*log(n+p)/(n*N))
  return(lambda_hat)
}

#----------------------------------------------------------------------------------------
#  KL-divergence
#  Input:
#                 X ------ n x p composition matrix (row/column is sample/variable)
#                 Y ------ n x p composition matrix (row/column is sample/variable)
#  Output:
#               val ------ KL divergence between X and Y
#-----------------------------------------------------------------------------------------
KL <- function(X, Y){
  ind <- which(X != 0, arr.ind = T)
  val <- sum((log(X[ind])-log(Y[ind]))*X[ind])
  return(val = val)
}

#----------------------------------------------------------------------------------------
#  split the data to training set and test set
#  Input:
#                 W ------ n x p count matrix (row/column is sample/variable)
#          fracTest ------ the fraction of test set
#  Output:
#           W_train ------ the training set (n x p count matrix), the elements in training index of W_train are non-zero
#    index_test_row ------ the row index of test set
#    index_test_col ------ the column index of test set
#-----------------------------------------------------------------------------------------
split_train_test <- function(W, fracTest = 1/4){
  
  non_zero_index <- which(W != 0, arr.ind = T)
  n_nonzero <- nrow(non_zero_index)
  test_index <- sample(n_nonzero, floor(n_nonzero * fracTest))
  
  index_test_row <- matrix(0, 1, length(test_index))
  index_test_col <- matrix(0, 1, length(test_index))
  
  while (TRUE){
    W_train <- W
    for (l in 1:length(test_index)){
      row_nonzero <- non_zero_index[test_index[l], 1]
      col_nonzero <- non_zero_index[test_index[l], 2]
      index_test_row[l] <- row_nonzero
      index_test_col[l] <- col_nonzero
      W_train[row_nonzero, col_nonzero] <- 0
    }
    
    # check whether there exists whole zero row
    W_b <- (W_train > 0)
    rowZero <- rowSums(W_b) 
    colZero <- colSums(W_b)
    
    if ((all(rowZero >0)) & (all(colZero > 0))){
      break
    }
  }
  return(list(W_train = W_train, index_test_row = index_test_row, index_test_col = index_test_col))
}

#----------------------------------------------------------------------------------------
#  Diversity index
#  Input:
#                 X ------ n x p composition matrix (row/column is sample/variable)
#  Output:
#                sh ------ shannon index (n x 1 vector)
#                sp ------ simpson index (n x 1 vector)
#                bc ------ bray-curtis index (n(n-1)/2 x 1 vector)
#-----------------------------------------------------------------------------------------
diversity <- function(X){
  # shannon index
  sh_ele <- -X * log(X)
  sh_ele[is.nan(sh_ele)] <- 0
  sh <- rowSums(sh_ele)
  
  # simpson index
  sp <- rowSums(X * X)
  
  # bray-curtis index
  n <- nrow(X)
  bc <- matrix(vegdist(X), n*(n-1)/2, 1)
  return(list(sh = sh, sp = sp, bc = bc))
}
