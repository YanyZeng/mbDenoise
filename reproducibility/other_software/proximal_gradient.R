require(corpcor)

#----------------------------------------------------------------------------------------
#  Proximal Gradient
#  minimize L(X) + lambda*||X||_nuc
#  s.t.   X*1 = 1, betaX/p >= X >= alphaX/p
#  Input:
#           W ------ n x p count data matrix (row/column is sample/variable)
#      lambda ------ regularized parameter
#      alphaX ------ lower bound of the constraint
#       betaX ------ upper bound of the constraint
#           r ------ friction parameter
#           L ------ the step size
#       gamma ------ the parameter used in line search
#     iterMax ------ the maximum number of iteration
#  Output:
#        XHat ------ estimated composition
#      n_iter ------ the number of iteration
#-----------------------------------------------------------------------------------------
proxGradient <- function(W, lambda, alphaX, betaX, r = 10, L = 1e-4, 
                         gamma = 3, iterMax = 1e+4, epsilon = 1e-6){
  
  p <- ncol(W)
  X <- W/(rowSums(W)%*%matrix(1,1,p))
  Y <- X
  k <- 1 # the number of iteration
  findStepSize <- FALSE
  diff <- matrix(0, iterMax, 1)
  while (k <= iterMax){
    # proximal gradient step
    while (!findStepSize){
      # singular value thresholding
      # S <- fast.svd(Y - 1/L * gradFun(W, Y))
      S <- svd(Y - 1/L * gradFun(W, Y))      
      d <- pmax(0, S$d - lambda/L)
      # project the updating variable onto the simplex constraint
      XTmp <- projSimplexMat(S$u %*% diag(d) %*% t(S$v), alphaX, betaX)
      # find a step size
      if (criteriaFun(W, XTmp, Y, L) > 0){
        L <- gamma * L
      }else{
        findStepSize <- TRUE
      }
    }
    findStepSize <- FALSE
    # exit criteria
    diff[k] <- abs(criteriaFun(W, XTmp, X, L))
    if (diff[k] < epsilon){
      return(list(XHat = XTmp, n_iter = k))
    }
    Y <- XTmp + (k - 1) / (k + r - 1) * (XTmp - X)
    Y <- projSimplexMat(Y, alphaX, betaX)
    X <- XTmp
    k <- k + 1
  }
  return(list(XHat = X, n_iter = k,d=sum(d!=0),d_vec=d,S=S))
}

#----------------------------------------------------------------------------------------
#  Objection Function
#  L(X) + lambda * nuclear norm of X
#  Input:
#           W ------ n x p count data matrix (row/column is sample/variable)
#           X ------ n x p estimated composition matrix (row/column is sample/variable)
#      lambda ------ regularization parameter
#  Output:
#         val ------ function value
#-----------------------------------------------------------------------------------------
objFun <- function(W, X, lambda){
  # negative likelihood
  ind <- which(W != 0, arr.ind = T)
  LX <- -sum((log(X[ind])) * W[ind])/sum(W)
  # nuclear norm penalty
  S <- fast.svd(X)
  nuclearNorm <- sum(S$d)
  # objective function
  val <- LX + lambda*nuclearNorm
  return(val = val)
}

#----------------------------------------------------------------------------------------
#  Criteria Function
#  L(X1) - L(X2) - <X1-X2, gradFun(X2)> - L*||X1-X2||_F^2/2
#  Input:
#           W ------ n x p count data matrix (row/column is sample/variable)
#          X1 ------ n x p composition data matrix (row/column is sample/variable)
#          X2 ------ n x p composition data matrix (row/column is sample/variable)
#           L ------ the step size
#  Output:
#         val ------ function value
#-----------------------------------------------------------------------------------------
criteriaFun <- function(W, X1, X2, L){
  ind <- which(W != 0, arr.ind = T)
  val <- sum((log(X2[ind]) - log(X1[ind])) * W[ind]) / sum(W) - sum((X1 - X2) * gradFun(W, X2)) - L/2 * (norm(X1 - X2, type = "F"))^2
  return(val = val)
}

#----------------------------------------------------------------------------------------
#  Gradient Function of L
#  Input:
#           W ------ n x p count data matrix (row/column is sample/variable)
#           X ------ n x p composition data matrix (row/column is sample/variable)
#  Output:
#         val ------ function value
#-----------------------------------------------------------------------------------------
gradFun <- function(W, X){
  val <- matrix(0, dim(W)[1], dim(W)[2])
  val[which(W > 1e-10, arr.ind = T)] = -1 / sum(W) * W[which(W > 1e-10, arr.ind = T)] / X[which(W > 1e-10, arr.ind = T)]
  return(val = val)
}

#----------------------------------------------------------------------------------------
#  Project a matrix X into the constrained simplex space XHat:
#  XHat * 1 = 1 and alpha / p <= XHat <= beta / p
#  Input:
#           X ------ n x p data matrix (row/column is sample/variable)
#       alpha ------ the lower bound (alpha <= 1 <= beta)
#        beta ------ the upper bound (beta >= 1 >= alpha)
#  Output:
#        XHat ------ the projected composition 
#-----------------------------------------------------------------------------------------
projSimplexMat <- function(X, alpha, beta){
  XHat <- X
  p <- ncol(X)
  beta <- max(beta, 1)
  if (abs(beta/p - 1) < 1e-10){
    XHat <- projSimplexMatOnlyAlpha(X = X, alpha = alpha)
    return (XHat = XHat)
  }
  alpha <- min(alpha, 1)
  V <- cbind(X-alpha/p, X-beta/p)
  V <- t(apply(V, 1, sort))
  for (i in 1:nrow(X)){
    XHat[i,] <- projSimplexVec(x = X[i,], alpha = alpha, beta = beta, v = V[i,])
  }
  return (XHat = XHat)
}

#----------------------------------------------------------------------------------------
#  Project a matrix X into the constrained simplex space XHat:
#  XHat * 1 = 1 and alpha / p <= XHat <= 1
#  Input:
#           X ------ n x p data matrix (row/column is sample/variable)
#       alpha ------ the lower bound (alpha <= 1 <= beta)
#  Output:
#        XHat ------ the projected composition 
#-----------------------------------------------------------------------------------------
projSimplexMatOnlyAlpha <- function(X, alpha){
  XHat <- X
  p <- ncol(X)
  V <- t(apply(X, 1, sort, decreasing = TRUE))
  for (i in 1:nrow(X)){
    XHat[i,] <- projSimplexVecOnlyAlpha(x = X[i,], alpha = alpha, v = V[i,])
  }
  return (XHat = XHat)
}

#----------------------------------------------------------------------------------------
#  Project a vector x into the constrained simplex space xhat:
#  xhat * 1 = 1 and alpha / p <= xhat <= beta / p
#  Input:
#           x ------ 1 x p data vector 
#       alpha ------ the lower bound (alpha <= 1 <= beta)
#        beta ------ the upper bound (beta >= 1 >= alpha)
#           v ------ an increasing 1 x 2p vector
#  Output:
#        xhat ------ the projected composition 
#-----------------------------------------------------------------------------------------
projSimplexVec <- function(x, alpha, beta, v){
  p <- length(x)
  x <- matrix(x, 1, p)
  alpha <- min(alpha, 1)
  beta <- max(beta, 1)
  d <- sum(pmax(pmin(x - v[1], beta / p), alpha / p)) - 1
  j <- 2
  while (j < length(v)){
    dnew <- sum(pmax(pmin(x - v[j], beta / p), alpha / p)) - 1
    if ((dnew <= 0) & (d >= 0)){
      break
    }else{
      d <- dnew
      j <- j + 1
    }
  }
  jstar <- j - 1
  xtmp <- x - v[jstar]
  lambda <- v[jstar] + d / length(which((xtmp <= beta /p) & (xtmp > alpha / p)))
  xhat <- pmax(pmin(x - lambda, beta / p), alpha / p)
  return (xhat = xhat)
}

#----------------------------------------------------------------------------------------
#  Project a vector x into the constrained simplex space xhat:
#  xhat * 1 = 1 and alpha / p <= xhat <= 1
#  Input:
#           x ------ 1 x p data vector 
#       alpha ------ the lower bound (alpha <= 1 <= beta)
#           v ------ the decreasing sorted x (1 x p vector)
#  Output:
#        xhat ------ the projected composition 
#-----------------------------------------------------------------------------------------
projSimplexVecOnlyAlpha <- function(x, alpha, v){
  p <- length(x)
  x <- matrix(x, 1, p)
  alpha <- min(alpha, 1)
  j <- p
  while (j > 1){
    if (v[j] + (1-alpha-sum(v[1:j]))/j > 0){
      break
    }else{
      j <- j - 1
    }
  }
  theta <- (1-alpha-sum(v[1:j]))/j + alpha/p
  xhat <- pmax(x + theta, alpha/p)
  return (xhat = xhat)
}

