# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - mean(Y)
  # [ToDo] Center and scale X
  Xmeans <- colMeans(X)
  Xcentered <- X - matrix(colMeans(X), nrow(X), ncol(X), byrow = TRUE)
  Xnorms <- colSums(Xcentered ^ 2) / nrow(X)
  Xtilde <- Xcentered %*% diag(1 / sqrt(Xnorms))
  
  weights <- sqrt(Xnorms)
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  return(sign(a) * max(abs(a) - lambda, 0))
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  # Get the row of X
  n <- nrow(Xtilde)
  # Calculate objective
  fobj <- (2*n)^(-1) * sum((Ytilde - Xtilde %*% beta)^2) + lambda * sum(abs(beta))
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  if (n != length(Ytilde)) {
    stop("n has to be the same between Xtilde and Ytilde")
  }
  #[ToDo]  Check that lambda is non-negative
  if (lambda < 0) {
    stop("Lambda must be non-negative")
  }
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  if (is.null(beta_start)) {
    beta_start <- rep(0, p)
  } else if (length(beta_start) != p) {
    stop("p should match between Xtilde and beta")
  }
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  
  # Initilize beta, fobj
  beta <- beta_start
  fcurrent <- lasso(Xtilde, Ytilde, beta, lambda)
  
  # For efficiency, initialize Ytilde - Xtilde %*% beta or the residual
  resid <- Ytilde - Xtilde %*% beta
  
  while(TRUE) {
    # Assign previous objective value to compare
    fprevious <- fcurrent
    
    for (i in 1:p) {
      beta_old_i <- beta[i]
      # the term below is equivalent to calculating partial residual, but we avoid calculating
      # Ytilde - Xtilde[, -i] %*% beta[-i] in the loop
      a <- (1 / n) * crossprod(Xtilde[, i], resid) + beta_old_i
      beta[i] <- soft(a, lambda)
      
      # Update the residual with updated beta
      if (beta[i] != beta_old_i) {
        resid <- resid - Xtilde[, i] * (beta[i] - beta_old_i)
      }
    }
    
    fcurrent <- lasso(Xtilde, Ytilde, beta, lambda)
    
    # Test if the difference between fprevious and fcurrent is lower than eps
    if (abs(fprevious - fcurrent) < eps) {
      fmin <- fcurrent
      break
    }
    
  }
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  if (n != length(Ytilde)) {
    stop("n has to be the same between Xtilde and Ytilde")
  }
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  
  # A switch for generating lambda sequence
  switch <- FALSE
  
  if (is.null(lambda_seq)) {
    switch <- TRUE
  } else {
    lambda_seq <- sort(lambda_seq[lambda_seq >= 0], decreasing = TRUE)
    
    if (length(lambda_seq) == 0) {
      message("None of the supplied values of lambdas satisfy the requirement; 
              proceeding as if the values were not supplied")
      switch <- TRUE
    }
  }
  
  if (switch) {
    lambda_max <- max(abs((1 / n) * crossprod(Xtilde, Ytilde)))
    lambda_min <- 0.0001 * lambda_max
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length = n_lambda))
  }
  
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  
  # Initialize
  beta_mat <- matrix(NA, nrow = p, ncol = length(lambda_seq))
  fmin_vec <- numeric(length(lambda_seq))
  
  # Warmstart
  beta_current <- rep(0, p) 
  
  for (i in 1:length(lambda_seq)) {
    
    currentfit <- fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[i], 
                                       beta_start = beta_current, eps = eps)
    
    beta_mat[, i] <- currentfit$beta
    fmin_vec[i] <- currentfit$fmin
    
    beta_current <- currentfit$beta
  }
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  standardized <- standardizeXY(X, Y)
  Xtilde <- standardized$Xtilde
  Ytilde <- standardized$Ytilde
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  fitlasso <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq, 
                                       n_lambda = n_lambda, eps = eps)
  
  # Returns actual lambda sequence used
  lambda_seq <- fitlasso$lambda_seq
  
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  
  # Backscale the X
  beta_mat <- fitlasso$beta_mat
  weights <- standardized$weights
  beta_mat <- beta_mat / weights
  
  meanY <- standardized$Ymean
  meansX <- standardized$Xmeans
  
  # Backscale the intercept
  beta0_vec <- rep(mean(Y), ncol(beta_mat)) - colSums(beta_mat * meansX)
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  # [ToDo] Fit Lasso on original data using fitLASSO
  fit <- fitLASSO(X, Y, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  
  # Output from the whole data
  lambda_seq <- fit$lambda_seq
  beta_mat <- fit$beta_mat
  beta0_vec <- fit$beta0_vec
  
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  n <- nrow(X)
  if (is.null(fold_ids)) {
    fold_seq <- sample(rep(1:k, length.out = n))
  } else{
    fold_seq <- fold_ids
  }
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  
  # Initialize
  cvm = rep(NA, length(lambda_seq)) # want to have CV(lambda)
  cvse = rep(NA, length(lambda_seq)) # want to have SE_CV(lambda)
  cv_folds = matrix(NA, k, length(lambda_seq)) # fold-specific errors
  
  for (fold in 1:k) {
    test_index <- which(fold_seq == fold, arr.ind = TRUE)
    
    Xtrain <- X[-test_index, ]
    Ytrain <- Y[-test_index]
    
    Xtest <- X[test_index, ]
    Ytest <- Y[test_index]
    
    fitKfold <- fitLASSO(Xtrain, Ytrain, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
    Xbeta <- Xtest %*% fitKfold$beta_mat
    
    for (i in 1:ncol(Xbeta)) {
      cv_folds[fold, i] <- sum((Ytest - fitKfold$beta0_vec[i] - Xbeta[, i])^2)
    }
  }
  
  cvm <- colMeans(cv_folds)
  cvse <- apply(cv_folds, 2, \(cv_i) {sd(cv_i) / sqrt(k)})
  
  # [ToDo] Find lambda_min
  lambda_min <- lambda_seq[which.min(cvm)]
  # [ToDo] Find lambda_1SE
  min_cvm <- min(cvm)
  min_se <- cvse[which.min(cvm)]
  error_threshold <- min_cvm + min_se
  lambdas <- lambda_seq[cvm <= error_threshold]
  lambda_1se <- max(lambdas)
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_seq, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

