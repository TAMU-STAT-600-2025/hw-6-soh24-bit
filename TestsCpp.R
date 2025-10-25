
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)

# Libraries for bench
library(bench)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
a1 <- 20
lambda1 <- 2.1

mark(
  soft(a1, lambda1),
  soft_c(a1,lambda1),
  check = TRUE
)

a2 <- -20
lambda2 <- 2.1

mark(
  soft(a2, lambda2),
  soft_c(a2,lambda2),
  check = TRUE
)

# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################

n1 <- 2
p1 <- 2
Xtilde1 <- matrix(c(1, 0, 0, 1), nrow = n1, ncol = p1) # 2x2 Identity matrix
Ytilde1 <- c(3, 5)
beta1 <- c(1, 1)
lambda1 <- 0.5

mark(
  lasso(Xtilde1, Ytilde1, beta1, lambda1),
  lasso_c(Xtilde1, Ytilde1, beta1, lambda1),
  check = TRUE
)

n2 <- 3
p2 <- 2
Xtilde2 <- matrix(c(1.5, 2.1, -0.5, 0.2, 1.0, 3.0), nrow = n2, ncol = p2)
Ytilde2 <- c(5.5, 1.2, 3.0)
beta2 <- c(0.5, -1.2)
lambda2 <- 1.1

mark(
  lasso(Xtilde2, Ytilde2, beta2, lambda2),
  lasso_c(Xtilde2, Ytilde2, beta2, lambda2),
  check = TRUE
)

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)
