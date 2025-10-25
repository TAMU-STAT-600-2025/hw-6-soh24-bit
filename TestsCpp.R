
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)

# Libraries for bench and load in example data set
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

data1 <- mtcars
X1 <- as.matrix(data1[,3:7])
Y1 <- as.matrix(data1[,1])
standardized_data1 <- standardizeXY(X1, Y1)
Xtilde1 <- standardized_data1$Xtilde # 2x2 Identity matrix
Ytilde1 <- standardized_data1$Ytilde
beta1 <- rep(1.5, 5)
lambda1 <- 0.5

mark(
  lasso(Xtilde1, Ytilde1, beta1, lambda1),
  lasso_c(Xtilde1, Ytilde1, beta1, lambda1),
  check = TRUE
)

data2 <- iris
X2 <- as.matrix(data2[,2:4])
Y2 <- as.matrix(data2[,1])
standardized_data2 <- standardizeXY(X2, Y2)
Xtilde2 <- standardized_data2$Xtilde # 2x2 Identity matrix
Ytilde2 <- standardized_data2$Ytilde
beta2 <- rep(2.5, 3)
lambda2 <- 0.5

mark(
  lasso(Xtilde2, Ytilde2, beta2, lambda2),
  lasso_c(Xtilde2, Ytilde2, beta2, lambda2),
  check = TRUE
)

# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

beta_start1 <- rep(-1.5, 5)

# cpp version only includes beta so set check = FALSE. However, I verified the beta values are equal
mark(
  fitLASSOstandardized(Xtilde1, Ytilde1, lambda1, beta_start1),
  fitLASSOstandardized_c(Xtilde1, Ytilde1, lambda1, beta_start1),
  check = FALSE
)

beta_start2 <- c(0.1, -0.2, 0.5)

# cpp version only includes beta so set check = FALSE. However, I verified the beta values are equal
mark(
  fitLASSOstandardized(Xtilde2, Ytilde2, lambda2, beta_start2),
  fitLASSOstandardized_c(Xtilde2, Ytilde2, lambda2, beta_start2),
  check = FALSE
)

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

microbenchmark::microbenchmark(
  fitLASSOstandardized(Xtilde1, Ytilde1, lambda1, beta_start1),
  fitLASSOstandardized_c(Xtilde1, Ytilde1, lambda1, beta_start1),
  times = 10
)

microbenchmark::microbenchmark(
  fitLASSOstandardized(Xtilde2, Ytilde2, lambda2, beta_start2),
  fitLASSOstandardized_c(Xtilde2, Ytilde2, lambda2, beta_start2),
  times = 10
)

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

# We are allowed to assume lambda_seq is already sorted

lambda

mark(
  fitLASSOstandardized(Xtilde1, Ytilde1, lambda1, beta_start1),
  fitLASSOstandardized_c(Xtilde1, Ytilde1, lambda1, beta_start1),
  check = FALSE
)


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
