library(igraph)
library(kernlab)
library(blockmodels)
library(graphon)
library(fossil)
library(mixer)
library(rmutil)
library(optimx)
library(ggplot2)
require(gridExtra)
library(randnet)
library(RMTstat)

llfunc_diag <- function(par, K, X, n){
  alpha <- par[1]
  beta <- par[2]
  if(alpha < 0 | beta < 0) 
    return(Inf)
  else{
    ll <- -K*lbeta(alpha, beta)
    for(i in 1:K)
      ll <- ll + lbeta(alpha + X[i,i], beta + n[i,i] - X[i,i])
    return(-ll)
  }
}

llfunc_offdiag <- function(par, K, X, n){
  alpha <- par[1]
  beta <- par[2]
  if(alpha < 0 | beta < 0) 
    return(Inf)
  else{
    ll <- -K*(K-1)*lbeta(alpha, beta)/2
    for(i in 1:K)
      for(j in 1:i)
        if(j != i)
          ll <- ll + lbeta(alpha + X[i,j], beta + n[i,j] - X[i,j])
    return(-ll)
  }
}

SBM_gof <- function(Theta, Y){
  B = Theta
  A = Y
  n = dim(A)[1]
  A_tilde = matrix(0, n, n) 
  for(i in 1:n){
    for(j in 1:n){
      if(i != j)
        A_tilde[i, j] = (A[i, j] - B[i, j])/sqrt((n - 1) * B[i, j] * (1 - B[i, j]))
    }
  }
  lambdas = eigen(A_tilde)$values
  lambda_1 = Re(lambdas[1])
  lambda_n = Re(tail(lambdas, 1))
  
  T_k0 = max(n^(2/3)*(lambda_1 - 2), n^(2/3)*(-lambda_n - 2))

  M = 100  
  lambda_m_1 <- rep(0, M)
  lambda_m_n <- rep(0, M)
  for(m in 1:M){
    A_m = matrix(0, n, n)
    for(i in 1:n){
      for(j in 1:n){
        A_m[i, j] = rbinom(1, 1, B[i, j])
      }
    }
    A_m_tilde = matrix(0, n, n)
    for(i in 1:n){
      for(j in 1:n){
        if(i != j)
          A_m_tilde[i, j] = (A_m[i, j] - B[i, j])/sqrt((n - 1) * B[i, j] * (1 - B[i, j]))
      }
    }
    lambdas = eigen(A_m_tilde)$values
    lambda_m_1[m] = Re(lambdas[1])
    lambda_m_n[m] = Re(tail(lambdas, 1))
  }
  mu1 = mean(lambda_m_1)
  s1 = sd(lambda_m_1)
  mun = mean(lambda_m_n)
  sn = sd(lambda_m_n)

  T_boot =  -1.2065335745820 + sqrt(1.607781034581) * max((lambda_1 - mu1)/s1, (lambda_n - mun)/sn)
  return(c(T_boot, T_k0))
}
