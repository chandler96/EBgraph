library(igraph)
library(kernlab)
library(blockmodels)
library(graphon)
library(fossil)
library(mixer)
library(rmutil)
library(optimx)

# likelihood given estimated parameters from VBEM model
llfunc_mixer <- function(alpha_mat, beta_mat, K, X, n){
  ll_diag <- sum(lbeta(diag(alpha_mat) + diag(X), diag(beta_mat) + diag(n) - diag(X)) - lbeta(diag(alpha_mat), diag(beta_mat)))
  ll <- sum(lbeta(alpha_mat + X, beta_mat + n - X)) - sum(lbeta(alpha_mat, beta_mat))
  ll <- (ll + ll_diag)/2
  return(ll)
}

# likelihood of the diagonal part
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

# likelihood of the off-diagonal part
llfunc_offdiag <- function(par, K, X, n){
  alpha <- par[1]
  beta <- par[2]
  if(alpha < 0 | beta < 0) 
    return(Inf)
  else{
    ll <- -K*(K-1)*lbeta(alpha, beta)
    for(i in 1:K)
      for(j in 1:K)
        if(j != i)
          ll <- ll + lbeta(alpha + X[i,j], beta + n[i,j] - X[i,j])
        return(-ll/2)
  }
}

llfunc_sparse <- function(par, K, X, n){
  alpha <- par[1]
  beta <- par[2]
  if(alpha < 0 | beta < 0) 
    return(Inf)
  else{
    ll <- -(K+1)*lbeta(alpha, beta)
    for(i in 1:K)
      ll <- ll + lbeta(alpha + X[i, K+1], beta + n[i, K+1] - X[i, K+1])
    ll <- ll + lbeta(alpha + X[K+1,K+1], beta + n[K+1,K+1] - X[K+1,K+1])
    return(-ll)
  }
}

binormal <- function(x,y,m1,m2,s1,s2,rho){
  exp(-1/(2*(1-rho^2))*((x-m1)^2/s1^2+(y-m2)^2/s2^2-2*rho*(x-m1)*(y-m2)/(s1*s2)))
}

alpha_from_posterior <- function(Y, posterior_matrix, K){
  alpha_mat <- matrix(0, K, K)
  for(i in 1:K){
    for(j in 1:i){
      if(i == j){
        alpha_mat[i, j] = 1 + posterior_matrix[i, ,drop = FALSE] %*% Y %*% as.matrix(posterior_matrix[j, ])/2
      }
      else{
        alpha_mat[i, j] = 1 + posterior_matrix[i, ] %*% Y %*% as.matrix(posterior_matrix[j, ])
        alpha_mat[j, i] = alpha_mat[i, j]
      }
    }
  }
  return(alpha_mat)
}


beta_from_posterior <- function(Y, posterior_matrix, K){
  beta_mat <- matrix(0, K, K)
  Y <- 1 - Y
  diag(Y) <- 0
  for(i in 1:K){
    for(j in 1:i){
      if(i == j){
        beta_mat[i, j] = 1 + posterior_matrix[i, ] %*% Y %*% as.matrix(posterior_matrix[j, ])/2
      }
      else{
        beta_mat[i, j] = 1 + posterior_matrix[i, ] %*% Y %*% as.matrix(posterior_matrix[j, ])
        beta_mat[j, i] = beta_mat[i, j]
      }
    }
  }
  return(beta_mat)
}
