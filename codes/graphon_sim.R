rm(list = ls())
source("lliks.R")

# args = commandArgs(trailingOnly = TRUE)
# N = as.numeric(args[1])
# logrho = as.numeric(args[2])
# lambda = as.numeric(args[3])

N = 100
logrho = -1
lambda = 2  

rho = 10^logrho

W <- function(x, y) {
  rho*lambda^2*(x*y)^(lambda-1)
}

SEED.seq <- 1:100
K.seq <- 1:10
graph_density <- c()

for(SEED in SEED.seq){
  set.seed(SEED)
  
  Theta <-  matrix(0, N, N)
  Y <- matrix(0, N, N)
  U <- sort(runif(N, 0, 1))
  for(i in 1:N){
    for(j in 1:i){
      Theta[j, i] <- min(W(U[i], U[j]), 1)
      Theta[i, j] <- Theta[j, i]
      Y[i, j] <- rbinom(1, 1, Theta[i, j])
      Y[j, i] <- Y[i, j]
    }
  }
  
  diag(Theta) <- 0
  diag(Y) <- 0
  
  graph_density <- c(graph_density, 2*sum(Y)/(N*(N-1)))
  
  res_mixer <- mixer(Y, qmin = min(K.seq), qmax = max(K.seq), method = 'bayesian', directed = FALSE, nbiter = 50, fpnbiter = 5, improve = TRUE)
  par_upper <- c(Inf, Inf)
  
  for(K in K.seq){
    pi_mixer <- res_mixer$output[[K]]$a
    alpha_mixer <- res_mixer$output[[K]]$eta
    beta_mixer <- res_mixer$output[[K]]$zeta
    theta_mixer <- res_mixer$output[[K]]$Pis
    LAB_ALL <- apply(res_mixer$output[[K]]$Taus, 2, which.max)
    LAB_ALL <- as.numeric(factor(LAB_ALL))
    grpsize <- rep(0, K)
    grpsize[1:length(table(LAB_ALL))] <- table(LAB_ALL)
    LLvb <- res_mixer$output[[K]]$criterion
    
    theta_mixer_sort <- as.matrix(theta_mixer[order(apply(theta_mixer, 1, sum)),order(apply(theta_mixer, 1, sum))])
    pi_mixer_sort <- pi_mixer[order(apply(theta_mixer, 1, sum))]
    intv <- c(0, cumsum(pi_mixer_sort)/sum(pi_mixer_sort))
    # MSE of VBEM
    sum_MSE_mixer <- 0
    for(i in 1:K){
      for(j in 1:K){
        MSE_cal <- function(u, v){
          return((W(u, v) - theta_mixer_sort[i, j])^2)
        }
        if(intv[i] < intv[i+1] & intv[j] < intv[j+1]){
          sum_MSE_mixer <- sum_MSE_mixer + int2(MSE_cal, a = c(intv[i], intv[j]), b = c(intv[i+1], intv[j+1]))
        }
      }
    }

    # MLE
    X <- matrix(NA, K, K)
    n <- matrix(NA, K, K)
    
    for(i in 1:K){
      for(j in 1:i){
        # nodes in cluster i and j
        ni <- which(LAB_ALL == i)
        nj <- which(LAB_ALL == j)
        if(length(ni) > 1 && length(nj) > 1){
          if(i == j){
            X[i,j] <- sum(Y[ni,nj])/2
            n[i,j] <- length(ni)*(length(nj) - 1)/2}
          else{
            X[i,j] <- sum(Y[ni,nj])
            n[i,j] <- length(ni)*length(nj)}
          n[j,i] <- n[i,j]
          X[j,i] <- X[i,j]
        }
      }
    }
    
    theta_MLE <- X/n
    for(i in 1:K){
      for(j in 1:K){
        if(is.na(n[i,j]))
          theta_MLE[i,j] <- 1/2
        if(is.na(X[i,j]))
          X[i,j] <- 0
        if(is.na(n[i,j]))
          n[i,j] <- 0
      }
    }
    
    pi_MLE <- grpsize/N
    theta_MLE_sort <- as.matrix(theta_MLE[order(apply(theta_MLE, 1, sum)), order(apply(theta_MLE, 1, sum))])
    pi_MLE <- pi_MLE[order(apply(theta_MLE, 1, sum))]
    intv <- c(0, cumsum(pi_MLE))
    
    sum_MSE_MLE <- 0
    for(i in 1:K){
      for(j in 1:K){
        MSE_cal <- function(u, v){
          return((W(u, v) - theta_MLE_sort[i, j])^2)
        }
        if(intv[i] < intv[i+1] & intv[j] < intv[j+1]){
          sum_MSE_MLE <- sum_MSE_MLE + int2(MSE_cal, a = c(intv[i], intv[j]), b = c(intv[i+1], intv[j+1]))
        }
      }
    }

    # EB
    for(kk in 1:100){
      res <- optimx(c(.001, .001), llfunc_diag, K=K, X=X, n=n,
                    method = "L-BFGS-B", lower = c(1e-5, 1e-5),
                    upper = par_upper)
      llik_diag <- -res$value
      alpha_diag <- res$p1
      beta_diag <- res$p2
      states_diag <- res$convcode
      if(states_diag == 0)
        break
    }
  
    for(kk in 1:100){
      res <- optimx(c(.001, .001), llfunc_offdiag, K=K,  X=X, n=n,
                    method = "L-BFGS-B", lower = c(1e-5, 1e-5),
                    upper = par_upper)
      llik_offdiag <- -res$value
      alpha_offdiag <- res$p1
      beta_offdiag <- res$p2
      states_offdiag <- res$convcode
      if(states_offdiag == 0)
        break
    }
    
    LLIK <- llik_diag + llik_offdiag
    
    pi_prior <- rep(1/2, K)
    # pi_prior <- grpsize/sum(grpsize)
    # pi_prior[pi_prior == 0] <- 1e-20
    zsc <- LLIK + lgamma(sum(pi_prior)) - lgamma(N + sum(pi_prior)) + sum(lgamma(sort(grpsize) + sort(pi_prior)) - lgamma(pi_prior))
    zsc_penal <- zsc - 1/2*(K-1)*log(N) - 1/4*K*(K+1)*log(N*(N-1)/2)
    
    h <- 1/K
    J <- 2/(h*(N-1)) - (N+1)/(h*(N-1))*sum((grpsize/N)^2)
    
    eta_diag <- (alpha_diag + beta_diag)/(alpha_diag + beta_diag + n) 
    eta_offdiag <- (alpha_offdiag + beta_offdiag)/(alpha_offdiag + beta_offdiag + n) 
    eta <- eta_offdiag
    diag(eta) <- diag(eta_diag)
    
    theta_est_diag <- eta_diag*alpha_diag/(alpha_diag+beta_diag) + (1-eta_diag)*theta_MLE
    theta_est_offdiag <- eta_offdiag*alpha_offdiag/(alpha_offdiag+beta_offdiag) + (1-eta_offdiag)*theta_MLE
    diag(theta_est_offdiag) <- diag(theta_est_diag)
    theta_est <- theta_est_offdiag
    
    pi_EB <- grpsize/N
    theta_EB_sort <- as.matrix(theta_est[order(apply(theta_est, 1, sum)), order(apply(theta_est, 1, sum))])
    pi_EB <- pi_EB[order(apply(theta_est, 1, sum))]
    intv <- c(0, cumsum(pi_EB))
    
    sum_MSE_EB <- 0
    for(i in 1:K){
      for(j in 1:K){
        MSE_cal <- function(u, v){
          return((W(u, v) - theta_EB_sort[i, j])^2)
        }
        if(intv[i] < intv[i+1] & intv[j] < intv[j+1]){
          sum_MSE_EB <- sum_MSE_EB + int2(MSE_cal, a = c(intv[i], intv[j]), b = c(intv[i+1], intv[j+1]))
        }
      }
    }
    
    # credible interval
    theta_est_diag_upper <- as.matrix(qbeta(0.975, alpha_diag + X, beta_diag + n - X))
    theta_est_diag_lower <- as.matrix(qbeta(0.025, alpha_diag + X, beta_diag + n - X))
    theta_est_upper <- as.matrix(qbeta(0.975, alpha_offdiag + X, beta_offdiag + n - X))
    theta_est_lower <- as.matrix(qbeta(0.025, alpha_offdiag + X, beta_offdiag + n - X))
    diag(theta_est_upper) <- diag(theta_est_diag_upper)
    diag(theta_est_lower) <- diag(theta_est_diag_lower)

    theta_mixer_upper <- as.matrix(qbeta(0.975, alpha_mixer, beta_mixer))
    theta_mixer_lower <- as.matrix(qbeta(0.025, alpha_mixer, beta_mixer))
    
    Theta_BAYES_UPPER <- matrix(0, N, N)
    Theta_BAYES_LOWER <- matrix(0, N, N)
    Theta_mixer_UPPER <- matrix(0, N, N)
    Theta_mixer_LOWER <- matrix(0, N, N)
    for(i in 1:N){
      for(j in 1:N){
        if(j != i){
          Theta_BAYES_UPPER[i,j] <- theta_est_upper[LAB_ALL[i], LAB_ALL[j]]
          Theta_BAYES_LOWER[i,j] <- theta_est_lower[LAB_ALL[i], LAB_ALL[j]]
          Theta_mixer_UPPER[i,j] <- theta_mixer_upper[LAB_ALL[i], LAB_ALL[j]]
          Theta_mixer_LOWER[i,j] <- theta_mixer_lower[LAB_ALL[i], LAB_ALL[j]]
        }
      }
    }

    cor_EB <- mean((Theta <= Theta_BAYES_UPPER) & (Theta >= Theta_BAYES_LOWER))
    cor_mixer <- mean((Theta <= Theta_mixer_UPPER) & (Theta >= Theta_mixer_LOWER))

    res <- list(grpsize = grpsize, 
                err = c(sum_MSE_MLE, sum_MSE_mixer, sum_MSE_EB),
                llik = LLIK,
                zsc = c(zsc, zsc_penal),
                J = J, 
                ICL = LLvb,
                labels = LAB_ALL,
                eta = eta,
                coverage = c(cor_EB, cor_mixer),
                states = c(states_diag, states_offdiag))
    
    assign(paste0("graphon_", N, "_", -logrho, "_", lambda, "_", K, "_", SEED), res)
    
  }
}

save.image(paste0("./Graphon_N", N, "_logrho", -logrho, "_lambda", lambda, ".RData"))
