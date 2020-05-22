rm(list = ls())
source("lliks.R")

# args = commandArgs(trailingOnly = TRUE)
# N = as.numeric(args[1])
# tK = as.numeric(args[2])
# rho = as.numeric(args[3])
N = 100
tK = 10
rho = 1
B = matrix(0.1, tK, tK)
diag(B) <- 0.9
B = B * rho
  
SEED.seq <- 1:100
K.seq <- 1:20
  
set.seed(0)
LAB_true <- sample(x = 1:tK, size = N , replace = TRUE, prob = rep(1/tK, tK))
  
graph_density <- c()
  
for(SEED in SEED.seq){
  set.seed(SEED)
  Theta <- matrix(0, N, N)
  Y <- matrix(0, N, N)
  for(i in 1:N){
    for(j in 1:i){
      Theta[j, i] <- B[LAB_true[i], LAB_true[j]]
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
    alpha_mixer <- res_mixer$output[[K]]$eta
    beta_mixer <- res_mixer$output[[K]]$zeta
    theta_mixer <- res_mixer$output[[K]]$Pis
    LAB_ALL <- apply(res_mixer$output[[K]]$Taus, 2, which.max)
    grpsize <- rep(0, K)
    grpsize[1:length(table(LAB_ALL))] <- table(LAB_ALL)
    LLvb <- res_mixer$output[[K]]$criterion
    rand_index <- rand.index(LAB_ALL, LAB_true)
      
    X <- matrix(NA, K, K)
    n <- matrix(NA, K, K)
    theta <- matrix(NA, K, K)
      
    for(i in 1:K){
      for(j in 1:i){
        # nodes in cluster i and j
        ni <- which(LAB_ALL == i)
        nj <- which(LAB_ALL == j)
        if(length(ni) > 1 && length(nj) > 1){
          if(i == j){
            X[i,j] <- sum(Y[ni,nj])/2
            n[i,j] <- length(ni)*(length(nj) - 1)/2
            }
          else{
            X[i,j] <- sum(Y[ni,nj])
            n[i,j] <- length(ni)*length(nj)
            }
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
    
    # for(kk in 1:10){
    #   res <- optimx(c(.001, .001), llfunc, K=K, X=X, n=n,
    #                 method = "L-BFGS-B", lower = c(1e-5, 1e-5),
    #                 upper = par_upper)
    #   llik_WHOLE <- -res$value
    #   alpha_WHOLE <- res$p1
    #   beta_WHOLE <- res$p2
    #   states_WHOLE <- res$convcode
    #   if(states_WHOLE == 0)
    #     break
    # }
      
    # diagonal part  
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
      
    # off-diagonal part
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
      
    LLIK_mle <- 0
    for(i in 1:K){
      for(j in 1:i){
        LLIK_mle <- LLIK_mle + X[i,j]*log(theta_MLE[i, j]+1e-20) + (n[i,j] - X[i,j])*(log(1-theta_MLE[i, j]+1e-20))
      }
    }
      
    LLIK <- llik_diag + llik_offdiag
      
    pi_prior <- rep(1/2, K)
    # pi_prior <- grpsize/sum(grpsize)
    # pi_prior[pi_prior == 0] <- 1e-20
    zsc <- LLIK + lgamma(sum(pi_prior)) - lgamma(N + sum(pi_prior)) + sum(lgamma(sort(grpsize) + sort(pi_prior)) - lgamma(pi_prior))
    zsc_penal <- zsc - 1/2*(K-1)*log(N) - 1/4*K*(K+1)*log(N*(N-1)/2)
      
    h <- 1/K
    J <- 2/(h*(N-1)) - (N+1)/(h*(N-1))*sum((grpsize/N)^2)
      
    eta_WHOLE <- (alpha_WHOLE + beta_WHOLE)/(alpha_WHOLE + beta_WHOLE + n)
    theta_est_WHOLE <- (eta_WHOLE*alpha_WHOLE)/(alpha_WHOLE + beta_WHOLE) + (1- eta_WHOLE)*theta_MLE
    eta_diag <- (alpha_diag + beta_diag)/(alpha_diag + beta_diag + n) 
    eta_offdiag <- (alpha_offdiag + beta_offdiag)/(alpha_offdiag + beta_offdiag + n) 
    eta <- eta_offdiag
    diag(eta) <- diag(eta_diag)
      
    theta_est_diag <- eta_diag*alpha_diag/(alpha_diag+beta_diag) + (1-eta_diag)*theta_MLE
    theta_est_offdiag <- eta_offdiag*alpha_offdiag/(alpha_offdiag+beta_offdiag) + (1-eta_offdiag)*theta_MLE
    diag(theta_est_offdiag) <- diag(theta_est_diag)
    theta_est <- theta_est_offdiag
      
    # the upper and lower bond
    theta_est_diag_upper <- as.matrix(qbeta(0.975, alpha_diag + X, beta_diag + n - X))
    theta_est_diag_lower <- as.matrix(qbeta(0.025, alpha_diag + X, beta_diag + n - X))
    theta_est_upper <- as.matrix(qbeta(0.975, alpha_offdiag + X, beta_offdiag + n - X))
    theta_est_lower <- as.matrix(qbeta(0.025, alpha_offdiag + X, beta_offdiag + n - X))
    diag(theta_est_upper) <- diag(theta_est_diag_upper)
    diag(theta_est_lower) <- diag(theta_est_diag_lower)
      
    # lower and upper bound for mixer
    theta_mixer_upper <- as.matrix(qbeta(0.975, alpha_mixer, beta_mixer))
    theta_mixer_lower <- as.matrix(qbeta(0.025, alpha_mixer, beta_mixer))
      
    # the upper and lower bond of estimates for each pair of the nodes
    PROB_DENSE <- Theta
    Theta_MLE <- matrix(0, N, N)
    Theta_mixer <- matrix(0, N, N)
    Theta_BAYES <- matrix(0, N, N)
    Theta_BAYES_UPPER <- matrix(0, N, N)
    Theta_BAYES_LOWER <- matrix(0, N, N)
    Theta_mixer_UPPER <- matrix(0, N, N)
    Theta_mixer_LOWER <- matrix(0, N, N)
    Theta_WHOLE <- matrix(0, N, N)
    for(i in 1:N){
      for(j in 1:N){
        if(j != i){ 
          Theta_MLE[i,j] <- theta_MLE[LAB_ALL[i], LAB_ALL[j]]
          Theta_mixer[i,j] <- theta_mixer[LAB_ALL[i], LAB_ALL[j]]
          Theta_BAYES[i,j] <- theta_est[LAB_ALL[i], LAB_ALL[j]]
          Theta_BAYES_UPPER[i,j] <- theta_est_upper[LAB_ALL[i], LAB_ALL[j]]
          Theta_BAYES_LOWER[i,j] <- theta_est_lower[LAB_ALL[i], LAB_ALL[j]]
          Theta_mixer_UPPER[i,j] <- theta_mixer_upper[LAB_ALL[i], LAB_ALL[j]]
          Theta_mixer_LOWER[i,j] <- theta_mixer_lower[LAB_ALL[i], LAB_ALL[j]]
          Theta_WHOLE[i,j] <- theta_est_WHOLE[LAB_ALL[i], LAB_ALL[j]]
        }
      }
    }
      
    err.MLE <- mean((Theta_MLE - PROB_DENSE)^2)
    err.BAYES <- mean((Theta_BAYES - PROB_DENSE)^2)
    err.WHOLE <- mean((Theta_WHOLE - PROB_DENSE)^2)
    err.mixer <- mean((Theta_mixer - PROB_DENSE)^2)
      
    coverage <- mean((PROB_DENSE <= Theta_BAYES_UPPER) & (PROB_DENSE >= Theta_BAYES_LOWER))
    coverage.mixer <- mean((PROB_DENSE <= Theta_mixer_UPPER) & (PROB_DENSE >= Theta_mixer_LOWER))
    
    res <- list(grpsize = grpsize, 
                err = c(err.MLE, err.BAYES, err.WHOLE, err.mixer),
                llik = c(LLIK, LLIK_mle),
                zsc = c(zsc, zsc_penal),
                J = J, 
                ICL = LLvb,
                labels = LAB_ALL,
                recovery = rand_index,
                eta = eta,
                coverage = c(coverage, coverage.mixer),
                states = c(states_diag, states_offdiag))
      
      assign(paste0("SBM_", N, "_", SEED, "_", K), res)
    }
  }
  
save.image(paste0("./SBM_N", N, "_tK", tK, "_rho", rho, ".RData"))
