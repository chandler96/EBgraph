source("lliks.R")
eucore_net <- read_graph("../datas/EuCore/email-Eu-core.txt",
                         format = "ncol",
                         directed = TRUE)
eucore_net <- as.undirected(eucore_net, "collapse")
# ground truth: 42 communities
LAB_all <- read.table("../datas/EuCore/email-Eu-core-department-labels.txt", row.names = 1)
LAB_all <- LAB_all[,1] + 1
set_vertex_attr(eucore_net, "label", index = V(eucore_net), LAB_all)
Y_all <- as_adjacency_matrix(eucore_net, type = "both",
                             sparse = FALSE)

diag(Y_all) <- 0
N_all <- dim(Y_all)[1]

SEED.seq <- 1:100
for(SEED in SEED.seq){
  set.seed(SEED)
  index <- sample(1:1005, N_all*0.7, replace = FALSE)
  test_index <- (1:N_all)[!(1: N_all) %in% index]
  N <- length(index)
  Y_train <- Y_all[index, index]
  LAB_train <- LAB_all[index]
  tK <- length(table(LAB_train))
  
  LAB_test <- LAB_all[test_index]
  tK_test <- length(table(LAB_test))
  Y_test <- Y_all[test_index, test_index]
  
  K.seq <- 31:50
  res_mixer <- mixer(Y_all, qmin = min(K.seq), qmax = max(K.seq), method = 'bayesian', directed = FALSE, nbiter = 100, fpnbiter = 10, improve = TRUE)
  
  tX <- matrix(NA, 42, 42)
  tn <- matrix(NA, 42, 42)
  
  for(i in 1:42){
    for(j in 1:i){
      # nodes in cluster i and j
      ni <- which(LAB_train == i)
      nj <- which(LAB_train == j)
      if(length(ni) > 1 && length(nj) > 1){
        if(i == j){
          tX[i,j] <- sum(Y[ni,nj])/2
          tn[i,j] <- length(ni)*(length(nj) - 1)/2
        }
        else{
          tX[i,j] <- sum(Y[ni,nj])
          tn[i,j] <- length(ni)*length(nj)
        }
        tX[j,i] <- tX[i,j]
        tn[j,i] <- tn[i,j]
      }
    }
  }
  
  theta_true <- tX/tn
  for(i in 1:K){
    for(j in 1:K){
      if(is.na(n[i,j]))
        theta_true[i,j] <- 1/2
      if(is.na(X[i,j]))
        tX[i,j] <- 0
      if(is.na(n[i,j]))
        tn[i,j] <- 0
    }
  }
  
  PROB_DENSE <- matrix(0, N, N)
  for(i in 1:N){
    for(j in 1:N){
      if(j != i)
        PROB_DENSE[i,j] <- theta_true[LAB_train[i], LAB_train[j]]
    }
  }
  
  for(K in K.seq){
    pi_mixer <- res_mixer$output[[K-min(K.seq)+1]]$a
    alpha_mixer <- res_mixer$output[[K-min(K.seq)+1]]$eta
    beta_mixer <- res_mixer$output[[K-min(K.seq)+1]]$zeta
    theta_mixer <- res_mixer$output[[K-min(K.seq)+1]]$Pis
    posterior_mixer <- res_mixer$output[[K-min(K.seq)+1]]$Taus
    LAB_ALL <- apply(res_mixer$output[[K-min(K.seq)+1]]$Taus, 2, which.max)
    LAB_training <- LAB_ALL[index]
    LAB_testing <- LAB_ALL[test_index]
    grpsize <- rep(0, K)
    grpsize[1:length(table(LAB_ALL))] <- table(LAB_ALL)
    LLvb <- res_mixer$output[[K-min(K.seq)+1]]$criterion
    
    theta_mixer_sort <- as.matrix(theta_mixer[order(apply(theta_mixer, 1, sum)),order(apply(theta_mixer, 1, sum))])
    pi_mixer_sort <- pi_mixer[order(apply(theta_mixer, 1, sum))]
    
    ## MLE
    X <- matrix(NA, K, K)
    n <- matrix(NA, K, K)
    
    for(i in 1:K){
      for(j in 1:i){
        ni <- which(LAB_training == i)
        nj <- which(LAB_training == j)
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

    ## EB
    par_upper <- c(Inf, Inf)
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
    
    # MLE likelihood on testing data
    X_test <- matrix(0, K, K)
    n_test <- matrix(0, K, K)
    
    for(i in 1:K){
      for(j in 1:i){
        # nodes in cluster i and j
        ni <- which(LAB_testing == i)
        nj <- which(LAB_testing == j)
        if(length(ni) > 1 && length(nj) > 1){
          if(i == j){
            X_test[i,j] <- sum(Y[ni,nj])/2
            n_test[i,j] <- length(ni)*(length(nj) - 1)/2}
          else{
            X_test[i,j] <- sum(Y[ni,nj])
            n_test[i,j] <- length(ni)*length(nj)}
          n_test[j,i] <- n_test[i,j]
          X_test[j,i] <- X_test[i,j]
        }
      }
    }
    
    # testing likelihood of MLE
    LLIK_mle <- 0
    for(i in 1:K){
      for(j in 1:i){
        LLIK_mle <- LLIK_mle + X_test[i,j]*log(theta_MLE[i, j]+1e-20) + (n_test[i,j] - X_test[i,j])*(log(1-theta_MLE[i, j]+1e-20))
      }
    }
    
    # testing likelihood of EB
    llik_diag <- -llfunc_diag(c(alpha_diag, beta_diag), K, X_test, n_test)
    llik_offdiag <- -llfunc_diag(c(alpha_offdiag, beta_offdiag), K, X_test, n_test)
    
    LLIK <- llik_diag + llik_offdiag
    
    # testing likelihood of VBEM
    posterior_mixer_train <- posterior_mixer[, index, drop = FALSE]
    alpha_mixer_training <- alpha_from_posterior(as.matrix(Y_train), posterior_mixer_train, K)
    beta_mixer_training <- beta_from_posterior(as.matrix(Y_train), posterior_mixer_train, K)
    
    LLIK_vbem <- llfunc_mixer(alpha_mixer_training, beta_mixer_training, K, X_test, n_test)
    
    pi_prior <- rep(1/2, K)
    zsc <- LLIK + lgamma(sum(pi_prior)) - lgamma(N + sum(pi_prior)) + sum(lgamma(sort(grpsize) + sort(pi_prior)) - lgamma(pi_prior))
    zsc_penal2 <- zsc - 1/2*(K-1)*log(N) - 1/4*K*(K+1)*log(N*(N-1)/2)
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
    
    # theta_est_diag_upper <- as.matrix(qbeta(0.975, alpha_diag + X, beta_diag + n - X))
    # theta_est_diag_lower <- as.matrix(qbeta(0.025, alpha_diag + X, beta_diag + n - X))
    # theta_est_upper <- as.matrix(qbeta(0.975, alpha_offdiag + X, beta_offdiag + n - X))
    # theta_est_lower <- as.matrix(qbeta(0.025, alpha_offdiag + X, beta_offdiag + n - X))
    # diag(theta_est_upper) <- diag(theta_est_diag_upper)
    # diag(theta_est_lower) <- diag(theta_est_diag_lower)
    
    # theta_mixer_upper <- as.matrix(qbeta(0.975, alpha_mixer, beta_mixer))
    # theta_mixer_lower <- as.matrix(qbeta(0.025, alpha_mixer, beta_mixer))
    
    # pi_EB <- grpsize/N
    # theta_EB_sort_upper <- as.matrix(theta_est_upper[order(apply(theta_est, 1, sum)), order(apply(theta_est, 1, sum))])
    # theta_EB_sort_lower <- as.matrix(theta_est_lower[order(apply(theta_est, 1, sum)), order(apply(theta_est, 1, sum))])
    # pi_EB <- pi_EB[order(apply(theta_est, 1, sum))]
    
    # theta_mixer_sort_upper <- as.matrix(theta_mixer_upper[order(apply(theta_mixer, 1, sum)), order(apply(theta_mixer, 1, sum))])
    # theta_mixer_sort_lower <- as.matrix(theta_mixer_lower[order(apply(theta_mixer, 1, sum)), order(apply(theta_mixer, 1, sum))])
    
    Theta_MLE <- matrix(0, N, N)
    Theta_mixer <- matrix(0, N, N)
    Theta_BAYES <- matrix(0, N, N)
    for(i in 1:N){
      for(j in 1:N){
        if(j != i){
          Theta_MLE[i,j] <- theta_MLE[LAB_ALL[i], LAB_ALL[j]]
          Theta_mixer[i,j] <- theta_mixer[LAB_ALL[i], LAB_ALL[j]]
          Theta_BAYES[i,j] <- theta_est[LAB_ALL[i], LAB_ALL[j]]
        }
      }
    }

    err.MLE <- mean((Theta_MLE - PROB_DENSE)^2)
    err.BAYES <- mean((Theta_BAYES - PROB_DENSE)^2)
    err.mixer <- mean((Theta_mixer - PROB_DENSE)^2)

    
    res <- list(grpsize = grpsize, 
                err = c(err.MLE, err.BAYES, err.mixer),
                llik = c(LLIK, LLIK_mle, LLIK_vbem),
                zsc = c(zsc, zsc_penal),
                J = J,
                ICL = LLvb,
                labels = LAB_ALL,
                eta = eta,
                states = c(states_diag, states_offdiag))
    
    assign(paste0("eucore_", K, "_", SEED), res)
  }
}

save.image(paste0("./SBMres/real/eucore.RData"))



