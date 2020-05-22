rm(list = ls())
source("common_functions.R")
try(data(package = "mixer"))
data(blog)
links <- blog$links
politicalParty <- blog$politicalParty
politicalParty[165] <- politicalParty[151]
politicalParty <- factor(politicalParty)
LAB_all<- as.numeric(politicalParty)
N_all <- dim(links)[1]
Y_all <- links
K.seq <- 1:20
SEED.seq <- 1:100
lliks <- c()
ps <- c()
for(SEED in SEED.seq){
  set.seed(SEED)
  index <- sort(sample(1:N_all, 138, replace = FALSE))
  test_index <- (1:N_all)[!(1: N_all) %in% index]
  N <- length(index)
  Y_train <- Y_all[index, index]
  LAB_train <- LAB_all[index]
  tK <- length(table(LAB_train))
  
  LAB_test <- LAB_all[test_index]
  tK_test <- length(table(LAB_test))
  Y_test <- Y_all[test_index, test_index]

  res_mixer <- mixer(Y_all, qmin = min(K.seq), qmax = max(K.seq), method = 'bayesian', directed = FALSE, nbiter = 100, fpnbiter = 10, improve = TRUE)

  # find the ground truth for training only
  tX <- matrix(0, 10, 10)
  tn <- matrix(0, 10, 10)
  
  for(i in 1:10){
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
      }
    }
  }
  
  theta_true <- tX/tn
  theta_true[is.na(theta_true)] <- 0
  PROB_DENSE <- matrix(0, N, N)
  
  for(i in 1:N){
    for(j in 1:N){
      if(j != i)
        PROB_DENSE[i,j] <- theta_true[LAB_train[i], LAB_train[j]]
    }
  }
  
  for(K in K.seq){
    ## results of mixer
    pi_mixer <- res_mixer$output[[K]]$a
    alpha_mixer <- res_mixer$output[[K]]$eta
    beta_mixer <- res_mixer$output[[K]]$zeta
    theta_mixer <- res_mixer$output[[K]]$Pis
    posterior_mixer <- res_mixer$output[[K]]$Taus
    LAB_ALL <- apply(res_mixer$output[[K]]$Taus, 2, which.max)
    LAB_training <- LAB_ALL[index]
    LAB_testing <- LAB_ALL[test_index]
    grpsize <- rep(0, K)
    grpsize[1:length(table(LAB_ALL))] <- table(LAB_ALL)
    LLvb <- res_mixer$output[[K]]$criterion
    rand_index <- rand.index(LAB_ALL, LAB_all)
    
    prop <- mean(!LAB_testing %in% LAB_training)
    print(prop)
    ps <- c(prop, ps)
    
    ## results of MLE
    ## estimate on the dense part
    X <- matrix(NA, K, K)
    n <- matrix(NA, K, K)
    
    for(i in 1:K){
      for(j in 1:i){
        # nodes in cluster i and j
        ni <- which(LAB_training == i)
        nj <- which(LAB_training == j)
        if(length(ni) > 1 && length(nj) > 1){
          if(i == j){
            X[i,j] <- sum(Y_train[ni,nj])/2
            n[i,j] <- length(ni)*(length(nj) - 1)/2}
          else{
            X[i,j] <- sum(Y_train[ni,nj])
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
    
    ## Results of EB
    par_upper <- c(Inf, Inf)
    for(kk in 1:10){
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
    
    # for the off-diagonal part
    for(kk in 1:10){
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
            X_test[i,j] <- sum(Y_test[ni,nj])/2
            n_test[i,j] <- max(length(ni)*(length(nj) - 1)/2, 1)}
          else{
            X_test[i,j] <- sum(Y_test[ni,nj])
            n_test[i,j] <- max(length(ni)*length(nj) , 1)}
          n_test[j,i] <- n_test[i,j]
          X_test[j,i] <- X_test[i,j]
        }
      }
    }
    
    # training likelihood of MLE
    eta_diag <- (alpha_diag + beta_diag)/(alpha_diag + beta_diag + n)
    eta_offdiag <- (alpha_offdiag + beta_offdiag)/(alpha_offdiag + beta_offdiag + n)
    eta <- eta_offdiag
    diag(eta) <- diag(eta_diag)
    
    theta_est_diag <- eta_diag*alpha_diag/(alpha_diag+beta_diag) + (1-eta_diag)*theta_MLE
    theta_est_offdiag <- eta_offdiag*alpha_offdiag/(alpha_offdiag+beta_offdiag) + (1-eta_offdiag)*theta_MLE
    diag(theta_est_offdiag) <- diag(theta_est_diag)
    theta_est <- theta_est_offdiag
    
    # testing likelihood of MLE
    LLIK_mle <- 0 # -10652.76
    LLIK2_EB <- 0 # -4102.95
    LLIK2_VBEM <- 0 # -3705.615
    
    N_test = N_all - N
    for(i in test_index){
      # connection to training data
      for(j in index){
        c1 <- LAB_ALL[i] # community 1
        c2 <- LAB_ALL[j] # community 2
        LLIK_mle <- LLIK_mle + Y_all[i, j]*log(theta_MLE[c1, c2] + 1e-20) + (1-Y_all[i, j])*log(1-theta_MLE[c1, c2] + 1e-20)
        LLIK2_VBEM <- LLIK2_VBEM + Y_all[i, j]*log(theta_mixer[c1, c2] + 1e-20) + (1-Y_all[i, j])*log(1-theta_mixer[c1, c2] + 1e-20)
        LLIK2_EB <- LLIK2_EB + Y_all[i, j]*log(theta_est[c1, c2] + 1e-20) + (1-Y_all[i, j])*log(1-theta_est[c1, c2] + 1e-20)
      }
      for(j in test_index){
        if(i != j){
          c1 <- LAB_ALL[i] # community 1
          c2 <- LAB_ALL[j] # community 2
          LLIK_mle <- LLIK_mle + (Y_all[i, j]*log(theta_MLE[c1, c2] + 1e-20) + (1-Y_all[i, j])*log(1-theta_MLE[c1, c2] + 1e-20))/2
          LLIK2_VBEM <- LLIK2_VBEM + (Y_all[i, j]*log(theta_mixer[c1, c2] + 1e-20) + (1-Y_all[i, j])*log(1-theta_mixer[c1, c2] + 1e-20))/2
          LLIK2_EB <- LLIK2_EB + (Y_all[i, j]*log(theta_est[c1, c2] + 1e-20) + (1-Y_all[i, j])*log(1-theta_est[c1, c2] + 1e-20))/2
        }
      }
    }
    
    #testing likelihood of EB
    llik_diag <- -llfunc_diag(c(alpha_diag, beta_diag), K, X_test, n_test)
    llik_offdiag <- -llfunc_offdiag(c(alpha_offdiag, beta_offdiag), K, X_test, n_test)
    #
    LLIK <- llik_diag + llik_offdiag
    #
    # testing likelihood of VBEM
    posterior_mixer_train <- posterior_mixer[, index, drop = FALSE]
    alpha_mixer_training <- alpha_from_posterior(as.matrix(Y_train), posterior_mixer_train, K)
    beta_mixer_training <- beta_from_posterior(as.matrix(Y_train), posterior_mixer_train, K)

    LLIK_vbem <- llfunc_mixer(alpha_mixer_training, beta_mixer_training, K, X_test, n_test)
    LLvb <- LLvb - LLIK_vbem
    LLIK_train <- llfunc_mixer(alpha_mixer_training, beta_mixer_training, K, X, n)

    #
    pi_prior <- rep(1/2, K)
    #pi_prior <- grpsize/N
    #pi_prior[pi_prior] <- 1e-20
    #zsc_WHOLE <- LLIK_WHOLE + sum(lgamma(grpsize+1/K)) - K*lgamma(1/K) - lgamma(N+1)
    zsc <- LLIK + lgamma(sum(pi_prior)) - lgamma(N + sum(pi_prior)) + sum(lgamma(sort(grpsize) + sort(pi_prior)) - lgamma(pi_prior))
    #zsc <- LLIK + + sum(lgamma(grpsize+1/K)) - K*lgamma(1/K) - lgamma(N+1)
    zsc_penal <- zsc - 1/2*(K-1)*log(N)
    zsc_penal2 <- zsc - 1/2*(K-1)*log(N) - 1/4*K*(K+1)*log(N*(N-1)/2)

    # # dense part estimator
    h <- 1/K
    J <- 2/(h*(N-1)) - (N+1)/(h*(N-1))*sum((grpsize/N)^2)

    # eta_diag <- (alpha_diag + beta_diag)/(alpha_diag + beta_diag + n)
    # eta_offdiag <- (alpha_offdiag + beta_offdiag)/(alpha_offdiag + beta_offdiag + n)
    # eta <- eta_offdiag
    # diag(eta) <- diag(eta_diag)

    # theta_est_diag <- eta_diag*alpha_diag/(alpha_diag+beta_diag) + (1-eta_diag)*theta_MLE
    # theta_est_offdiag <- eta_offdiag*alpha_offdiag/(alpha_offdiag+beta_offdiag) + (1-eta_offdiag)*theta_MLE
    # diag(theta_est_offdiag) <- diag(theta_est_diag)
    # theta_est <- theta_est_offdiag
    #
    # ## RMSE and MAE for EB
    # 
    # pi_EB <- grpsize/N
    # theta_EB_sort <- as.matrix(theta_est[order(apply(theta_est, 1, sum)), order(apply(theta_est, 1, sum))])
    # pi_EB <- pi_EB[order(apply(theta_est, 1, sum))]
    # intv <- c(0, cumsum(pi_EB))
    # 
    # # credible interval
    # theta_est_diag_upper <- as.matrix(qbeta(0.975, alpha_diag + X, beta_diag + n - X))
    # theta_est_diag_lower <- as.matrix(qbeta(0.025, alpha_diag + X, beta_diag + n - X))
    # theta_est_upper <- as.matrix(qbeta(0.975, alpha_offdiag + X, beta_offdiag + n - X))
    # theta_est_lower <- as.matrix(qbeta(0.025, alpha_offdiag + X, beta_offdiag + n - X))
    # diag(theta_est_upper) <- diag(theta_est_diag_upper)
    # diag(theta_est_lower) <- diag(theta_est_diag_lower)
    # 
    # ## lower and upper bound for mixer
    # theta_mixer_upper <- as.matrix(qbeta(0.975, alpha_mixer, beta_mixer))
    # theta_mixer_lower <- as.matrix(qbeta(0.025, alpha_mixer, beta_mixer))
    # 
    # ## RMSE and MAE for EB
    # pi_EB <- grpsize/N
    # ort_upper <- as.matrix(theta_est_upper[order(apply(theta_est, 1, sum)), order(apply(theta_est, 1, sum))])
    # theta_EB_sort_lower <- as.matrix(theta_est_lower[order(apply(theta_est, 1, sum)), order(apply(theta_est, 1, sum))])
    # pi_EB <- pi_EB[order(apply(theta_est, 1, sum))]
    # 
    # theta_mixer_sort_upper <- as.matrix(theta_mixer_upper[order(apply(theta_mixer, 1, sum)), order(apply(theta_mixer, 1, sum))])
    # theta_mixer_sort_lower <- as.matrix(theta_mixer_lower[order(apply(theta_mixer, 1, sum)), order(apply(theta_mixer, 1, sum))])
    # 
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

    err.MLE <- sqrt(mean((Theta_MLE - PROB_DENSE)^2))
    err.BAYES <- sqrt(mean((Theta_BAYES - PROB_DENSE)^2))
    err.mixer <- sqrt(mean((Theta_mixer - PROB_DENSE)^2))
    # 
    res <- list(grpsize = grpsize, err = c(err.MLE, err.BAYES, err.mixer), llik_marg = c(LLIK, LLIK_mle, LLIK_vbem),
                llik = c(LLIK_mle, LLIK2_EB, LLIK2_VBEM) , zsc = c(zsc, zsc_penal, zsc_penal2),
                J = J,
                ICL = LLvb,
                #labels = LAB_ALL,
                #theta = theta_EB_sort,
                #Theta = Theta_BAYES,
                #rand_index = tK,
                #eta = eta,
                states = c(states_diag, states_offdiag))
    #posterior_mixer_train <- posterior_mixer[, index, drop = FALSE]
    #alpha_mixer_training <- alpha_from_posterior(as.matrix(Y_train), posterior_mixer_train, K)
    #beta_mixer_training <- beta_from_posterior(as.matrix(Y_train), posterior_mixer_train, K)
    
    #LLIK_vbem <- llfunc_mixer(alpha_mixer_training, beta_mixer_training, K, X_test, n_test)
    #res <- c(LLIK_mle, LLIK2_EB, LLIK2_VBEM)
    #lliks <- rbind(lliks, res)

    print(c(K, SEED, res$err))
    # 
    #res <- c(tK, LLIK, LLIK_mle, zsc, zsc_penal, zsc_penal2, LLvb, LLIK_train, J, LLIK_mle)
    #print(c(K, res))
    assign(paste0("blog_", K, "_", SEED), res)
  }
}

getstat <- function(dat){
  return(list(err = dat$err,
              llik = dat$llik,
              zscs = dat$zsc,
              J = dat$J,
              llik = dat$llik,
              icl = dat$ICL,
              label = dat$labels,
              tK = dat$rand_index,
              eta_mat = dat$eta,
              states = dat$states,
              theta = dat$theta,
              Theta = dat$Theta,
              grp = dat$grpsize))
}

SEED.seq <- 1:100
E1.seed <- matrix(0, length(SEED.seq), length(K.seq))
E2.seed <- matrix(0, length(SEED.seq), length(K.seq))
E3.seed <- matrix(0, length(SEED.seq), length(K.seq))
LLIK.seed <- matrix(0, length(SEED.seq), length(K.seq))
LLIK.mle.seed <- matrix(0, length(SEED.seq), length(K.seq))
LLIK.mixer.seed <- matrix(0, length(SEED.seq), length(K.seq))
L1.seed <- matrix(0, length(SEED.seq), length(K.seq)) # EB
L2.seed <- matrix(0, length(SEED.seq), length(K.seq)) # MLE
L3.seed <- matrix(0, length(SEED.seq), length(K.seq)) # VBEM
zsc.seed <- matrix(0, length(SEED.seq), length(K.seq))
zsc.penal.seed <- matrix(0, length(SEED.seq), length(K.seq))
zsc.penal2.seed <- matrix(0, length(SEED.seq), length(K.seq))
nzsc.seed <- matrix(0, length(SEED.seq), length(K.seq))
nzsc.penal.seed <- matrix(0, length(SEED.seq), length(K.seq))
nzsc.penal2.seed <- matrix(0, length(SEED.seq), length(K.seq))
J.seed <- matrix(0, length(SEED.seq), length(K.seq))
tK.seed <- matrix(0, length(SEED.seq), length(K.seq))
ICL.seed <- matrix(0, length(SEED.seq), length(K.seq))
grp.seed <- array(dim = c(length(SEED.seq), 8, length(K.seq)))
theta.seed <- array(dim = c(length(SEED.seq), 8, 8, length(K.seq)))
#Theta.seed <- array(dim = c(length(SEED.seq), N, N, length(K.seq)))

for(SEED in SEED.seq){
  for(K in K.seq){
    res <- getstat(eval(parse(text = paste0("blog_", K, "_", SEED))))
    index <- SEED
    offset <- 0
    E1.seed[index, K-offset] <- res$err[1]^2 # MLE
    E2.seed[index, K-offset] <- res$err[3]^2 # mixer
    E3.seed[index, K-offset] <- res$err[2]^2 # EB
    # LLIK.seed[index, K-offset] <- res$llik[3] # Bayes likelihood
    # LLIK.mle.seed[index, K-offset] <- res$llik[1]
    # LLIK.mixer.seed[index, K-offset] <- res$llik[2]
    #L1.seed[index, K-offset] <- res$llik[1] #MLE
    #L2.seed[index, K-offset] <- res$llik[2] #EB
    #L3.seed[index, K-offset] <- res$llik[3] #VBEM
    #zsc.seed[index, K-offset] <- res$zsc[1] # Bayes joint likelihood
    #zsc.penal.seed[index, K-offset] <- res$zsc[2] # Penalized Bayes likelihood
    #zsc.penal2.seed[index, K-offset] <- res$zsc[3] # icl Penalized Bayes likelihood
    # tK.seed[index, K-offset] <- res$tK
    #J.seed[index, K-offset] <- res$J
    #ICL.seed[index, K-offset] <- res$icl
    # theta.seed[index, , ,K-offset] <- res$theta
    # #Theta.seed[index, , ,K-offset] <- res$Theta
    # grp.seed[index, K-offset] <- res$grp
  }
}


for(SEED in SEED.seq){
  for(K in K.seq){
    res <- eval(parse(text = paste0("blog_", K, "_", SEED)))
    index <- SEED
    offset <- 0
    zsc.seed[index, K-offset] <- res[1] # zsc
    zsc.penal.seed[index, K-offset] <- res[2] # zsc.penal
    zsc.penal2.seed[index, K-offset] <- res[3] # zsc.penal2
    nzsc.seed[index, K-offset] <- res[4] # nzsc
    nzsc.penal.seed[index, K-offset] <- res[5] # nzsc.penal
    nzsc.penal2.seed[index, K-offset] <- res[6] # nzsc.penal2
    J.seed[index, K-offset] <- res[7] # J
    ICL.seed[index, K-offset] <- res[8] # LLvb
  }
}

apply(zsc.penal2.seed, 1, which.max)
apply(ICL.seed, 1, which.max)
which.max(apply(zsc.penal2.seed, 2, mean))
boxplot(zsc.penal2.seed)
apply(zsc.penal2.seed,1,which.max)


theta_EB_est <- apply(theta.seed[,,,8], c(2,3), mean)
# SD1 <- apply(L1.seed, 2, sd)
# SD2 <- apply(L2.seed, 2, sd)
# SD3 <- apply(L3.seed, 2, sd)
library(plotly)
p <- plot_ly(z = ~theta_EB_est) %>% add_surface()
p

setEPS()
accuracy_length = 300
postscript("./figs_real_data/poliblog_llik3.eps")
dat2 <- data.frame(x = K.seq, L1 = apply(L2.seed, 2, mean),
                   L2 = apply(L1.seed, 2, mean),
                   L3 = apply(L3.seed, 2, mean))
cols <- c("EB" = "#619CFF", "MLE" = "#F8766D", "VBEM" = "#00BA38")
p2 <- ggplot(dat2, aes(factor(x), group = 1)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  geom_line(aes(y = L1, col = "EB"), size = 1.8) +
  geom_line(aes(y = L2, col = "MLE"), size = 1.8) +
  geom_line(aes(y = L3, col = "VBEM"), size = 1.8) +
  #geom_line(aes(y = L4, col = "EB2"), size = 1.8, linetype = "dashed") +
  #geom_line(aes(y = L5, col = "VBEM2"), size = 1.8,linetype = "dashed") +
  theme(axis.title = element_text(size = 25)) + 
  theme(axis.text = element_text(size = 18)) + 
  theme(axis.text.y = element_text(face = "italic", color = "black")) +
  theme(axis.text.x = element_text(face = "italic", color = "black")) +
  #theme(axis.text.x = element_text(face = ifelse(dat2$x == 10, "bold", "italic"),
  #geom_vline(xintercept = which.min(dat2$E3), linetype = "dashed") +
  # geom_errorbar(aes(ymin = L1-SD1,
  #                   ymax = L1+SD1, col = "EB"),
  #               width = .2, position = position_dodge(0.05)) +
  # geom_errorbar(aes(ymin = L2-SD2,
  #                   ymax = L2+SD2, col = "MLE"),
  #               width = .2, position = position_dodge(0.05)) +
  # geom_errorbar(aes(ymin = L3-SD3,
  #                   ymax = L3+SD3, col = "VBEM"),
  #               width = .2, position = position_dodge(0.05)) +
  scale_color_manual(values = cols) +
  labs(x = "") +
  labs(y = "") 
  # annotate(geom = 'text', x = 17, y = -150, label = 'EB', color = 'black', size = 5) +
  # annotate(geom = 'text', x = 17, y = -350, label = 'VBEM', color = 'black', size = 5) +
  # annotate(geom = 'text', x = 17, y = -630, label = 'MLE', color = 'black', size = 5)
p2
dev.off()

setEPS()
postscript(paste0("./plots0830/blog_marg_llik.eps"))
p2
dev.off()


R31 <- apply(E3.seed,2,mean)/apply(E1.seed,2,mean)*100
R32 <- apply(E3.seed,2,mean)/apply(E2.seed,2,mean)*100
# sd_cal <- function(S){
#   q1 <- quantile(S, 0.25)
#   q3 <- quantile(S, 0.75)
#   S_sel <- S[S >= q1 & S <= q3]
#   return(sd(S_sel))
# }
# 
# m_cal <- function(S){
#   q1 <- quantile(S, 0.25)
#   q3 <- quantile(S, 0.75)
#   S_sel <- S[S >= q1 & S <= q3]
#   return(mean(S_sel))
# }
# 
# sd31 <- apply(R31, 2, sd_cal)
# sd32 <- apply(R32, 2, sd_cal)

shows <- c(5,10,15,20)
xlabels <- rep("", 20)
xlabels[shows] <- shows
setEPS()
accuracy_length = 400
postscript("./figs_real_data/0509/poliblog_MSE_ratio.eps")
dat1 <- data.frame(x = factor(1:20), 
                   I1 = R31,
                   I2 = R32)
cols <- c("EB/MLE" = "#F8766D", "EB/VBEM" = "#619CFF")
p1 <- ggplot(dat1, aes(x, group = 1)) + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") + 
  geom_line(aes(y = I1, col = "EB/MLE"),size=1.8) +
  geom_line(aes(y = I2, col = "EB/VBEM"),size=1.8) +
  geom_point(aes(y = I1, col = "EB/MLE"),size=1.8) +
  geom_point(aes(y = I2, col = "EB/VBEM"),size=1.8) +
  scale_x_discrete(breaks = 1:20,
                   labels = xlabels) +
  theme(axis.title = element_text(size = 25)) + 
  theme(axis.text = element_text(size = 18)) + 
  theme(axis.text.y = element_text(face = "italic", color = "black",size = 40)) +
  theme(axis.text.x = element_text(face = "italic", color = "black",size = 40)) +
  #                                 color = ifelse(dat1$x == 10, "red", "black"))) + 
  # geom_errorbar(aes(ymin = I1-sd31,
  #                   ymax = I1+sd31, col = "EB/MLE"),
  #               width = .2, position = position_dodge(0.05)) +
  # geom_errorbar(aes(ymin = I2-sd32,
  #                   ymax = I2+sd32, col = "EB/VBEM"),
  #               width = .2, position = position_dodge(0.05)) +
  theme(axis.ticks.length.x.bottom = unit(0.5, "cm")) +
  #theme(axis.ticks.x.bottom = element_line(colour = ifelse(dat1$x ==tK, "red", "black"))) +
  theme(axis.ticks.x.bottom = element_line(size = ifelse(dat1$x %in% c(5,10,15,20), 2, .5))) +
  theme(axis.ticks.length.y.left = unit(0.5, "cm")) +
  scale_color_manual(values = cols) + 
  labs(x = "") + 
  labs(y = "") +
  theme(plot.margin=unit(c(0,0.5,0,0),"cm"))+
  # annotate(geom = 'text', x = 15, y = 87, label = 'EB/VBEM', color = 'black') +
  # annotate(geom = 'text', x = 15, y = 97.5, label = 'EB/MLE', color = 'black') +
  geom_hline(yintercept=100, linetype="dashed")
p1
dev.off()


setEPS()
accuracy_length = 300
postscript("./figs_real_data/poliblog_MSE.eps")
dat2 <- data.frame(x = 1:20, E1 = apply(E1.seed, 2, mean), 
                   E2 = apply(E2.seed, 2, mean), 
                   E3 = apply(E3.seed, 2, mean))
cols <- c("MLE" = "#F8766D", "VBEM" = "#00BA38", "EB" = "#619CFF")
p2 <- ggplot(dat2, aes(factor(x), group = 1)) +
  theme_classic() +
  #theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  geom_line(aes(y = E1, col = "MLE")) +
  geom_line(aes(y = E2, col = "VBEM")) +
  geom_line(aes(y = E3, col = "EB")) +
  scale_color_manual(values = cols) + 
  labs(x = "") + 
  labs(y = "")
p2
dev.off()

setEPS()
postscript("./plots0830/blog_mse.eps")
grid.arrange(p2, p1, ncol = 1)
dev.off()

myscale <- function(x){
  avg <- x
  (avg-min(avg))/(max(avg)-min(avg))
}


dat3 <- data.frame(x = 1:20, E1 = myscale(apply(ICL.seed, 2, mean)), 
                   E2 = myscale(apply(J.seed, 2, mean)), 
                   E3 = myscale(apply(zsc.penal2.seed, 2, mean)))
cols <- c("Ilvb" = "#F8766D", "J" = "#00BA38", "Llik" = "#619CFF")
p3 <- ggplot(dat3, aes(factor(x), group = 1)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  #theme(legend.position = "none")
  geom_line(aes(y = E1, col = "Ilvb")) +
  geom_line(aes(y = E2, col = "J")) +
  geom_line(aes(y = E3, col = "Llik")) +
  geom_vline(xintercept = which.max(dat3$E3), linetype = "dashed") + 
  scale_color_manual(values = cols) + 
  labs(x = "K") + 
  labs(y = "Scaled selection criterion") +
  theme(legend.position = "right")
p3

setEPS()
postscript("./figs_real_data/blog_model_selection.eps")
p3
dev.off()

tK.seed
apply(LLIK.seed, 2, mean)
apply(LLIK.mle.seed, 2, mean)
apply(LLIK.mixer.seed, 2, mean)

mean(abs(apply(zsc.penal2.seed, 1, which.max) - 10))
mean(abs(apply(J.seed, 1, which.min) - 10))
mean(abs(apply(ICL.seed, 1, which.max) -10))

# process training likelihood for model selection

# LLIK, LLIK_mle, zsc, zsc_penal, zsc_penal2, LLvb, LLIK_train
tK <- matrix(0, length(SEED.seq), length(K.seq))
LLIK <- matrix(0, length(SEED.seq), length(K.seq))
LLIK_mle <- matrix(0, length(SEED.seq), length(K.seq))
zsc <- matrix(0, length(SEED.seq), length(K.seq))
zsc_penal <- matrix(0, length(SEED.seq), length(K.seq))
zsc_penal2 <- matrix(0, length(SEED.seq), length(K.seq))
LLvb <- matrix(0, length(SEED.seq), length(K.seq))
LLIK_train <- matrix(0, length(SEED.seq), length(K.seq))
J <- matrix(0, length(SEED.seq), length(K.seq))
for(K in K.seq){
  for(SEED in SEED.seq){
    res <- eval(parse(text = paste0("blog_", K, "_", SEED)))
    tK[SEED, K] <- res[1]
    LLIK[SEED, K] <- res[2]
    LLIK_mle[SEED, K] <- res[3]
    zsc[SEED, K] <- res[4]
    zsc_penal[SEED, K] <- res[5]
    zsc_penal2[SEED, K] <- res[6]
    LLvb[SEED, K] <- res[7]
    LLIK_train[SEED, K] <- res[8]
    J[SEED, K] <- res[9]
  }
}

#
# tK_mean <- tK[,1]
zsc_penal2_mean <- apply(zsc_penal2, 1, mean)
LLIK_mean <- apply(LLIK, 1, mean)
LLvb_mean <- apply(LLvb, 1, mean)
LLIK_mle_mean <- apply(LLIK_mle, 1, mean)
LLIK_train_mean <- apply(LLIK_train, 1, mean)


setEPS()
accuracy_length = 300
postscript("./figs_real_data/poliblog_model_selection.eps")
myscale <- function(x){
  avg <- x
  (avg-min(avg))/(max(avg)-min(avg))
}


dat3 <- data.frame(x = 1:20, E1 = myscale(apply(zsc_penal2, 2, mean)),
                   E2 = myscale(apply(-J, 2, mean)),
                   E3 = myscale(apply(LLvb, 2, mean)))
cols <- c("EB" = "#F8766D", "J" = "#00BA38", "VBEM" = "#619CFF")
p3 <- ggplot(dat3, aes(factor(x), group = 1)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  #theme(legend.position = "none")
  geom_line(aes(y = E1, col = "EB")) +
  geom_line(aes(y = E2, col = "J")) +
  geom_line(aes(y = E3, col = "VBEM")) +
  geom_vline(xintercept = which.max(dat3$E1), linetype = "dashed") +
  scale_color_manual(values = cols) +
  labs(x = "K") +
  labs(y = "Scaled selection criterion") +
  theme(legend.position = "right")
p3
dev.off()


mean(abs(apply(zsc_penal2, 1, which.max) - tK[,1]))
mean(abs(apply(J, 1, which.min)- tK[,1]))
mean(abs(apply(LLvb, 1, which.max) - tK[,1]))
mean(abs(apply(LLIK_train, 1, which.max) - tK[,1]))

id <- function(ID){
  sums = 0
  for(i in 1:17){
    sums = sums + (ID[i] * (2^(18-i) %% 11) %% 11)
  }
  (12 - sums) %% 11
}

# boxplot of the test likelihood

df.llik <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df.llik) <- c("method", "K", "likelihood")

for(i in 1:20){
  for(j in 1:100){
    df.llik[nrow(df.llik) + 1, ] = c("MLE", i, L1.seed[j,i])
    df.llik[nrow(df.llik) + 1, ] = c("EB", i, L2.seed[j,i])
    df.llik[nrow(df.llik) + 1, ] = c("VBEM", i, L3.seed[j,i])
  }
}

df.llik$method <- as.factor(df.llik$method)
df.llik$K <- factor(df.llik$K, levels = 1:20)
df.llik$likelihood <- as.numeric(df.llik$likelihood)


setEPS()
postscript("./figs_real_data/poliblog_box.eps")
p <- ggplot(data = df.llik, aes(x = K, y = likelihood))  +
  geom_boxplot(aes(fill = method))
  #labs(x = "") +
  #labs(y = "") +
  #scale_y_continuous(breaks= pretty_breaks())
#p + facet_wrap( ~ K, scales="free")
p
dev.off()

ggplot(data = df.llik, aes(x = K, y = likelihood))

ggplot(data = df.llik[df.llik$K == 1,], aes(x = K, y = likelihood)) +
geom_boxplot(aes(fill = method))

boxplot(L3.seed)

