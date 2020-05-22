library(ggpubr)
library(cowplot)

getstat <- function(dat){
  return(list(err = dat$err,
              llik = dat$llik,
              zscs = dat$zsc,
              J = dat$J,
              icl = dat$ICL,
              label = dat$labels,
              cover_rate = dat$coverage,
              eta_mat = dat$eta,
              states = dat$states))
}

E1.seed <- matrix(0, length(SEED.seq), length(K.seq))
E2.seed <- matrix(0, length(SEED.seq), length(K.seq))
E3.seed <- matrix(0, length(SEED.seq), length(K.seq))
zsc.seed <- matrix(0, length(SEED.seq), length(K.seq))
zsc.penal.seed <- matrix(0, length(SEED.seq), length(K.seq))
J.seed <- matrix(0, length(SEED.seq), length(K.seq))
ICL.seed <- matrix(0, length(SEED.seq), length(K.seq))
labels.seed <- array(dim = c(length(K.seq), N, length(SEED.seq)))
cover.seed <- matrix(0, length(SEED.seq), length(K.seq))
cover.mixer.seed <- matrix(0, length(SEED.seq), length(K.seq))

for(SEED in SEED.seq){
  for(K in K.seq){
    res <- getstat(eval(parse(text = paste0("SBM_", N, "_", SEED, "_", K))))
    index <- SEED
    offset <- min(K.seq) - 1
    E1.seed[index, K-offset] <- res$err[1] # MLE
    E2.seed[index, K-offset] <- res$err[2] # EB
    E3.seed[index, K-offset] <- res$err[3] # VBEM
    zsc.seed[index, K-offset] <- res$zsc[1] # Bayes joint likelihood
    zsc.penal.seed[index, K-offset] <- res$zsc[2] # Penalized Bayes likelihood
    J.seed[index, K-offset] <- res$J
    ICL.seed[index, K-offset] <- res$icl
    labels.seed[K-offset, , index] <- res$label
    cover.seed[index, K-offset] <- res$cover_rate[1]
    cover.mixer.seed[index, K-offset] <- res$cover_rate[2]
  }
}

R1 <- (apply(E2.seed,2,mean)/apply(E1.seed,2,mean))*100
R2 <- (apply(E2.seed,2,mean)/apply(E3.seed,2,mean))*100

shows <- c(5,10,15,tK,20)
xlabels <- rep("", 20)
xlabels[shows] <- shows
# setEPS()
# accuracy_length = 400
# postscript(paste0("./figs_SBM/0509/SBM_sparse_", N, "_", tK, "_MSE_ratio.eps"))
dat1 <- data.frame(x = factor(1:20), 
                   I1 = R1,
                   I2 = R2)
cols <- c("EB/MLE" = "#F8766D", "EB/VBEM" = "#619CFF")
p1 <- ggplot(dat1, aes(x, group = 1)) + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +  
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140), limits = c(0, 140)) + 
  scale_x_discrete(breaks = 1:20,
                   labels = xlabels) + 
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 40)) +
  theme(axis.text.x = element_text(face = "italic",
                                   color = ifelse(dat1$x ==tK, "red", "black"),
                                   size = 40)) + 
  theme(axis.ticks.length.x.bottom = unit(0.5, "cm")) +
  theme(axis.ticks.x.bottom = element_line(colour = ifelse(dat1$x ==tK, "red", "black"))) +
  theme(axis.ticks.x.bottom = element_line(size = ifelse(dat1$x %in% c(10,15,shows,20), 2, .5))) +
  theme(axis.ticks.length.y.left = unit(0.5, "cm")) +
  geom_line(aes(y = I1, col = "EB/MLE"), size = 1.8) +
  geom_line(aes(y = I2, col = "EB/VBEM"), size = 1.8) +
  geom_point(aes(y = I1, col = "EB/MLE"), size = 1.8) +
  geom_point(aes(y = I2, col = "EB/VBEM"), size = 1.8) +
  scale_color_manual(values = cols) + 
  labs(x = "") + 
  labs(y = "") +
  theme(plot.margin=unit(c(0,0.5,0,0),"cm"))+
  geom_hline(yintercept=100, linetype="dashed")
p1
# dev.off()


# plot legend
# setEPS()
# postscript("legend_MSE.eps", height = 0.5)
# cols <- c("MLE" = "#F8766D", "VBEM" = "#00BA38", "EB" = "#619CFF")
# pleg <- ggplot(dat2, aes(x, group = 1)) + 
#   theme(legend.title = element_blank()) +
#   theme(legend.position = "bottom") + 
#   theme(legend.key.width = unit(3, "line")) + 
#   theme(legend.key = element_rect(fill = "white")) +
#   theme(legend.text = element_text(margin = margin(r = 50, unit = "pt"))) + 
#   theme(legend.text = element_text(size = 15, face = "italic")) + 
#   geom_line(aes(y = E1, col = "MLE"), size = 1.8) +
#   geom_line(aes(y = E2, col = "VBEM"), size = 1.8) +
#   geom_line(aes(y = E3, col = "EB"), size = 1.8) +
#   scale_color_manual(values = cols) 
# as_ggplot(get_legend(pleg))
# dev.off()

myscale <- function(x){
  avg <- apply(x, 2, mean)
  (avg-min(avg))/(max(avg)-min(avg))
}

sel_ICL <- myscale(ICL.seed)
sel_J <- myscale(J.seed)
sel_EB <- myscale(zsc.penal.seed)

# setEPS()
# accuracy_length = 400
# postscript(paste0("./figs_SBM/0509/SBM_sparse_",  N, "_", tK, "_model_selection.eps"))
dat3 <- data.frame(x = 1:20, E1 = sel_ICL, 
                   E2 = sel_J, 
                   E3 = sel_EB)
cols <- c("Ilvb" = "#F8766D", "J" = "#00BA38", "Llik" = "#619CFF")
p3 <- ggplot(dat3, aes(factor(x), group = 1)) +
  theme_classic() +
  theme(legend.position = "none") + 
  # scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125), limits = c(0,125)) + 
  scale_x_discrete(breaks = 1:20,
                   labels = xlabels) + 
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 40)) +
  theme(axis.text.x = element_text(face = "italic",
                                   color = ifelse(dat1$x ==tK, "red", "black"),
                                   size = 40)) + 
  theme(axis.ticks.length.x.bottom = unit(0.5, "cm")) +
  theme(axis.ticks.x.bottom = element_line(colour = ifelse(dat1$x ==tK, "red", "black"))) +
  theme(axis.ticks.x.bottom = element_line(size = ifelse(dat1$x %in% c(10,15,shows,20), 2, .5))) +
  theme(axis.ticks.length.y.left = unit(0.5, "cm")) +
  geom_line(aes(y = E1, col = "Ilvb"), size = 1.8) +
  geom_line(aes(y = E2, col = "J"), size = 1.8) +
  geom_line(aes(y = E3, col = "Llik"), size = 1.8) +
  geom_vline(xintercept = which.max(dat3$E3), linetype = "dashed") + 
  scale_color_manual(values = cols) + 
  labs(x = "") + 
  labs(y = "") +
  theme(plot.margin=unit(c(0,0.5,0,0),"cm"))
p3
# dev.off()

# setEPS()
# postscript(paste0("./figs_SBM/0509/legend_llik.eps"), height = 0.5)
# cols <- c("VBEM" = "#F8766D", "EB" = "#619CFF", "CVRP" = "#00BA38")
# pleg <- ggplot(dat3, aes(x, group = 1)) +
#   theme(legend.title = element_blank()) +
#   theme(legend.position = "bottom") +
#   theme(legend.key = element_rect(fill = "white")) +
#   theme(legend.key.width = unit(3, "line")) +
#   theme(legend.text = element_text(margin = margin(r = 50, unit = "pt"))) +
#   theme(legend.text = element_text(size = 15)) +
#   geom_line(aes(y = E3, col = "EB"), size = 1.8) +
#   geom_line(aes(y = E1, col = "VBEM"), size = 1.8) +
#   geom_line(aes(y = E2, col = "CVRP"), size = 1.8) +
#   scale_color_manual(values = cols)
# as_ggplot(get_legend(pleg))
# dev.off()

library(xtable)
Krange <- 1:20
EB <- table(factor(apply(zsc.penal.seed,1, which.max), Krange))
J <- table(factor(apply(J.seed, 1, which.min), Krange))
Ilvb <- table(factor(apply(ICL.seed, 1, which.max), Krange))
