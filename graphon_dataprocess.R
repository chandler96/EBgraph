source("../../common_functions.R")

getstat <- function(dat){
  return(list(err = dat$err,
              #mae = dat$mae,
              llik = dat$llik,
              zscs = dat$zsc,
              J = dat$J,
              icl = dat$ICL,
              label = dat$labels,
              #cover_rate = dat$coverage,
              eta_mat = dat$eta,
              states = dat$states))
}

E1.seed <- matrix(0, length(SEED.seq), length(K.seq))
E2.seed <- matrix(0, length(SEED.seq), length(K.seq))
E3.seed <- matrix(0, length(SEED.seq), length(K.seq))
A1.seed <- matrix(0, length(SEED.seq), length(K.seq))
A2.seed <- matrix(0, length(SEED.seq), length(K.seq))
A3.seed <- matrix(0, length(SEED.seq), length(K.seq))
LLIK.seed <- matrix(0, length(SEED.seq), length(K.seq))
zsc.seed <- matrix(0, length(SEED.seq), length(K.seq))
zsc.penal.seed <- matrix(0, length(SEED.seq), length(K.seq))
zsc.penal2.seed <- matrix(0, length(SEED.seq), length(K.seq))
J.seed <- matrix(0, length(SEED.seq), length(K.seq))
ICL.seed <- matrix(0, length(SEED.seq), length(K.seq))
labels.seed <- array(dim = c(length(K.seq), N, length(SEED.seq)))
cover.seed <- matrix(0, length(SEED.seq), length(K.seq))
cover.mixer.seed <- matrix(0, length(SEED.seq), length(K.seq))

for(SEED in SEED.seq){
  for(K in K.seq){
    res <- getstat(eval(parse(text = paste0("graphon_", N, "_", -logrho, "_", lambda, "_", K, "_", SEED))))
    index <- SEED
    offset <- 0
    E1.seed[index, K-offset] <- res$err[1] # MLE
    E2.seed[index, K-offset] <- res$err[2] # mixer
    E3.seed[index, K-offset] <- res$err[3] # EB
    #A1.seed[index, K-offset] <- res$mae[1] # MLE
    #A2.seed[index, K-offset] <- res$mae[2] # mixer
    #A3.seed[index, K-offset] <- res$mae[3] # EB
    LLIK.seed[index, K-offset] <- res$llik[1] # Bayes likelihood
    zsc.seed[index, K-offset] <- res$zsc[1] # Bayes joint likelihood
    zsc.penal.seed[index, K-offset] <- res$zsc[2] # Penalized Bayes likelihood
    zsc.penal2.seed[index, K-offset] <- res$zsc[3] # icl Penalized Bayes likelihood
    J.seed[index, K-offset] <- res$J
    ICL.seed[index, K-offset] <- res$icl
    #labels.seed[K-offset, , index] <- res$label
    #cover.seed[index, K-offset] <- res$cover_rate[1]
    #cover.mixer.seed[index, K-offset] <- res$cover_rate[2]
  }
}

E1.seed <- E1.seed # MLE
temp.seed <- E3.seed
E3.seed <- E2.seed  # VBEM
E2.seed <- temp.seed # EB

# #png("./newestplots/picG4.png")
# theme_update(plot.title = element_text(hjust = 0.5)) 
# dat1 <- data.frame(x = 1:10, 
#                    I1 = (1 - apply(10^E2.seed, 2, mean)/apply(10^E1.seed, 2, mean))*100,
#                    I2 = (1 - apply(10^E2.seed, 2, mean)/apply(10^E3.seed, 2, mean))*100)
# cols <- c("EB/MLE" = "#F8766D", "EB/VBEM" = "#619CFF")
# p1 <- ggplot(dat1, aes(x)) +
#   theme_classic() +
#   theme(legend.title = element_blank(), axis.text.x=element_blank()) +
#   #theme(legend.position = "none")
#   geom_line(aes(y = I1, col = "EB/MLE")) +
#   geom_line(aes(y = I2, col = "EB/VBEM")) + 
#   #geom_vline(xintercept = which.min(dat1$E2), linetype = "dashed") + 
#   scale_color_manual(values = cols) + 
#   labs(x = "") + 
#   labs(y = "mean RMSE reduction") +
#   geom_hline(yintercept=0, linetype="dashed") + 
#   theme(legend.position = "bottom")
# p1
# ## boxplot
# dat2 <- data.frame(x = rep(rep(c(1:10), each = 100), 3), 
#                    y = c(c(E1.seed), c(E2.seed), c(E3.seed)),
#                    gp = rep(1:3, each = 1000))
# p2 <- ggplot(dat2, aes(x = factor(x), y = y)) +
#   theme_classic() + 
#   geom_boxplot(aes(fill = factor(gp)), outlier.shape = NA, show.legend = FALSE) + 
#   labs(x = "K") + 
#   labs(y = "log10(RMSE)") +
#   scale_fill_discrete(name = "", labels = c("MLE", "EB", "VBEM")) +
#   theme(legend.position = "bottom") 
# p2
# grid.arrange(p2, p1, ncol = 1)
# #dev.off()


#R21 <- (mean(E2.seed^2)/mean(E1.seed^2))*100
#R23 <- (mean(E2.seed^2)/mean(E3.seed^2))*100

# m_cal <- function(S){
#   q1 <- quantile(S, 0.05)
#   q3 <- quantile(S, 0.95)
#   S_sel <- S[S >= q1 & S <= q3]
#   return(mean(S_sel))
# }

R21 <- apply(E2.seed^2,2 , mean)/apply(E1.seed^2,2 , mean)*100
R23 <- apply(E2.seed^2,2 , mean)/apply(E3.seed^2,2 , mean)*100

# sd_cal <- function(S){
#   q1 <- quantile(S, 0.025)
#   q3 <- quantile(S, 0.095)
#   S_sel <- S[S >= q1 & S <= q3]
#   return(sd(S_sel))
# }

# 
# sd21 <- apply(R21, 2, sd_cal)
# sd23 <- apply(R23, 2, sd_cal)

#shows <- c(5,10,15,tK,20)
xlabels <- 1:10
#xlabels[shows] <- shows
dat1 <- data.frame(x = factor(1:10), 
                   I1 = R21,
                   I2 = R23)
cols <- c("EB/MLE" = "#F8766D", "EB/VBEM" = "#619CFF")
p1 <- ggplot(dat1, aes(x, group = 1)) + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") + 
  theme(axis.title = element_text(size = 25)) + 
  theme(axis.text = element_text(size = 18)) + 
  theme(axis.text.y = element_text(face = "italic", color = "black")) +
  theme(axis.text.x = element_text(face = "italic",
                                   color = "black")) + 
  geom_line(aes(y = I1, col = "EB/MLE"), size = 1.8) +
  geom_line(aes(y = I2, col = "EB/VBEM"), size = 1.8) +
  geom_point(aes(y = I1, col = "EB/MLE"), size = 1.8) +
  geom_point(aes(y = I2, col = "EB/VBEM"), size = 1.8) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140), limits = c(0, 140)) + 
  scale_x_discrete(breaks = 1:10,
                   labels = xlabels)+
  theme(axis.ticks.length.x.bottom = unit(0.5, "cm")) +
  #theme(axis.ticks.x.bottom = element_line(colour = ifelse(dat1$x ==tK, "red", "black"))) +
  #theme(axis.ticks.x.bottom = element_line(size = ifelse(dat1$x %in% c(10,15,shows,20), 2, .5))) +
  theme(axis.ticks.length.y.left = unit(0.5, "cm")) +
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 40)) +
  theme(axis.text.x = element_text(face = "italic",
                                   color = "black",
                                   size = 40)) + 
  #geom_errorbar(aes(ymin = I1-sd21,
  #                  ymax = I1+sd21, col = "EB/MLE"),
  #                  width = .5, position = position_dodge(0.05),size = 1.2) +
  #geom_errorbar(aes(ymin = I2-sd23,
  #                  ymax = I2+sd23, col = "EB/VBEM"),
  #              width = .5, position = position_dodge(0.05), size = 1.2) +
  scale_color_manual(values = cols) + 
  labs(x = "") + 
  labs(y = "") +
  #ylim(-1, 150) +
  theme(plot.margin=unit(c(0,0.5,0,0),"cm"))+
  geom_hline(yintercept=100, linetype="dashed")
p1

# setEPS()
# accuracy_length = 400
# postscript(paste0("../0509/G_", N, "_", -logrho, "_", lambda, "_MSE_ratio.eps"))
# p1
# dev.off()

# dat2 <- data.frame(x = 1:10, E1 = apply(E1.seed^2, 2, mean), 
#                    E2 = apply(E2.seed^2, 2, mean), 
#                    E3 = apply(E3.seed^2, 2, mean))
# cols <- c("MLE" = "#F8766D", "VBEM" = "#00BA38", "EB" = "#619CFF")
# p2 <- ggplot(dat2, aes(factor(x), group = 1)) +
#   theme_classic() +
#   theme(legend.title = element_blank()) +
#   theme(legend.position = "none") +
#   theme(axis.title = element_text(size = 25)) + 
#   theme(axis.text = element_text(size = 18)) + 
#   theme(axis.text.y = element_text(face = "italic", color = "black")) +
#   theme(axis.text.x = element_text(face = "italic",
#                                    color = "black")) +
#   geom_line(aes(y = E1, col = "MLE"), size = 1.8) +
#   geom_line(aes(y = E2, col = "EB"), size = 1.8) +
#   geom_line(aes(y = E3, col = "VBEM"), size = 1.8) +
#   #geom_vline(xintercept = which.min(dat2$E2), linetype = "dashed") + 
#   scale_color_manual(values = cols) + 
#   labs(x = "") + 
#   labs(y = "MSE")
# 
# setEPS()
# accuracy_length = 400
# postscript(paste0("./figs_graphon/0201/G_", N, "_", -logrho, "_", lambda, "_MSE.eps"))
# p2
# dev.off()


# mean(graph_density)
# 
# graph_density
# print("RMSE improvement")
# (1 - apply(10^E2.seed, 2, mean)/apply(10^E1.seed, 2, mean))*100
# (1 - apply(10^E2.seed, 2, mean)/apply(10^E3.seed, 2, mean))*100
# #print("MAE improvement")
# #(1 - apply(A3.seed, 2, mean)/apply(A1.seed, 2, mean))*100
# # (1 - apply(A3.seed, 2, mean)/apply(A2.seed, 2, mean))*100

# print("Estimate Coverage Rate")
# apply(cover.seed, 2, mean)
# apply(cover.mixer.seed, 2, mean)

# model selection
# table(apply(zsc.penal2.seed,1, which.max))
# table(apply(zsc.penal.seed,1, which.max))
# table(apply(J.seed,1, which.max))
# table(apply(ICL.seed,1, which.max))

# model selection error
# e1 <- mean(abs(apply(zsc.penal2.seed,1, which.max) - apply(pmin(E1.seed, E2.seed, E3.seed), 1, which.min)))
# e2 <- mean(abs(apply(zsc.penal.seed,1, which.max) - apply(pmin(E1.seed, E2.seed, E3.seed), 1, which.min)))
# e3 <- mean(abs(apply(J.seed,1, which.max) - apply(pmin(E1.seed, E2.seed, E3.seed), 1, which.min)))
# e4 <- mean(abs(apply(ICL.seed,1, which.max) - apply(pmin(E1.seed, E2.seed, E3.seed), 1, which.min)))


myscale <- function(x){
  avg <- apply(x, 2, mean)
  (avg-min(avg))/(max(avg)-min(avg))
}

sel_ICL <- myscale(ICL.seed)
sel_J <- myscale(J.seed)
sel_llk2 <- myscale(zsc.penal2.seed)

dat3 <- data.frame(x = 1:10, E1 = sel_ICL, 
                   E2 = sel_J, 
                   E3 = sel_llk2)
cols <- c("Ilvb" = "#F8766D", "J" = "#00BA38", "Llik" = "#619CFF")
p3 <- ggplot(dat3, aes(factor(x), group = 1)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 25)) + 
  theme(axis.text = element_text(size = 18)) + 
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 40)) +
  theme(axis.text.x = element_text(face = "italic",
                                   color = "black", size = 40)) +
  theme(axis.ticks.length.x.bottom = unit(0.5, "cm")) +
  #theme(axis.ticks.x.bottom = element_line(colour = ifelse(dat1$x ==tK, "red", "black"))) +
  #theme(axis.ticks.x.bottom = element_line(size = ifelse(dat1$x %in% c(10,15,shows,20), 2, .5))) +
  geom_line(aes(y = E1, col = "Ilvb"), size = 1.8) +
  geom_line(aes(y = E2, col = "J"), size = 1.8) +
  geom_line(aes(y = E3, col = "Llik"), size = 1.8) +
  geom_vline(xintercept = which.max(dat3$E3), linetype = "dashed") + 
  scale_color_manual(values = cols) + 
  labs(x = "") + 
  labs(y = "") +
  theme(plot.margin=unit(c(0,0.5,0,0),"cm"))
p3

# setEPS()
# accuracy_length = 400
# postscript(paste0("../0509/G_", N, "_", -logrho, "_", lambda, "_model_selection.eps"))
# p3
# dev.off()

# 
# 
# png("./newestplots/picP4.png")
# rho = 10^-1
# lambda = 5
# 
# W <- function(x, y) {
#   exp(-rho*((x-0.5)^2 + (y-0.5)^2))
# }
# 
# x <- seq(0, 1, length = 100)
# y <- x
# z <- outer(x, y, W)
# zfacet <- z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)]
# cols <- colorRampPalette( c("white", "red") )
# nbcol <- 100
# cols <- cols(nbcol)
# facetcol <- cut(zfacet, nbcol)
# #persp(x, y, z, xlab = "x", ylab = "y", zlab = "W(x,y)",
# #      main = "C=5", col = cols[facetcol], theta = 90, phi = 90, border = NA)
# persp(x, y, z, xlab = "x", ylab = "y", zlab = "W(x,y)",
#       main = "", col = cols[facetcol], theta = -30, phi = 30, border = NA)
# #dev.off()


# library(xtable)
# Krange <- 1:20
# penal2 <- table(factor(apply(zsc.penal2.seed,1, which.max), Krange))
# penal <- table(factor(apply(zsc.penal.seed, 1, which.max), Krange))
# J <- table(factor(apply(J.seed,1, which.max), Krange))
# Ilvb <- table(factor(apply(ICL.seed,1, which.max), Krange))
# J
# Ilvb
# penal2
# 
# n1<-mean(abs(apply(zsc.penal2.seed,1, which.max) - apply(E1.seed, 1 , which.min)))
# n2<-mean(abs(apply(J.seed,1, which.max) - apply(E1.seed, 1 , which.min)))
# n3<-mean(abs(apply(ICL.seed,1, which.max) - apply(E1.seed, 1 , which.min)))
# #mean(abs(apply(zsc.penal2.seed,1, which.max) - apply(pmin(E1.seed), 1 , which.min)))
# #mean(abs(apply(J.seed,1, which.max) - apply(pmin(E1.seed), 1 , which.min)))
# #mean(abs(apply(ICL.seed,1, which.max) - apply(pmin(E1.seed), 1 , which.min)))
# c(n2,n3,n1)

#units <- 1e5
#paste(round(apply(E1.seed^2, 2, mean)*units), collapse = "&")
#paste(round(apply(E3.seed^2, 2, mean)*units), collapse = "&")
#paste(rep("\bf ", 10), round(apply(E2.seed^2, 2, mean)*units), sep = "", collapse = "&")

# library(latex2exp)
# library(grDevices)
# windowsFonts(A = windowsFont("Helvetica"))
# setEPS()
# postscript(paste0("../0509/texts.eps"))
# accuracy_length = 10
# plot(NA, xlim = c(0,0.1), ylim = c(0,0.1), bty = 'n', xaxt = "n", yaxt = "n",
#      xlab = "", ylab = "", family ="A")
# text(0.05,0.05, TeX('$\\rho = 10^{-1}$, $\\lambda = 2$'),
#      col = "red")
# dev.off()

mean(abs(apply(J.seed, 1, which.max) - apply(E1.seed,1,which.min)))
mean(abs(apply(ICL.seed, 1, which.max) - apply(E1.seed,1,which.min)))
mean(abs(apply(zsc.penal2.seed, 1, which.max) - apply(E1.seed,1,which.min)))
mean(abs(apply(zsc.penal.seed, 1, which.max) - apply(E1.seed,1,which.min)))
lambda

