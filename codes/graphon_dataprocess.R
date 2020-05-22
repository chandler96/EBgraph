source("lliks.R")

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
    res <- getstat(eval(parse(text = paste0("graphon_", N, "_", -logrho, "_", lambda, "_", K, "_", SEED))))
    index <- SEED
    offset <- min(K.seq) - 1
    E1.seed[index, K-offset] <- res$err[1] # MLE
    E2.seed[index, K-offset] <- res$err[2] # VBEM
    E3.seed[index, K-offset] <- res$err[3] # EB
    zsc.seed[index, K-offset] <- res$zsc[1] # Bayes joint likelihood
    zsc.penal.seed[index, K-offset] <- res$zsc[2] # Penalized Bayes likelihood
    J.seed[index, K-offset] <- res$J
    ICL.seed[index, K-offset] <- res$icl
    labels.seed[K-offset, , index] <- res$label
    cover.seed[index, K-offset] <- res$cover_rate[1]
    cover.mixer.seed[index, K-offset] <- res$cover_rate[2]
  }
}

R1 <- apply(E3.seed, 2, mean)/apply(E1.seed, 2, mean)*100
R2 <- apply(E3.seed, 2, mean)/apply(E2.seed, 2, mean)*100

xlabels <- 1:10
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
  theme(axis.ticks.length.y.left = unit(0.5, "cm")) +
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 40)) +
  theme(axis.text.x = element_text(face = "italic",
                                   color = "black",
                                   size = 40)) + 
  scale_color_manual(values = cols) + 
  labs(x = "") + 
  labs(y = "") +
  theme(plot.margin=unit(c(0,0.5,0,0),"cm"))+
  geom_hline(yintercept=100, linetype="dashed")
p1

# dat2 <- data.frame(x = 1:10, E1 = apply(E1.seed, 2, mean), 
#                    E2 = apply(E2.seed, 2, mean), 
#                    E3 = apply(E3.seed, 2, mean))
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
#   geom_line(aes(y = E2, col = "VBEM"), size = 1.8) +
#   geom_line(aes(y = E3, col = "EB"), size = 1.8) +
#   #geom_vline(xintercept = which.min(dat2$E2), linetype = "dashed") + 
#   scale_color_manual(values = cols) + 
#   labs(x = "") + 
#   labs(y = "MSE")

# model selection
table(apply(zsc.penal.seed, 1, which.max))
table(apply(J.seed, 1, which.max))
table(apply(ICL.seed, 1, which.max))

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

### contour plot ###

# setEPS()
# accuracy_length = 400
# postscript(paste0("G_", N, "_", -logrho, "_", lambda, "_model_selection.eps"))
# p3
# dev.off()

# png("pic.png")
# rho = 10^-1
# lambda = 5

# W <- function(x, y) {
#   exp(-rho*((x-0.5)^2 + (y-0.5)^2))
# }

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

mean(abs(apply(J.seed, 1, which.min) - apply(E1.seed, 1, which.min)))
mean(abs(apply(ICL.seed, 1, which.max) - apply(E1.seed, 1, which.min)))
mean(abs(apply(zsc.penal.seed, 1, which.max) - apply(E1.seed, 1, which.min)))

