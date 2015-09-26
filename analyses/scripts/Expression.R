#This script analyzes gene expression vs. CpG O/E

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Apalm")

Apalm_exp_CpG <- read.delim("Apalm_exp_cpg2", header = F)

Apalm_mean_exp <- rowMeans(Apalm_exp_CpG[,2:37])

Apalm_exp_CpG2 <- cbind(Apalm_exp_CpG, Apalm_mean_exp)

Apalm_exp_CpG2_sub <- subset(Apalm_exp_CpG2, V38 >= 0.001 & V38 <= 1.5)

x <- Apalm_exp_CpG2_sub$V38

y <- Apalm_exp_CpG2_sub$Apalm_mean_exp

log_y <- log(y)

plot(x, y)

smoothScatter(x, y, nrpoints = 0, ylim = c(9, 9.5))

lm.1 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7))

lm.s <- step(lm.1)

lm1 <- lm(y ~x)
abline(lm1)
lm2 <- lm(y ~ x + I(x^2))
lines(sort(x), fitted(lm2)[order(x)], col='red', type='l') 
lm3 <- lm(y ~ x + I(x^2) + I(x^3))
lines(sort(x), fitted(lm3)[order(x)], col='orange', type='l') 
lm4 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
lines(sort(x), fitted(lm3)[order(x)], col='purple', type='l') 

x_cut <- cut_number(x, n = 30)

z <- cbind(x_cut, x)

cuts <- z[,1]

Apalm_y_means <- tapply(y, INDEX = cuts, FUN = mean)

Apalm_y_stdev <- tapply(y, INDEX = cuts, FUN = sd)

Apalm_y_stderr <- Apalm_y_stdev/sqrt(4331)

Apalm_x_means <- tapply(x, INDEX = cuts, FUN = mean)

plot(Apalm_x_means, Apalm_y_means)
segments(Apalm_x_means, Apalm_y_means-Apalm_y_stderr, Apalm_x_means, Apalm_y_means+Apalm_y_stderr)

Apalm_lm.1 <- lm(Apalm_y_means ~ Apalm_x_means + I(Apalm_x_means^2) + I(Apalm_x_means^3) + I(Apalm_x_means^4) + I(Apalm_x_means^5) + I(Apalm_x_means^6))

Apalm_lm.s <- step(lm.1)

#Best fit model appears to be 4th order polynomial

Apalm_lm1 <- lm(Apalm_y_means ~ Apalm_x_means)
abline(Apalm_lm1)
Apalm_lm2 <- lm(Apalm_y_means ~ Apalm_x_means + I(Apalm_x_means^2))
lines(sort(Apalm_x_means), fitted(Apalm_lm2)[order(Apalm_x_means)], col='red', type='l') 
Apalm_lm3 <- lm(Apalm_y_means ~ Apalm_x_means + I(Apalm_x_means^2) + I(Apalm_x_means^3))
lines(sort(Apalm_x_means), fitted(Apalm_lm3)[order(Apalm_x_means)], col='orange', type='l') 
Apalm_lm4 <- lm(Apalm_y_means ~ Apalm_x_means + I(Apalm_x_means^2) + I(Apalm_x_means^3) + I(Apalm_x_means^4))
lines(sort(Apalm_x_means), fitted(Apalm_lm4)[order(Apalm_x_means)], col='purple', type='l') 

#Now onto Acropora millepora data

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Amil")

Amil_exp_CpG <- read.delim("Amil_exp_CpG2", header = F)

cols <- c(3, 5, 7, 9, 11, 13, 15, 17, 19)

Amil_mean_exp <- rowMeans(Amil_exp_CpG[,cols])

Amil_exp_CpG2 <- cbind(Amil_exp_CpG, Amil_mean_exp)

Amil_exp_CpG2_sub <- subset(Amil_exp_CpG2, V20 >= 0.001 & V20 <= 1.5)

x <- Amil_exp_CpG2_sub$V20

y <- Amil_exp_CpG2_sub$Amil_mean_exp

log_y <- log(y)

plot(x, log_y) #Using log transform here for y

smoothScatter(x, log_y, nrpoints = 0)

lm.1 <- lm(log_y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7))

lm.s <- step(lm.1)

lm1 <- lm(y ~x)
abline(lm1)
lm2 <- lm(y ~ x + I(x^2))
lines(sort(x), fitted(lm2)[order(x)], col='red', type='l') 
lm3 <- lm(y ~ x + I(x^2) + I(x^3))
lines(sort(x), fitted(lm3)[order(x)], col='orange', type='l') 
lm4 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
lines(sort(x), fitted(lm3)[order(x)], col='purple', type='l') 

x_cut <- cut_number(x, n = 30)

z <- cbind(x_cut, x)

cuts <- z[,1]

Amil_y_means <- tapply(log_y, INDEX = cuts, FUN = mean) #Using log transform here for y

Amil_y_stdev <- tapply(log_y, INDEX = cuts, FUN = sd)

Amil_y_stderr <- y_stdev/sqrt(1638) 

Amil_x_means <- tapply(x, INDEX = cuts, FUN = mean)

plot(Amil_x_means, Amil_y_means)
segments(Amil_x_means, Amil_y_means-Amil_y_stderr, Amil_x_means, Amil_y_means+Amil_y_stderr)

Amil_lm.1 <- lm(Amil_y_means ~ Amil_x_means + I(Amil_x_means^2) + I(Amil_x_means^3) + I(Amil_x_means^4) + I(Amil_x_means^5) + I(Amil_x_means^6))

Amil_lm.s <- step(Amil_lm.1)

#Best fit model appears to be 4th order polynomial

Amil_lm1 <- lm(Amil_y_means ~ Amil_x_means)
abline(lm1)
Amil_lm2 <- lm(Amil_y_means ~ Amil_x_means + I(Amil_x_means^2))
lines(sort(Amil_x_means), fitted(Amil_lm2)[order(Amil_x_means)], col='red', type='l') 
Amil_lm3 <- lm(Amil_y_means ~ Amil_x_means + I(Amil_x_means^2) + I(Amil_x_means^3))
lines(sort(Amil_x_means), fitted(Amil_lm3)[order(Amil_x_means)], col='orange', type='l') 
Amil_lm4 <- lm(Amil_y_means ~ Amil_x_means + I(Amil_x_means^2) + I(Amil_x_means^3) + I(Amil_x_means^4))
lines(sort(Amil_x_means), fitted(Amil_lm4)[order(Amil_x_means)], col='purple', type='l') 

#Now onto Acropora hyacinthus data

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Ahya")

Ahya_exp_CpG <- read.delim("Ahya_exp_CpG2", header = F)

Ahya_mean_exp <- rowMeans(Ahya_exp_CpG[,2:33])

Ahya_exp_CpG2 <- cbind(Ahya_exp_CpG, Ahya_mean_exp)

Ahya_exp_CpG2_sub <- subset(Ahya_exp_CpG2, V34 >= 0.001 & V34 <= 1.5)

x <- Ahya_exp_CpG2_sub$V34

y <- Ahya_exp_CpG2_sub$Ahya_mean_exp

log_y <- log(y)

plot(x, log_y)

smoothScatter(x, log_y, nrpoints = 0)

lm.1 <- lm(log_y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7))

lm.s <- step(lm.1)

lm1 <- lm(log_y ~x)
abline(lm1)
lm2 <- lm(log_y ~ x + I(x^2))
lines(sort(x), fitted(lm2)[order(x)], col='red', type='l') 
lm3 <- lm(log_y ~ x + I(x^2) + I(x^3))
lines(sort(x), fitted(lm3)[order(x)], col='orange', type='l') 
lm4 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
lines(sort(x), fitted(lm3)[order(x)], col='purple', type='l') 

x_cut <- cut_number(x, n = 30)

z <- cbind(x_cut, x)

cuts <- z[,1]

Ahya_y_means <- tapply(log_y, INDEX = cuts, FUN = mean) #Using log transform of y

Ahya_y_stdev <- tapply(log_y, INDEX = cuts, FUN = sd)

Ahya_y_stderr <- y_stdev/sqrt(1064)

Ahya_x_means <- tapply(x, INDEX = cuts, FUN = mean)

plot(Ahya_x_means, Ahya_y_means)
segments(Ahya_x_means, Ahya_y_means-Ahya_y_stderr, Ahya_x_means, Ahya_y_means+Ahya_y_stderr)

Ahya_lm.1 <- lm(Ahya_y_means ~ Ahya_x_means + I(Ahya_x_means^2) + I(Ahya_x_means^3) + I(Ahya_x_means^4) + I(Ahya_x_means^5) + I(Ahya_x_means^6))

Ahya_lm.s <- step(Ahya_lm.1)

Ahya_lm1 <- lm(Ahya_y_means ~ Ahya_x_means)

#Best fit model according to AIC is 2nd order polynomial

abline(lm1)
Ahya_lm2 <- lm(Ahya_y_means ~ Ahya_x_means + I(Ahya_x_means^2))
lines(sort(Ahya_x_means), fitted(Ahya_lm2)[order(Ahya_x_means)], col='red', type='l') 
Ahya_lm3 <- lm(Ahya_y_means ~ Ahya_x_means + I(Ahya_x_means^2) + I(Ahya_x_means^3))
lines(sort(Ahya_x_means), fitted(Ahya_lm3)[order(Ahya_x_means)], col='orange', type='l') 
Ahya_lm4 <- lm(Ahya_y_means ~ Ahya_x_means + I(Ahya_x_means^2) + I(Ahya_x_means^3) + I(Ahya_x_means^4))
lines(sort(Ahya_x_means), fitted(Ahya_lm4)[order(Ahya_x_means)], col='purple', type='l') 
Ahya_lm5 <- lm(Ahya_y_means ~ Ahya_x_means + I(Ahya_x_means^2) + I(Ahya_x_means^3) + I(Ahya_x_means^4) + I(Ahya_x_means^5))
lines(sort(Ahya_x_means), fitted(Ahya_lm5)[order(Ahya_x_means)], col='blue', type='l') 
Ahya_lm6 <- lm(Ahya_y_means ~ Ahya_x_means + I(Ahya_x_means^2) + I(Ahya_x_means^3) + I(Ahya_x_means^4) + I(Ahya_x_means^5) + I(Ahya_x_means^6))
lines(sort(Ahya_x_means), fitted(Ahya_lm6)[order(Ahya_x_means)], col='green', type='l') 

#Now plotting all together

par(mfrow = c(1,3)) # 1 x 3 plots
par(xaxs="i", yaxs="i") 
par(mar = c(5, 5, 4, 1), oma = c(1, 2, 1, 1))

plot(Ahya_x_means, Ahya_y_means, ylab = "Mean expression level", xlab = " ", xlim = c(0,1.5), ylim = c(1,2), xaxp  = c(0, 1.5, 3), yaxp  = c(1, 2, 4), las = 1, ann = FALSE)
segments(Ahya_x_means, Ahya_y_means-Ahya_y_stderr, Ahya_x_means, Ahya_y_means+Ahya_y_stderr)
lines(sort(Ahya_x_means), fitted(Ahya_lm2)[order(Ahya_x_means)], col='red', type='l')
mtext(side = 2, text = "Mean expression level", line = 4, cex = 0.8)
mtext(side = 3, text = "Acropora hyacinthus", line = 2, cex = 0.8,  font  = 3)

plot(Amil_x_means, Amil_y_means, ylab = " ", xlab = "CpG O/E", xlim = c(0,1.5), ylim = c(3,5), xaxp  = c(0, 1.5, 3), yaxp  = c(3, 5, 4), las = 1)
segments(Amil_x_means, Amil_y_means-Amil_y_stderr, Amil_x_means, Amil_y_means+Amil_y_stderr)
lines(sort(Amil_x_means), fitted(Amil_lm4)[order(Amil_x_means)], col='red', type='l')
mtext(side = 3, text = "Acropora millepora", line = 2, cex = 0.8,  font  = 3)

plot(Apalm_x_means, Apalm_y_means, ylab = " ", xlab = " ", xlim = c(0,1.5), ylim = c(9.1,9.35), xaxp  = c(0, 1.5, 3), yaxp  = c(9.1, 9.35, 5), las = 1)
segments(Apalm_x_means, Apalm_y_means-Apalm_y_stderr, Apalm_x_means, Apalm_y_means+Apalm_y_stderr)
lines(sort(Apalm_x_means), fitted(Apalm_lm4)[order(Apalm_x_means)], col='red', type='l') 
mtext(side = 3, text = "Acropora palmata", line = 2, cex = 0.8,  font  = 3)
