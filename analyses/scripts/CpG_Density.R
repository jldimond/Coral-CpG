# This script plots and compares mixture models for CpG O/E data.
# Source files are derived from analyses in jupyter (.ipynb) notebooks.

#Read in CpG data and remove NA values

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Ahya")

Ahya<-read.delim("Ahya_cpg_anno", header=F)

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Amil")

Amil<-read.delim("Amil_cpg_anno", header=F)

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Apalm")

Apalm<-read.delim("Apalm_cpg_anno", header=F)

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Pdam")

Pdam<-read.delim("Pdam_cpg_anno", header=F)

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Spist")

Spist<-read.delim("Spist_cpg_anno", header=F)

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Past")

Past<-read.delim("Past_cpg_anno", header=F)

#Fitting mixture model with mixtools normalmixEM

library(mixtools)

par(mfrow = c(2, 3)) # 2 x 3 plots

data <- Ahya$V2
Ahya_data <- data[data >= 0.001 & data <= 1.5] #Cutting off high and low values
set.seed(1)
Ahya_mixmdl <- normalmixEM(Ahya_data, k = 2)
plot(Ahya_mixmdl, which = 2, col2 = c("red", "blue"), xlab2 = " ", main2 = "Acropora hyacinthus", font.main = 3)
lines(density(Ahya_data), lty=2, lwd=2)

data <- Amil$V2
Amil_data <- data[data >= 0.001 & data <= 1.5] #Cutting off high and low values
set.seed(1)
Amil_mixmdl <- normalmixEM(Amil_data)
plot(Amil_mixmdl, which = 2, col2 = c("red", "blue"), xlab2 = " ", main2 = "Acropora millepora", font.main = 3)
lines(density(Amil_data), lty=2, lwd=2)

data <- Apalm$V2
Apalm_data <- data[data >= 0.001 & data <= 1.5] #Cutting off high and low values
set.seed(101)
Apalm_mixmdl <- normalmixEM(Apalm_data)
plot(Apalm_mixmdl, which = 2, col2 = c("red", "blue"), xlab2 = " ", main2 = "Acropora palmata", font.main = 3)
lines(density(Apalm_data), lty=2, lwd=2)

data <- Pdam$V2
Pdam_data <- data[data >= 0.001 & data <= 1.5] #Cutting off high and low values
set.seed(100)
Pdam_mixmdl <- normalmixEM(Pdam_data)
plot(Pdam_mixmdl, which = 2, col2 = c("red", "blue"), xlab2 = " ", main2 = "Pocillopora damicornis", font.main = 3)
lines(density(Pdam_data), lty=2, lwd=2)

data <- Past$V2
Past_data <- data[data >= 0.001 & data <= 1.5] #Cutting off high and low values
set.seed(101)
Past_mixmdl <- normalmixEM(Past_data)
plot(Past_mixmdl, which = 2, col2 = c("red", "blue"), xlab2 = "CpG O/E", main2 = "Porites astreoides", font.main = 3)
lines(density(Past_data), lty=2, lwd=2)

data <- Spist$V2
Spist_data <- data[data >= 0.001 & data <= 1.5] #Cutting off high and low values
set.seed(101)
Spist_mixmdl <- normalmixEM(Spist_data)
plot(Spist_mixmdl, which = 2, col2 = c("red", "blue"), xlab2 = " ", main2 = "Stylophora pistillata", font.main = 3)
lines(density(Spist_data), lty=2, lwd=2)

##### Finds intersection point of two component model

intersect <- function(m1, s1, m2, s2, prop1, prop2){
  
  B <- (m1/s1^2 - m2/s2^2)
  A <- 0.5*(1/s2^2 - 1/s1^2)
  C <- 0.5*(m2^2/s2^2 - m1^2/s1^2) - log((s1/s2)*(prop2/prop1))
  
  (-B + c(1,-1)*sqrt(B^2 - 4*A*C))/(2*A)
}

Ahya_intersect <- intersect(0.30, 0.11, 0.75, 0.26, 0.14, 0.86)
Amil_intersect <- intersect(0.38, 0.11, 0.74, 0.20, 0.17, 0.83)
Apalm_intersect <- intersect(0.35, 0.12, 0.73, 0.22, 0.27, 0.73)
Pdam_intersect <- intersect(0.24, 0.10, 0.69, 0.25, 0.22, 0.78)
Past_intersect <- intersect(0.34, 0.12, 0.73, 0.23, 0.26, 0.74)
Spist_intersect <- intersect(0.33,0.11, 0.72, 0.23, 0.28, 0.72)

##### Tests fit of single component model

Ahya_1comp <- fitdistr(Ahya_data, "normal")
Amil_1comp <- fitdistr(Amil_data, "normal")
Apalm_1comp <- fitdistr(Apalm_data, "normal")
Pdam_1comp <- fitdistr(Pdam_data, "normal")
Past_1comp <- fitdistr(Past_data, "normal")
Spist_1comp <- fitdistr(Spist_data, "normal")
