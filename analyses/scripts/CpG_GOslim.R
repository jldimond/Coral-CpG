# This script deals with reading in a tab delimited file
# containing CpGo/e and GOSlim information.
# The script performs Fisher's exact tests and plots various types
# of figures. Note that some files are derived from analyses created in 
# CpG_Density.R, so that script should be run prior to and alongside this
# one.
# Source files are derived from analyses in jupyter (.ipynb) notebooks.

#Read in CpG data, remove NA values, and compute means and SEs

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Ahya")

Ahya<-na.omit(read.delim("Ahya_cpg_GOslim.tab", header=F))
Ahya2<-subset(Ahya, V2 >= 0.001 & V2 <= 1.5)
AhyaMean<-tapply(Ahya2$V2, Ahya2$V3, mean)
AhyaSD<-tapply(Ahya2$V2, Ahya2$V3, sd)
AhyaLength<-tapply(Ahya2$V2, Ahya2$V3, length)
AhyaSE<-AhyaSD/sqrt(AhyaLength)
AhyaMean<-AhyaMean[2:15]
AhyaSE<-AhyaSE[2:15]

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Amil")

Amil<-na.omit(read.delim("Amil_cpg_GOslim.tab", header=F))
Amil2<-subset(Amil, V2 >= 0.001 & V2 <= 1.5)
AmilMean<-tapply(Amil2$V2, Amil2$V3, mean)
AmilSD<-tapply(Amil2$V2, Amil2$V3, sd)
AmilLength<-tapply(Amil2$V2, Amil2$V3, length)
AmilSE<-AmilSD/sqrt(AmilLength)
AmilMean<-AmilMean[2:15]
AmilSE<-AmilSE[2:15]

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Apalm")

Apalm<-na.omit(read.delim("Apalm_cpg_GOslim.tab", header=F))
Apalm2<-subset(Apalm, V2 >= 0.001 & V2 <= 1.5)
ApalmMean<-tapply(Apalm2$V2, Apalm2$V3, mean)
ApalmSD<-tapply(Apalm2$V2, Apalm2$V3, sd)
ApalmLength<-tapply(Apalm2$V2, Apalm2$V3, length)
ApalmSE<-ApalmSD/sqrt(ApalmLength)
ApalmMean<-ApalmMean[2:15]
ApalmSE<-ApalmSE[2:15]

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Pdam")

Pdam<-na.omit(read.delim("Pdam_cpg_GOslim.tab", header=F))
Pdam2<-subset(Pdam, V2 >= 0.001 & V2 <= 1.5)
PdamMean<-tapply(Pdam2$V2, Pdam2$V3, mean)
PdamSD<-tapply(Pdam2$V2, Pdam2$V3, sd)
PdamLength<-tapply(Pdam2$V2, Pdam2$V3, length)
PdamSE<-PdamSD/sqrt(PdamLength)
PdamMean<-PdamMean[2:15]
PdamSE<-PdamSE[2:15]

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Spist")

Spist<-na.omit(read.delim("Spist_cpg_GOslim.tab", header=F))
Spist2<-subset(Spist, V2 >= 0.001 & V2 <= 1.5)
SpistMean<-tapply(Spist2$V2, Spist2$V3, mean)
SpistSD<-tapply(Spist2$V2, Spist2$V3, sd)
SpistLength<-tapply(Spist2$V2, Spist2$V3, length)
SpistSE<-SpistSD/sqrt(SpistLength)

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses/Past")

Past<-na.omit(read.delim("Past_cpg_GOslim.tab", header=F))
Past2<-subset(Past, V2 >= 0.001 & V2 <= 1.5)
PastMean<-tapply(Past2$V2, Past2$V3, mean)
PastSD<-tapply(Past2$V2, Past2$V3, sd)
PastLength<-tapply(Past2$V2, Past2$V3, length)
PastSE<-PastSD/sqrt(PastLength)
PastMean<-PastMean[2:15]
PastSE<-PastSE[2:15]

#Classify data as high- or low- CpG O/E based on intersection of components
#in mixture model, then apply Fisher's exact test.
#Intersection point file derived from another script (CpG_Density.R)

Ahya_class <- ifelse(Ahya2$V2 > Ahya_intersect[2] ,"high" ,"low")
Ahya_cpg_class <- cbind(Ahya2, Ahya_class)
Ahya_cpg_class2 <- subset(Ahya_cpg_class, select=c("V3", "Ahya_class"))
Ahya_class_table <- table(Ahya_cpg_class2)
Ahya_class_sum <- table(Ahya_class)
Ahya_table_matrix <- as.matrix(Ahya_class_table[1:15,1:2])
Ahya_high <- rep(c(Ahya_class_sum[1]), 15)
Ahya_low <- rep(c(Ahya_class_sum[2]), 15)
Ahya_all <- cbind(Ahya_table_matrix, Ahya_high, Ahya_low)
Ahya_fisher <- as.vector(apply(Ahya_all,1, function(x) fisher.test(matrix(x,nr=2))$p.value))

Amil_class <- ifelse(Amil2$V2 > Amil_intersect[2] ,"high" ,"low")
Amil_cpg_class <- cbind(Amil2, Amil_class)
Amil_cpg_class2 <- subset(Amil_cpg_class, select=c("V3", "Amil_class"))
Amil_class_table <- table(Amil_cpg_class2)
Amil_class_sum <- table(Amil_class)
Amil_table_matrix <- as.matrix(Amil_class_table[1:15,1:2])
Amil_high <- rep(c(Amil_class_sum[1]), 15)
Amil_low <- rep(c(Amil_class_sum[2]), 15)
Amil_all <- cbind(Amil_table_matrix, Amil_high, Amil_low)
Amil_fisher <- as.vector(apply(Amil_all,1, function(x) fisher.test(matrix(x,nr=2))$p.value))

Apalm_class <- ifelse(Apalm2$V2 > Apalm_intersect[2] ,"high" ,"low")
Apalm_cpg_class <- cbind(Apalm2, Apalm_class)
Apalm_cpg_class2 <- subset(Apalm_cpg_class, select=c("V3", "Apalm_class"))
Apalm_class_table <- table(Apalm_cpg_class2)
Apalm_class_sum <- table(Apalm_class)
Apalm_table_matrix <- as.matrix(Apalm_class_table[1:15,1:2])
Apalm_high <- rep(c(Apalm_class_sum[1]), 15)
Apalm_low <- rep(c(Apalm_class_sum[2]), 15)
Apalm_all <- cbind(Apalm_table_matrix, Apalm_high, Apalm_low)
Apalm_fisher <- as.vector(apply(Apalm_all,1, function(x) fisher.test(matrix(x,nr=2))$p.value))

Pdam_class <- ifelse(Pdam2$V2 > Pdam_intersect[2] ,"high" ,"low")
Pdam_cpg_class <- cbind(Pdam2, Pdam_class)
Pdam_cpg_class2 <- subset(Pdam_cpg_class, select=c("V3", "Pdam_class"))
Pdam_class_table <- table(Pdam_cpg_class2)
Pdam_class_sum <- table(Pdam_class)
Pdam_table_matrix <- as.matrix(Pdam_class_table[1:15,1:2])
Pdam_high <- rep(c(Pdam_class_sum[1]), 15)
Pdam_low <- rep(c(Pdam_class_sum[2]), 15)
Pdam_all <- cbind(Pdam_table_matrix, Pdam_high, Pdam_low)
Pdam_fisher <- as.vector(apply(Pdam_all,1, function(x) fisher.test(matrix(x,nr=2))$p.value))

Past_class <- ifelse(Past2$V2 > Past_intersect[2] ,"high" ,"low")
Past_cpg_class <- cbind(Past2, Past_class)
Past_cpg_class2 <- subset(Past_cpg_class, select=c("V3", "Past_class"))
Past_class_table <- table(Past_cpg_class2)
Past_class_sum <- table(Past_class)
Past_table_matrix <- as.matrix(Past_class_table[1:15,1:2])
Past_high <- rep(c(Past_class_sum[1]), 15)
Past_low <- rep(c(Past_class_sum[2]), 15)
Past_all <- cbind(Past_table_matrix, Past_high, Past_low)
Past_fisher <- as.vector(apply(Past_all,1, function(x) fisher.test(matrix(x,nr=2))$p.value))

Spist_class <- ifelse(Spist2$V2 > Spist_intersect[2] ,"high" ,"low")
Spist_cpg_class <- cbind(Spist2, Spist_class)
Spist_cpg_class2 <- subset(Spist_cpg_class, select=c("V3", "Spist_class"))
Spist_class_table <- table(Spist_cpg_class2)
Spist_class_sum <- table(Spist_class)
Spist_table_matrix <- as.matrix(Spist_class_table[1:14,1:2])
Spist_high <- rep(c(Spist_class_sum[1]), 14)
Spist_low <- rep(c(Spist_class_sum[2]), 14)
Spist_all <- cbind(Spist_table_matrix, Spist_high, Spist_low)
Spist_fisher <- as.vector(apply(Spist_all,1, function(x) fisher.test(matrix(x,nr=2))$p.value))


# Create vector of asterisks according to p-value

Ahya_asterisks = rep(0, length(Ahya_fisher))
for (i in 1:length(Ahya_fisher)) {
  if (Ahya_fisher[i] <= 0.001) {
    Ahya_asterisks[i] = "***"}
  else if (Ahya_fisher[i] >= 0.001 & Ahya_fisher[i] <= 0.01) {
    Ahya_asterisks[i] = "**"}
  else if (Ahya_fisher[i] >= 0.01 & Ahya_fisher[i] <= 0.05) {
    Ahya_asterisks[i] = "*"}
  else if (Ahya_fisher[i] >= 0.05) {  
    Ahya_asterisks[i] = " "}
}

Amil_asterisks = rep(0, length(Amil_fisher))
for (i in 1:length(Amil_fisher)) {
  if (Amil_fisher[i] <= 0.001) {
    Amil_asterisks[i] = "***"}
  else if (Amil_fisher[i] >= 0.001 & Amil_fisher[i] <= 0.01) {
    Amil_asterisks[i] = "**"}
  else if (Amil_fisher[i] >= 0.01 & Amil_fisher[i] <= 0.05) {
    Amil_asterisks[i] = "*"}
  else if (Amil_fisher[i] >= 0.05) {  
    Amil_asterisks[i] = " "}
}

Apalm_asterisks = rep(0, length(Apalm_fisher))
for (i in 1:length(Apalm_fisher)) {
  if (Apalm_fisher[i] <= 0.001) {
    Apalm_asterisks[i] = "***"}
  else if (Apalm_fisher[i] >= 0.001 & Apalm_fisher[i] <= 0.01) {
    Apalm_asterisks[i] = "**"}
  else if (Apalm_fisher[i] >= 0.01 & Apalm_fisher[i] <= 0.05) {
    Apalm_asterisks[i] = "*"}
  else if (Apalm_fisher[i] >= 0.05) {  
    Apalm_asterisks[i] = " "}
}

Pdam_asterisks = rep(0, length(Pdam_fisher))
for (i in 1:length(Pdam_fisher)) {
  if (Pdam_fisher[i] <= 0.001) {
    Pdam_asterisks[i] = "***"}
  else if (Pdam_fisher[i] >= 0.001 & Pdam_fisher[i] <= 0.01) {
    Pdam_asterisks[i] = "**"}
  else if (Pdam_fisher[i] >= 0.01 & Pdam_fisher[i] <= 0.05) {
    Pdam_asterisks[i] = "*"}
  else if (Pdam_fisher[i] >= 0.05) {  
    Pdam_asterisks[i] = " "}
}

Past_asterisks = rep(0, length(Past_fisher))
for (i in 1:length(Past_fisher)) {
  if (Past_fisher[i] <= 0.001) {
    Past_asterisks[i] = "***"}
  else if (Past_fisher[i] >= 0.001 & Past_fisher[i] <= 0.01) {
    Past_asterisks[i] = "**"}
  else if (Past_fisher[i] >= 0.01 & Past_fisher[i] <= 0.05) {
    Past_asterisks[i] = "*"}
  else if (Past_fisher[i] >= 0.05) {  
    Past_asterisks[i] = " "}
}

Spist_asterisks = rep(0, length(Spist_fisher))
for (i in 1:length(Spist_fisher)) {
  if (Spist_fisher[i] <= 0.001) {
    Spist_asterisks[i] = "***"}
  else if (Spist_fisher[i] >= 0.001 & Spist_fisher[i] <= 0.01) {
    Spist_asterisks[i] = "**"}
  else if (Spist_fisher[i] >= 0.01 & Spist_fisher[i] <= 0.05) {
    Spist_asterisks[i] = "*"}
  else if (Spist_fisher[i] >= 0.05) {  
    Spist_asterisks[i] = " "}
}

#Create single datafile with means, SEs, and asterisks for plotting

GO <- as.data.frame(cbind(PdamMean, SpistMean, ApalmMean, PastMean, AmilMean, AhyaMean, PdamSE, SpistSE, ApalmSE, PastSE, AmilSE, AhyaSE))
GO2 <- GO[order(row.names(GO),decreasing=FALSE),]
GO3 <- cbind(GO2, Pdam_asterisks[2:15], Spist_asterisks[1:14], Apalm_asterisks[2:15], Past_asterisks[2:15], Amil_asterisks[2:15], Ahya_asterisks[2:15])

#Plot multiple barplots

par(mar = c(5,13,2,1))
par(mfrow = c(2, 3)) # 2 x 3 plots

Ahya_plot<-GO3[order(GO[,6],decreasing=FALSE),]
Ahya_means<-Ahya_plot[,6]
Ahya_se<-Ahya_plot[,12]
Ahya_mp<-barplot(Ahya_plot[,6], plot=FALSE, beside=TRUE, horiz=TRUE, las=2, col=c(40,42,44,46,48,50))
barplot(Ahya_plot[,6],xlim=c(0.5,0.84), axes=FALSE, beside=TRUE, xpd=F, horiz=TRUE, las=2, names.arg=row.names(Ahya_plot), col=0, cex.axis = 0.7, cex.lab = 0.7, cex = 0.7, main = "Acropora hyacinthus", font.main = 3, cex.main = 0.8)
segments(Ahya_means - Ahya_se, Ahya_mp, Ahya_means + Ahya_se, Ahya_mp)
axis(side =1, cex.axis = 0.7)
end_segment <- Ahya_means + (Ahya_se + 0.02)
text(x = end_segment, y = Ahya_mp, label = Ahya_plot[,18], cex = 0.8)

Amil_plot<-GO3[order(GO[,5],decreasing=FALSE),]
Amil_means<-Amil_plot[,5]
Amil_se<-Amil_plot[,11]
Amil_mp<-barplot(Amil_plot[,5], plot=FALSE, beside=TRUE, horiz=TRUE, las=2, col=c(40,42,44,46,48,50))
barplot(Amil_plot[,5],xlim=c(0.5,0.84), axes=FALSE, beside=TRUE, xpd=F, horiz=TRUE, las=2, names.arg=row.names(Amil_plot), col=0, cex.axis = 0.7, cex.lab = 0.7, cex = 0.7, main = "Acropora millepora", font.main = 3, cex.main = 0.8)
segments(Amil_means - Amil_se, Amil_mp, Amil_means + Amil_se, Amil_mp)
axis(side =1, cex.axis = 0.7)
end_segment <- Amil_means + (Amil_se + 0.02)
text(x = end_segment, y = Amil_mp, label = Amil_plot[,17], cex = 0.8)

Apalm_plot<-GO3[order(GO[,3],decreasing=FALSE),]
Apalm_means<-Apalm_plot[,3]
Apalm_se<-Apalm_plot[,9]
Apalm_mp<-barplot(Apalm_plot[,3], plot=FALSE, beside=TRUE, horiz=TRUE, las=2, col=c(40,42,44,46,48,50))
barplot(Apalm_plot[,3],xlim=c(0.5,0.84), axes=FALSE, beside=TRUE, xpd=F, horiz=TRUE, las=2, names.arg=row.names(Apalm_plot), col=0, cex.axis = 0.7, cex.lab = 0.7, cex = 0.7, main = "Acropora palmata", font.main = 3, cex.main = 0.8)
segments(Apalm_means - Apalm_se, Apalm_mp, Apalm_means + Apalm_se, Apalm_mp)
axis(side =1, cex.axis = 0.7)
end_segment <- Apalm_means + (Apalm_se + 0.02)
text(x = end_segment, y = Apalm_mp, label = Apalm_plot[,15], cex = 0.8)

Pdam_plot<-GO3[order(GO[,1],decreasing=FALSE),]
Pdam_means<-Pdam_plot[,1]
Pdam_se<-Pdam_plot[,7]
Pdam_mp<-barplot(Pdam_plot[,1], plot=FALSE, beside=TRUE, horiz=TRUE, las=2, col=c(40,42,44,46,48,50))
barplot(Pdam_plot[,1],xlim=c(0.5,0.84), axes=FALSE, beside=TRUE, xpd=F, horiz=TRUE, las=2, names.arg=row.names(Pdam_plot), col=0, cex.axis = 0.7, cex.lab = 0.7, cex = 0.7, main = "Pocillopora damicornis", font.main = 3, cex.main = 0.8)
segments(Pdam_means - Pdam_se, Pdam_mp, Pdam_means + Pdam_se, Pdam_mp)
axis(side =1, cex.axis = 0.7)
end_segment <- Pdam_means + (Pdam_se + 0.02)
text(x = end_segment, y = Pdam_mp, label = Pdam_plot[,13], cex = 0.8)

Past_plot<-GO3[order(GO[,4],decreasing=FALSE),]
Past_means<-Past_plot[,4]
Past_se<-Past_plot[,10]
Past_mp<-barplot(Past_plot[,4], plot=FALSE, beside=TRUE, horiz=TRUE, las=2, col=c(40,42,44,46,48,50))
barplot(Past_plot[,4],xlim=c(0.5,0.84), axes=FALSE, beside=TRUE, xpd=F, horiz=TRUE, las=2, names.arg=row.names(Past_plot), col=0, xlab = "CpG O/E", cex.axis = 0.7, cex.lab = 0.7, cex = 0.7, main = "Porites astreoides", font.main = 3, cex.main = 0.8)
segments(Past_means - Past_se, Past_mp, Past_means + Past_se, Past_mp)
axis(side =1, cex.axis = 0.7)
end_segment <- Past_means + (Past_se + 0.02)
text(x = end_segment, y = Past_mp, label = Past_plot[,16], cex = 0.8)

Spist_plot<-GO3[order(GO[,2],decreasing=FALSE),]
Spist_means<-Spist_plot[,2]
Spist_se<-Spist_plot[,8]
Spist_mp<-barplot(Spist_plot[,2], plot=FALSE, beside=TRUE, horiz=TRUE, las=2, col=c(40,42,44,46,48,50))
barplot(Spist_plot[,2],xlim=c(0.5,0.84), axes=FALSE, beside=TRUE, xpd=F, horiz=TRUE, las=2, names.arg=row.names(Spist_plot), col=0, cex.axis = 0.7, cex.lab = 0.7, cex = 0.7, main = "Stylophora pistillata", font.main = 3, cex.main = 0.8)
segments(Spist_means - Spist_se, Spist_mp, Spist_means + Spist_se, Spist_mp)
axis(side =1, cex.axis = 0.7)
end_segment <- Spist_means + (Spist_se + 0.02)
text(x = end_segment, y = Spist_mp, label = Spist_plot[,14], cex = 0.8)


######Stats


#Levene's test of homogeneity of variances
library(car)
Ahya_levene<-leveneTest(V2 ~ V3, data = Ahya2) #Ahya data is homoscedastic
Ahya_levene

#ANOVA
Ahya_aov<-aov(V2 ~ V3, data = Ahya2) #Highly sig. different (p = 2e-16)
Ahya_aov
summary(Ahya_aov)

#Tukey HSD post-hoc test
Ahya_tukey<-TukeyHSD(Ahya_aov)
Ahya_tukey
Ahya_pvalues<-Ahya_tukey$V3[,4]

#Write output to excel spreadsheet
library(xlsx)
write.xlsx(Ahya_pvalues, "/Users/jd/Documents/Projects/Coral-CpG-ratio-MS/analyses/Ahya_pvalues.xlsx") 

#Repeat for next datasets

Amil_levene<-leveneTest(V2 ~ V3, data = Amil2) 
Amil_levene
Amil_aov<-aov(V2 ~ V3, data = Amil2) 
Amil_aov
summary(Amil_aov)
Amil_tukey<-TukeyHSD(Amil_aov)
Amil_tukey
Amil_pvalues<-Amil_tukey$V3[,4]
write.xlsx(Amil_pvalues, "/Users/jd/Documents/Projects/Coral-CpG-ratio-MS/analyses/Amil_pvalues.xlsx") 


Apalm_levene<-leveneTest(V2 ~ V3, data = Apalm2) 
Apalm_levene
Apalm_aov<-aov(V2 ~ V3, data = Apalm2) 
Apalm_aov
summary(Apalm_aov)
Apalm_tukey<-TukeyHSD(Apalm_aov)
Apalm_tukey
Apalm_pvalues<-Apalm_tukey$V3[,4]
write.xlsx(Apalm_pvalues, "/Users/jd/Documents/Projects/Coral-CpG-ratio-MS/analyses/Apalm_pvalues.xlsx") 


Pdam_levene<-leveneTest(V2 ~ V3, data = Pdam2) 
Pdam_levene
Pdam_aov<-aov(V2 ~ V3, data = Pdam2) 
Pdam_aov
summary(Pdam_aov)
Pdam_tukey<-TukeyHSD(Pdam_aov)
Pdam_tukey
Pdam_pvalues<-Pdam_tukey$V3[,4]
write.xlsx(Pdam_pvalues, "/Users/jd/Documents/Projects/Coral-CpG-ratio-MS/analyses/Pdam_pvalues.xlsx") 


Spist_levene<-leveneTest(V2 ~ V3, data = Spist2) 
Spist_levene
Spist_aov<-aov(V2 ~ V3, data = Spist2) 
Spist_aov
summary(Spist_aov)
Spist_tukey<-TukeyHSD(Spist_aov)
Spist_tukey
Spist_pvalues<-Spist_tukey$V3[,4]
write.xlsx(Spist_pvalues, "/Users/jd/Documents/Projects/Coral-CpG-ratio-MS/analyses/Spist_pvalues.xlsx") 


Past_levene<-leveneTest(Past3 ~ V2, data = Past5) 
Past_levene
Past_aov<-aov(Past3 ~ V2, data = Past5) 
Past_aov
summary(Past_aov)
Past_tukey<-TukeyHSD(Past_aov)
Past_tukey
Past_pvalues<-Past_tukey$V3[,4]
write.xlsx(Past_pvalues, "/Users/jd/Documents/Projects/Coral-CpG-ratio-MS/analyses/Past_pvalues.xlsx") 

###All datasets except Ahya violate assumption of equal variances, but not severely. Examination of variances
#by group indicates that the largest variance is never more than about 2-3 times larger than the smallest.

Ahya_var<-tapply(Ahya2$V2, Ahya2$V3, var)
Amil_var<-tapply(Amil2$V2, Amil2$V3, var)
Apalm_var<-tapply(Apalm2$V2, Apalm2$V3, var)
Pdam_var<-tapply(Pdam2$V2, Pdam2$V3, var)
Spist_var<-tapply(Spist2$V2, Spist2$V3, var)
Past_var<-tapply(Past5$Past3, Past5$V2, var)



##Plot single barplot ordered by grand mean

GO3<-GO2[order(GO2[,13],decreasing=FALSE),]
GO4<-t(GO3)
par(mar=c(5,18,2,2))
midpoints<-barplot(GO4[1:6,],xlim=c(0.5,0.8),plot=FALSE, beside=TRUE, xpd=F, horiz=TRUE, las=2, col=c(1,2,3,4,5,6))
means<-GO4[1:6,]
SE<-GO4[7:12]
library(RColorBrewer)
colors<-brewer.pal(6, "Set3")
barplot(GO4[1:6,],xlim=c(0.5,0.9), axes=FALSE, beside=TRUE, xpd=F, horiz=TRUE, las=2, col=colors)
segments(means - SE, midpoints, means + SE, midpoints)
axis(side =1)
title(xlab="CpG O/E")
names<-c("P. damicornis", "S. pistillata", "A. palmata", "P. astreoides", "A. millepora", "A. hyacinthus")
revnames<-rev(names)
revcolors<-rev(colors)
legend(x = .75, y = 28, legend = revnames, fill = revcolors, border = 1, bty = "n", text.font = 3, cex =0.75)

##Plot heatmap order by grand mean

GO5 <- GO[,1:6]
means <- rowMeans(GO5)
GO6 <- cbind(GO5, means)
GO7 <- GO6[order(GO6[,7],decreasing=TRUE),]
GOmatrix <- data.matrix(GO7)

library(RColorBrewer)

cols <- colorRampPalette(brewer.pal(8,"YlGnBu"))(length(GOmatrix))

GOmatrix2 <- GOmatrix[, c("AhyaMean", "AmilMean", "ApalmMean", "PastMean", "PdamMean", "SpistMean", "means")]

setwd("~/Documents/Projects/Coral-CpG-ratio-MS/analyses")

asterisks <- read.delim("asterisks.txt", header=T) #asterisks from results of Fisher tests

par(oma=c(3,5,2,2))

GO_heatmap <- heatmap.2(GOmatrix2[,1:6], scale="none", Rowv=FALSE, Colv=FALSE,
                        col = cols, dendrogram = "none", trace = "none", density.info = "none",
                        labCol = c("A. hyacinthus", "A. millepora", "A. palmata", "P. astreoides", 
                                   "P. damicornis", "S. pistillata"), srtCol = 45,
                        lmat=rbind(c(4, 2), c(1, 3)), lhei=c(2, 8), lwid=c(4, 1), key.xlab = "CpG O/E",
                        key.title = NA, cexCol = 0.95, cexRow = 0.95, cellnote = asterisks, notecol="black")


