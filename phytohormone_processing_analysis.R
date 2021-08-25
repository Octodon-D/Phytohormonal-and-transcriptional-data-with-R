# ------------------------------------------------------
# ------------------------------------------------------
# Phytohormone - SKRIPT
# ------------------------------------------------------
# ------------------------------------------------------

# ------------------------------------------------------
# SETUP
# ------------------------------------------------------
# load reqired libraries
library(psych)
library(car)
library(lme4)
library(multcomp)
library(lsmeans)
library(tibble)

# function for standard error
std.err <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x[!is.na(x)]))

# ------------------------------------------------------
# DATA PROCESSING
# ------------------------------------------------------
# function to process raw phytohormone data to phytohormone values per ng fresh weight
# within this function the raw values get cleared by the known amounts of the internal standards and 
# normalized by the known fresh weight of the sample
phy.pro<-function(Datei, ISTD_SA=20, ISTD_ABA=20, ISTD_JA=60.4, ISTD_JA.Ile=20) 
{
  DF<-read.table(Datei, header=TRUE)
  DF <- subset(DF,DF$type=="sample")
  D6_JA.Ile_ratio<-(DF$D6_JA.Ile/(DF$D6_JA.Ile.impurity + DF$D6_JA.Ile))
  mean_D6_JA.Ile<- mean(D6_JA.Ile_ratio,na.rm=TRUE)
  DF$SA_ng.fw <- ((((DF$SA/ DF$D4_SA)*ISTD_SA)*1000)/DF$fw_mg)
  DF$ABA_ng.fw <- ((((DF$ABA/ DF$D6_ABA)*ISTD_ABA)*1000) /DF$fw_mg)
  DF$JA_ng.fw <- ((((DF$JA/ DF$D6_JA)*ISTD_JA)*1000) /DF$fw_mg)
  DF$JA_Ile_ng.fw <- (((((DF$JA.Ile/ DF$D6_JA.Ile)*ISTD_JA.Ile)*1000) /DF$fw_mg)*mean_D6_JA.Ile)
  DF$JA_Ile_per_JA_ratio <- DF$JA_Ile_ng.fw/DF$JA_ng.fw
  return(DF)
}

# process the raw data
dataframe <- phy.pro("raw_input_data.txt")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# NOTE: Until now the skipt proceeded automatized without adaptations according to the data. From now 
# on the script will include certain adaptations according to the given data structure in this experiment 
# (e.g. different leaves, different treatments). 
# For this example skript we will only consider the phytohormones in leaf 0 to prevent confusion.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# change order of the feature leaf
dataframe$leaf <- as.factor(dataframe$leaf)
dataframe$leaf <- factor(dataframe$leaf,levels=c("m1","0","1","2","3","4","5","6","7","8"))

# combine treatment with leaf position
dataframe$treat.leaf <- (gsub(" ", "",paste(dataframe$treatment, dataframe$leaf)))

# change order of treat.leaf
dataframe$treat.leaf <- as.factor(dataframe$treat.leaf)
dataframe$treat.leaf <- factor(dataframe$treat.leaf, 
                               levels=c("Cm1", "C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8",
                                        "P0", "P5", "P8", 
                                        "T0", "T5", "T8", "PT0", "PT5", "PT8"))

# sort dataframe according to treat.leaf
dataframe <- dataframe[order(dataframe$treat.leaf),]

# save dataframe
write.table(dataframe, file = "dataframe.txt",sep=",")

# ------------------------------------------------------
# STATISTICS
# ------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# SA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# model
model_0_SA <- lmer((SA_ng.fw) ~ priming * triggering + (1|replicate), REML=FALSE, data=DF.0)

# get a diagnostic plot with estimated and actual random intercepts as well as residuals
qint <- ranef(model_0_SA)$replicate[["(Intercept)"]]
qres <- residuals(model_0_SA)
layout(matrix(c(1:2),1,2,byrow=T))
qqnorm(qint, ylab="Estimated random intercepts", main="Random intercepts")
qqline(qint)
qqnorm(qres, ylab="Estimated residuals", main="Residuals")
qqline(qres)

# get model results
summary(model_0_SA)
Anova(model_0_SA)
plot(model_0_SA)
cftest(model_0_SA)
lsmeans(model_0_SA, pairwise ~ priming * triggering, adjust="none")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# ABA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# model
model_0_ABA <- lmer((ABA_ng.fw) ~ priming * triggering + (1|replicate), REML=FALSE, data=DF.0)

# get a diagnostic plot with estimated and actual random intercepts as well as residuals
qint <- ranef(model_0_ABA)$replicate[["(Intercept)"]]
qres <- residuals(model_0_ABA)
layout(matrix(c(1:2),1,2,byrow=T))
qqnorm(qint, ylab="Estimated random intercepts", main="Random intercepts")
qqline(qint)
qqnorm(qres, ylab="Estimated residuals", main="Residuals")
qqline(qres)

# get model results
summary(model_0_ABA)
Anova(model_0_ABA)
plot(model_0_ABA)
cftest(model_0_ABA)
lsmeans(model_0_ABA, pairwise ~ priming * triggering, adjust="none")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# JA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# model
model_0_JA<- lmer((JA_ng.fw) ~ priming * triggering + (1|replicate), REML=FALSE, data=DF.0)

# get a diagnostic plot with estimated and actual random intercepts as well as residuals
qint <- ranef(model_0_JA)$replicate[["(Intercept)"]]
qres <- residuals(model_0_JA)
layout(matrix(c(1:2),1,2,byrow=T))
qqnorm(qint, ylab="Estimated random intercepts", main="Random intercepts")
qqline(qint)
qqnorm(qres, ylab="Estimated residuals", main="Residuals")
qqline(qres)

# get model results
summary(model_0_JA)
Anova(model_0_JA)
plot(model_0_JA)
cftest(model_0_JA)
lsmeans(model_0_JA, pairwise ~ priming * triggering, adjust="none")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# JA-Ile
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# model
model_0_JA_Ile<- lmer((JA_Ile_ng.fw) ~ priming * triggering + (1|replicate), REML=FALSE, data=DF.0)

# get a diagnostic plot with estimated and actual random intercepts as well as residuals
qint <- ranef(model_0_JA_Ile)$replicate[["(Intercept)"]]
qres <- residuals(model_0_JA_Ile)
layout(matrix(c(1:2),1,2,byrow=T))
qqnorm(qint, ylab="Estimated random intercepts", main="Random intercepts")
qqline(qint)
qqnorm(qres, ylab="Estimated residuals", main="Residuals")
qqline(qres)

# get model results
summary(model_0_JA_Ile)
Anova(model_0_JA_Ile)
plot(model_0_JA_Ile)
cftest(model_0_JA_Ile)
lsmeans(model_0_JA_Ile, pairwise ~ priming * triggering, adjust="none")

# ------------------------------------------------------
# VISUALIZATION: BOXPLOTS
# ------------------------------------------------------
# preperation: subset dataframe
DF.0 <- subset(dataframe, dataframe$leaf=="0")
DF.0 <- droplevels(DF.0)
DF.0$treatment <- factor(DF.0$treatment, levels=c("C", "P", "T", "PT"))

# boxplot combined with a scatterplott (each replicate with different icon) with the four different pythohormones
# settings
layout(matrix(c(1:6),3,2,byrow=T))
par(las =2)

# leaf 0
# SA
boxplot(SA_ng.fw ~ treatment , vertical = TRUE,
        main="SA (ng/g fw)",ylab="ng/g fw", data=DF.0)
for(z in (1:length(unique(DF.0$replicate)))) {
  x<-subset(DF.0, DF.0$replicate==unique(DF.0$replicate)[z])
  stripchart(SA_ng.fw ~ treatment, col="orangered", vertical=TRUE, add=T, pch=14+z, data=x, cex=2)
}
# ABA
boxplot(ABA_ng.fw ~ treatment , vertical=TRUE,
        main="ABA (ng/g fw)",ylab="ng/g fw", data=DF.0)
for(z in (1:length(unique(DF.0$replicate)))) {
  x<-subset(DF.0, DF.0$replicate==unique(DF.0$replicate)[z])
  stripchart(ABA_ng.fw ~ treatment, col="orangered", vertical=TRUE, add=T, pch=14+z, data=x, cex=2)
}
# JA
boxplot(JA_ng.fw ~ treatment, vertical=TRUE,
        main="JA (ng/g fw)",ylab="ng/g fw", data=DF.0)
for(z in (1:length(unique(DF.0$replicate)))) {
  x<-subset(DF.0, DF.0$replicate==unique(DF.0$replicate)[z])
  stripchart(JA_ng.fw ~ treatment, col="orangered", vertical=TRUE, add=T, pch=14+z, data=x, cex=2)
}
# JA_Ile
boxplot(JA_Ile_ng.fw ~ treatment, vertical=TRUE,
        main="JA.Ile (ng/g fw)",ylab="ng/g fw", data=DF.0)
for(z in (1:length(unique(DF.0$replicate)))) {
  x<-subset(DF.0, DF.0$replicate==unique(DF.0$replicate)[z])
  stripchart(JA_Ile_ng.fw ~ treatment, col="orangered", vertical=TRUE, add=T, pch=14+z, data=x, cex=2)
}
# JA_Ile_per_JA_ratio
boxplot(JA_Ile_per_JA_ratio ~ treatment, vertical=TRUE,
        main="JA.Ile_ng/JA_ng ratio",ylab="ratio JA-Ile/JA", data=DF.0, ylim=c(0,0.1))
for(z in (1:length(unique(DF.0$replicate)))) {
  x<-subset(DF.0,DF.0$replicate==unique(DF.0$replicate)[z])
  stripchart(JA_Ile_per_JA_ratio ~ treatment, col="orangered", vertical=TRUE, add=T, pch=14+z, data=x, cex=2)
}
# fw_mg
boxplot(fw_mg ~ treatment, vertical=TRUE,
        main="fresh weight used",ylab="fw_mg", data=DF.0)
for(z in (1:length(unique(DF.0$replicate)))) {
  x<-subset(DF.0,DF.0$replicate==unique(DF.0$replicate)[z])
  stripchart(fw_mg ~ treatment, col="orangered", vertical=TRUE, add=T, pch=14+z, data=x, cex=2)
}

# ------------------------------------------------------------------------------
# VISUALIZATION: BARPLOTS
# ------------------------------------------------------------------------------
# preperation to generate barplots with mean values
# create an empty vector
mean.0 <- matrix(data=NA, nrow=4, ncol=4, byrow=FALSE, dimnames=NULL)
SE.0 <- matrix(data=NA, nrow=4, ncol=4, byrow=FALSE, dimnames=NULL)

# loop to fill in empty vector with mean values of certain leaf and treatment
for(z in (1:length(unique(DF.0$treatment)))) {
  x<-subset(DF.0,DF.0$treatment==unique(DF.0$treatment)[z])
  mean.0[1,z] <- mean(x$SA_ng.fw, na.rm=T)
  SE.0[1,z] <- std.err(x$SA_ng.fw)
  mean.0[2,z] <- mean(x$ABA_ng.fw, na.rm=T)
  SE.0[2,z] <- std.err(x$ABA_ng.fw)
  mean.0[3,z] <- mean(x$JA_ng.fw, na.rm=T)
  SE.0[3,z] <- std.err(x$JA_ng.fw)
  mean.0[4,z] <- mean(x$JA_Ile_ng.fw, na.rm=T)
  SE.0[4,z] <- std.err(x$JA_Ile_ng.fw)
}

# barplot: generate a barplot and save it directly in a svg file
# settings for the barplots
angle<-c(0.6,1.7,2.8,3.9)
mycol<-c("gray100","#ffaaaaff","#333333ff","#ff0000ff")
par(mar =c(0.25,2.5,0.25,0.25), cex=1)

# SA leaf 0
svg(file="phytohormone_figureSA_L0.svg", width=1.5748, height=1.1811)
barplot(mean.0[1,], 
        xlab="", ylab="", main="",
        ylim=c(0,150),
        col=mycol,
        xaxt='n',
        las=2,
        space=0.1,
        lwd=0.5,
        tck=-0.05)
MP<-mean.0[1,]+SE.0[1,]
MM<-mean.0[1,]-SE.0[1,]
segments(angle, MP, angle,MM)
segments(angle-0.1, MP, angle+0.1, MP)
segments(angle-0.1, MM, angle+0.1, MM)
dev.off()

# ABA leaf 0
svg(file="phytohormone_figureABA_L0.svg", width=1.9685, height=1.37795)
barplot(mean.0[2,], 
        xlab="", ylab="", main="",
        ylim=c(0,500),
        col=mycol,
        xaxt='n',
        las=2,
        space=0.1,
        lwd=0.5,
        tck=-0.05)
MP<-mean.0[2,]+SE.0[2,]
MM<-mean.0[2,]-SE.0[2,]
segments(angle, MP, angle,MM)
segments(angle-0.1, MP, angle+0.1, MP)
segments(angle-0.1, MM, angle+0.1, MM)
dev.off()

# JA leaf 0
svg(file="phytohormone_figureJA_L0.svg", width=1.9685, height=1.37795)
barplot(mean.0[3,], 
        xlab="", ylab="", main="",
        ylim=c(0,5),
        col=mycol,
        xaxt='n',
        las=2,
        space=0.1,
        lwd=0.5,
        tck=-0.05)
MP<-mean.0[3,]+SE.0[3,]
MM<-mean.0[3,]-SE.0[3,]
segments(angle, MP, angle,MM)
segments(angle-0.1, MP, angle+0.1, MP)
segments(angle-0.1, MM, angle+0.1, MM)
dev.off()

# JA-Ile leaf 0
svg(file="phytohormone_figureJA-Ile_01L0.svg", width=1.9685, height=1.37795)
barplot(mean.0[4,], 
        xlab="", ylab="", main="",
        ylim=c(0,5),
        col=mycol,
        xaxt='n',
        las=2,
        space=0.1,
        lwd=0.5,
        tck=-0.05)
MP<-mean.0[4,]+SE.0[4,]
MM<-mean.0[4,]-SE.0[4,]
segments(angle, MP, angle,MM)
segments(angle-0.1, MP, angle+0.1, MP)
segments(angle-0.1, MM, angle+0.1, MM)
dev.off()