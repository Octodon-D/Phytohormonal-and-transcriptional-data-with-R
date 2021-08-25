# ------------------------------------------------------
# ------------------------------------------------------
# qPCR - MASTER - SKRIPT
# ------------------------------------------------------
# ------------------------------------------------------

# ------------------------------------------------------
# SETUP
# ------------------------------------------------------
# load required packages
library(plyr)
library(car)
library(lme4)
library(psych)
library(multcomp)
library(emmeans)
library(ComplexHeatmap)
library(dendextend)
library(circlize)

# function for standard error
std.err <- function(x) sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))

# set number of plates
number_of_plates <- 4

# set CRITICAL VALUE for outlier correction:
# define the critical value for the outlier correction when comparing the mesurements within each replicate. It is based on
# the distance of the datapoint (measurement) from the median of the tiplicate and distance of the other (most distant) 
# datapoint (measurement) from the median of the triplicate. If the distance is more than the defined value, the datapoint 
# is exluded from further calculations. The datapoint is kept if its distance to the median of the triplicate is smaller than
# the defined value. 
# In simple words: If one extreme point is more times (than defined by the value) further away from the middle than the 
# other one, it is kicked out.
outlier_critical <- 3

# ------------------------------------------------------
# DATA PROCESSING (PART 1)
# ------------------------------------------------------
# load parameters file & delete empty rows
parameters_frame <- read.table("parameters.txt", header=T)
parameters_frame <- parameters_frame[(is.na(parameters_frame$sample)==FALSE),]

# load output-files into a result dataframe & convert threshold and ct into numeric
raw_results <- read.table("input_data_R.txt", fill=TRUE, header=TRUE)
raw_results$Plate <- c(rep(1:number_of_plates, each=96))

# add position column (combine plate and Well)
raw_results$position <- (gsub(" ", "",paste(raw_results$Plate, raw_results$name)))

# read Linreg_grouping and assignment of Samples to wells and genes in
layout <- read.table("layout.txt", header=T)

# create wells
well_letter <- rep(c("A", "B", "C", "D", "E", "F", "G", "H"), (times=number_of_plates*12))
well_number <- rep(c(1:12), each=(number_of_plates*8))

# create plate numbers for new dataframe
plate <- rep(1:number_of_plates, each=8,12)

# create a common identifier "position" for results and layout 
# FIRST: plate, SECOND: well_letter, THIRD: well_number
position <- (gsub(" ", "",paste(plate, well_letter, well_number)))

# Combine results with sample_id and gen_id in "dataframe":
# create new dataframe with position as first column
dataframe <- data.frame(position)

# add plate number
dataframe$plate <- plate

# add sample_id
dataframe$samples <- layout$samples[match(dataframe$position, layout$position)]

# add type_id 
dataframe$gene <- layout$gene[match(dataframe$position, layout$position)]

# generate (sample + gene) and add replicate
dataframe$replicate <- as.factor(gsub(" ", ".",paste(dataframe$samples, dataframe$gene)))

# match results into dataframe
dataframe$ct <- raw_results$Cq[match(dataframe$position, raw_results$position)]
dataframe$threshold <- raw_results$threshold[match(dataframe$position, raw_results$position)]
dataframe$efficiency <- raw_results$mean_PCR_eff[match(dataframe$position, raw_results$position)]
dataframe$individual.efficiency <- raw_results$indiv_PCR_eff[match(dataframe$position, raw_results$position)]
dataframe$amplicon <- raw_results$Amplicon[match(dataframe$position, raw_results$position)]

# remove not measured wells
dataframe <- dataframe [(is.na(dataframe$samples)==FALSE),]
dataframe$replicate <- droplevels(dataframe$replicate)

# add treatment from parameters
dataframe$treatment <- parameters_frame$treatment[match(dataframe$samples, parameters_frame$sample)]
dataframe$factor1 <- parameters_frame$factor_1[match(as.factor(dataframe$samples), parameters_frame$sample)]
dataframe$factor2 <- parameters_frame$factor2[match(as.factor(dataframe$samples), parameters_frame$sample)]
dataframe$factor3 <- parameters_frame$factor_3[match(as.factor(dataframe$samples), parameters_frame$sample)]

# sort dataframe by replicate to match the output of describeBy for outlier correction
dataframe <- dataframe[order(dataframe$replicate),]

# ------------------------------------------------------
# OUTLIER CORRECTION
# ------------------------------------------------------
# get descriptive statistics for sortframe_measured_noblank used for later outlier correction
a <- describeBy(dataframe$ct, dataframe$replicate, na.rm=TRUE)

ct_mean <- vector("integer", length=length(unique(dataframe$replicate)))
sd.ct_mean <- vector("integer", length=length(unique(dataframe$replicate)))
result.number <- vector("integer", length=length(unique(dataframe$replicate)))
for(i in 1:length(unique(dataframe$replicate)))  {
  ct_mean[i] <- a[[0+i]][["mean"]][[1]]
  sd.ct_mean[i] <- a[[0+i]][["sd"]][[1]]
  result.number[i] <- a[[0+i]][["n"]][[1]]
}

ct_mean <- as.numeric(ct_mean)
sd.ct <- as.numeric(sd.ct_mean)

# get descriptive statistics for sortframe_measured_noblank used for later outlier correction
a <- describeBy(dataframe$ct, dataframe$replicate)

# get number frequency of technical replicates (should be 3 for every triplicate)
triplicates_left <- count(dataframe, vars="replicate")
triplicates_freq <- triplicates_left$"freq"

# get number of triplicates to calculate for the loops
triplicates_num <- nrow(triplicates_left)

# get number of actually measured ct-values to add "outlier" column later
# pick ranges for each triplicate
# pick medians for each triplicate
triplicate_count <- vector("integer", length=triplicates_num)
range_single <- vector("integer", length=triplicates_num)
median_single <- vector("integer", length=triplicates_num)  
for(i in (1:triplicates_num))  {
  triplicate_count[i] <- a[[0+i]][["n"]][[1]]
  range_single[i] <- a[[0+i]][["range"]][[1]]
  median_single[i] <- a[[0+i]][["median"]][[1]]
}

# repeat the range_single and median_single according to the number of reads left in each triplicate
range <- rep(range_single, times=triplicates_freq)
median <- rep(median_single, times=triplicates_freq)

# absolute values of ranges for outlier removal (need to get rid of possible minus-signs)
range_abs <- abs(range)

# calculate absolute values of distances between measured value and the median value within each triplicate  
# (need to get rid of possible minus-signs).
median_dist_abs <- abs(median - dataframe$ct)

# perform outlier correction for the slope values based on:
# comparison between 
# 1.distance of the datapoint from the median and 
# 2.distance of the other (distant) slope from the median
# default:value is kept if its distance to the median is smaller than the distance of the other slope to the median 
# times the number specified in (slope_critical)
ct_corrected <- ifelse(median_dist_abs > ((range_abs - median_dist_abs)*outlier_critical),"NA", dataframe$ct)

# the warning message in the R-output is necessarry... everything all right!
dataframe$ct_corrected <- as.numeric(ct_corrected)

# get mean values of the final measurements within each triplicate
b <- describeBy(dataframe$ct_corrected, dataframe$replicate)

# get mean values of the standard deviation of the ct_mean within each triplicate
# get number of ct_means that were left after the outlier correction in each triplicate 
# (required for outlier report in the result table)

ct_mean <- vector("integer", length=triplicates_num)
sd.ct_mean <- vector("integer", length=triplicates_num)
result.number <- vector("integer", length=triplicates_num)
for(i in 1:triplicates_num)  {
  ct_mean[i] <- b[[0+i]][["mean"]][[1]]
  sd.ct_mean[i] <- b[[0+i]][["sd"]][[1]]
  result.number[i] <- b[[0+i]][["n"]][[1]]
}

# report in the data file in which samples outlieres were removed
outlier <- ifelse((triplicate_count > result.number), "outlier", "ok")

# ------------------------------------------------------
# DATA PROCESSING (PART 2)
# ------------------------------------------------------
# get triplicate labels & create resultframe
replicate <- unique(dataframe$replicate)
resultframe <- data.frame(replicate)

# match type id´s to triplicate labels
resultframe$plate <- dataframe$plate[match(resultframe$replicate, dataframe$replicate)]
resultframe$sample <- dataframe$samples[match(resultframe$replicate, dataframe$replicate)]
resultframe$treatment <- parameters_frame$treatment[match(resultframe$sample, parameters_frame$sample)]
resultframe$factor1 <- parameters_frame$factor_1[match(resultframe$sample, parameters_frame$sample)]
resultframe$factor2 <- parameters_frame$factor2[match(resultframe$sample, parameters_frame$sample)]
resultframe$factor3 <- parameters_frame$factor_3[match(resultframe$sample, parameters_frame$sample)]
resultframe$gene <- dataframe$gene[match(resultframe$replicate, dataframe$replicate)]
resultframe$efficiency <- dataframe$efficiency[match(resultframe$replicate, dataframe$replicate)]
resultframe$ct_mean_corrected <- ct_mean
resultframe$sd.ct_mean_correc <- sd.ct
resultframe$outlier.removed <- outlier
# sort resultframe by gene
resultframe <- resultframe[order(resultframe$gene),]

# include primer efficiency in calculation of ct-value
resultframe$eff.ct <- resultframe$efficiency^-resultframe$ct_mean

# save resultframe as csv
write.table(resultframe, file = "results.txt",sep=",")

# delete row with NTC (no information) 
resultframe_noNTC <- resultframe[resultframe$sample!="NTC",]
resultframe_noNTC <- droplevels(resultframe_noNTC)

# new dataframe summary
sample <- unique(resultframe_noNTC$sample) 
summary1 <- data.frame(sample)

# get parameters in summary
summary1$treatment <- parameters_frame$treatment[match(summary1$sample, parameters_frame$sample)]
summary1$block <- parameters_frame$block[match(summary1$sample, parameters_frame$sample)]
summary1$factor1 <- parameters_frame$factor_1[match(summary1$sample, parameters_frame$sample)]
summary1$factor2 <- parameters_frame$factor2[match(summary1$sample, parameters_frame$sample)]
summary1$factor3 <- parameters_frame$factor_3[match(summary1$sample, parameters_frame$sample)]
summary1$block <- parameters_frame$block[match(summary1$sample, parameters_frame$sample)]

# create three new and empty matrix with length sample x numbers of gens+1 
# to put in ct - means and eff.ct values from the resultframe
eff <- array(data=NA, dim=c(length(sample), length(unique(resultframe$gene))+1))
eff[,1] <- as.character(sample)
means <- array(data=NA, dim=c(length(sample), length(unique(resultframe$gene))+1))
means[,1] <- as.character(sample)
mean.sd <- array(data=NA, dim=c(length(sample), length(unique(resultframe$gene))+1))
mean.sd[,1] <- as.character(sample)
rqs <- array(data=NA,dim=c(length(sample), length(unique(resultframe$gene))+1))
rqs[,1] <- as.character(sample)

# pick mean-ct/eff.ct values for each gene out of resultframe and put into the two matrix 
# from one below the other to one beside the other
for(i in (1:length(unique(resultframe$gene))))  {
  x <- resultframe[resultframe$gene==unique(resultframe$gene)[i],]
  eff[,1+i] <- x$efficiency[match(eff[,1], x$sample)]
  means[,1+i] <- x$ct_mean[match(means[,1], x$sample)]
  mean.sd[,1+i] <- x$sd.ct_mean[match(mean.sd[,1], x$sample)]
  rqs[,1+i] <- x$eff.ct[match(rqs[,1], x$sample)]
}

# convert the two matrix into dataframes
eff <- as.data.frame(eff)
means <- as.data.frame(means)
mean.sd <- as.data.frame(mean.sd)
rqs <- as.data.frame(rqs)

# give the dataframes names
names(eff) <- c("sample", paste(unique(resultframe$gene), "_mean.eff", sep=""))
names(means) <- c("sample",paste(unique(resultframe$gene), "_mean.ct", sep=""))
names(mean.sd) <- c("sample",paste(unique(resultframe$gene), "_mean.ct.sd", sep=""))
names(rqs) <- c("sample",paste(unique(resultframe$gene), "_eff.ct", sep=""))

# combine all three dataframes 
summary2 <- merge(summary1, eff, by="sample")
summary3 <- merge(summary2, means, by="sample")
summary4 <- merge(summary3, mean.sd, by="sample")
summary <- merge(summary4, rqs, by="sample")

# convert factors into numeric
for(i in 1:(length(names(summary))-6)) {
  summary[,6+i] <- as.numeric(as.character(summary[,6+i]))
}

# ------------------------------------------------------
# PCR QUALITY EVALUATION
# ------------------------------------------------------
# make a boxplot comparing the NTC and sample ct values
par(mar = c(10,4,4,3))
boxplot(resultframe$ct_mean_corrected ~ resultframe$gene, 
        col="gray50",
        pch=20,
        las=2,
        ylim=c(0, max(resultframe$ct_mean_corrected)+1))
points(resultframe$ct_mean_corrected[resultframe$sample=="NTC"] ~ resultframe$gene[resultframe$sample=="NTC"],
       col="red",
       pch=19)
legend("topright",
       pch=c(20,19),
       col=c("black", "red"), 
       legend=c("samples", "NTC"),
       bty="n",
       inset=c(0.1, -0.1, 0.1, 0.1),
       xpd=TRUE,
       horiz=TRUE)

# check the primer efficiency for each gene
print("primer efficiency")
print(sapply(split(resultframe$efficiency, resultframe$gene), mean))

# ------------------------------------------------------
# CALCULATIONS
# ------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# NOTE: Until now the skipt proceeded automatized without adaptations. From now on the script will include 
# certain adaptations according to genes. As this is an example skript, we will have ELF1 as reference gene
# and PR1, PG4, HCT1, PRX2 and CDH as target genes/GOI´s (gene of intrest)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# put treatments in a defined order 
summary$treatment <- factor(summary$treatment, levels=c("C","P","T","PT"))

# set reference eff.ct 
summary$REF.eff.ct <- summary$ELF1_eff.ct

# calculate RQ forfor GOI´s (eff.ct relative to REF.eff.ct)
summary$PR1_RQ <- summary$PR1_eff.ct / summary$REF.eff.ct
summary$PR1_log2.RQ <- log2(summary$PR1_RQ)
summary$PG4_RQ <- summary$PG4_eff.ct / summary$REF.eff.ct
summary$PG4_log2.RQ <- log2(summary$PG4_RQ)
summary$HCT1_RQ <- summary$HCT1_eff.ct / summary$REF.eff.ct
summary$HCT1_log2.RQ <- log2(summary$HCT1_RQ)
summary$PRX2_RQ <- summary$PRX2_eff.ct / summary$REF.eff.ct
summary$PRX2_log2.RQ <- log2(summary$PRX2_RQ)
summary$CDH_RQ<-summary$CDH_eff.ct / summary$REF.eff.ct
summary$CDH_log2.RQ<-log2(summary$CDH_RQ)

# calculate log2.NRQ relative to control 
# 1. calculate mean control log2.RQ for each interesting factor
PR1_log2.RQ.mean.c <- mean(summary$PR1_log2.RQ[summary$treatment=="C"], na.rm=T)
PG4_log2.RQ.mean.c <- mean(summary$PG4_log2.RQ[summary$treatment=="C"], na.rm=T)
HCT1_log2.RQ.mean.c <- mean(summary$HCT1_log2.RQ[summary$treatment=="C"], na.rm=T)
PRX2_log2.RQ.mean.c <- mean(summary$PRX2_log2.RQ[summary$treatment=="C"], na.rm=T)
CDH_log2.RQ.mean.c <- mean(summary$CDH_log2.RQ[summary$treatment=="C"], na.rm=T)

# 2. calculate log2.NRQ relative to mean control for each interesting factor
summary$PR1_log2.NRQ <- summary$PR1_log2.RQ - PR1_log2.RQ.mean.c
summary$PG4_log2.NRQ <- summary$PG4_log2.RQ - PG4_log2.RQ.mean.c
summary$HCT1_log2.NRQ <- summary$HCT1_log2.RQ - HCT1_log2.RQ.mean.c
summary$PRX2_log2.NRQ <- summary$PRX2_log2.RQ - PRX2_log2.RQ.mean.c
summary$CDH_log2.NRQ <- summary$CDH_log2.RQ - CDH_log2.RQ.mean.c

# save summary dataframe as csv
write.csv(summary, file = "summary.txt")

# ------------------------------------------------------
# STATISTICS
# ------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# CDH
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# assumption: normal distribution
# check qq-plots with different data transformations
qqp(summary$CDH_log2.NRQ, "norm") 
qqp(log(summary$CDH_log2.NRQ), "norm") 
qqp(log10(summary$CDH_log2.NRQ), "norm") 
qqp(sqrt(summary$CDH_log2.NRQ), "norm")  
qqp(1/sqrt(summary$CDH_log2.NRQ), "norm")
qqp(1/(summary$CDH_log2.NRQ), "norm")

shapiro.test((summary$CDH_log2.NRQ))
ks.test((summary$CDH_log2.NRQ), "pnorm", mean=mean((summary$CDH_log2.NRQ), na.rm=T), sd=sd((summary$CDH_log2.NRQ), na.rm = T))

# assumption: equal variances
var.test((summary$CDH_log2.NRQ[summary$treatment=="C"]), (summary$CDH_log2.NRQ[summary$treatment=="P"])) 
var.test((summary$CDH_log2.NRQ[summary$treatment=="C"]), (summary$CDH_log2.NRQ[summary$treatment=="T"])) 
var.test((summary$CDH_log2.NRQ[summary$treatment=="C"]), (summary$CDH_log2.NRQ[summary$treatment=="PT"])) 
var.test((summary$CDH_log2.NRQ[summary$treatment=="P"]), (summary$CDH_log2.NRQ[summary$treatment=="T"])) 
var.test((summary$CDH_log2.NRQ[summary$treatment=="P"]), (summary$CDH_log2.NRQ[summary$treatment=="PT"])) 
var.test((summary$CDH_log2.NRQ[summary$treatment=="T"]), (summary$CDH_log2.NRQ[summary$treatment=="PT"])) 

# model
model_CDH <- lmer(CDH_log2.NRQ ~ factor2 * factor3 + (1|block), REML=FALSE, data=summary)

# get a diagnostic plot with estimated and actual random intercepts as well as residuals
qint <- ranef(model_CDH)$replicate[["(Intercept)"]]
qres <- residuals(model_CDH)
layout(matrix(c(1:2),1,2,byrow=T))
qqnorm(qint, ylab="Estimated random intercepts", main="Random intercepts")
qqline(qint)
qqnorm(qres, ylab="Estimated residuals", main="Residuals")
qqline(qres)

# get model results
summary(model_CDH)
Anova(model_CDH)
plot(model_CDH)
cftest(model_CDH)
lsmeans(model_CDH, pairwise ~ factor2 + factor3, adjust="none")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# HCT1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# assumption: normal distribution
# check qq-plots with different data transformations
qqp(summary$HCT1_log2.NRQ, "norm") 
qqp(log(summary$HCT1_log2.NRQ+10), "norm") 
qqp(log10(summary$HCT1_log2.NRQ+10), "norm") 
qqp(sqrt(summary$HCT1_log2.NRQ+10), "norm")  
qqp(1/sqrt(summary$HCT1_log2.NRQ+10), "norm")
qqp(1/(summary$HCT1_log2.NRQ+10), "norm")

shapiro.test(1/(summary$HCT1_log2.NRQ+10))
ks.test(1/(summary$HCT1_log2.NRQ+10), "pnorm", mean=mean(1/(summary$HCT1_log2.NRQ+10), na.rm=T), sd=sd(1/(summary$HCT1_log2.NRQ+10), na.rm=T))

# assumption: equal variances
var.test(1/(summary$HCT1_log2.NRQ[summary$treatment=="C"]+10), 1/(summary$HCT1_log2.NRQ[summary$treatment=="P"]+10)) 
var.test(1/(summary$HCT1_log2.NRQ[summary$treatment=="C"]+10), 1/(summary$HCT1_log2.NRQ[summary$treatment=="T"]+10)) 
var.test(1/(summary$HCT1_log2.NRQ[summary$treatment=="C"]+10), 1/(summary$HCT1_log2.NRQ[summary$treatment=="PT"]+10)) 
var.test(1/(summary$HCT1_log2.NRQ[summary$treatment=="P"]+10), 1/(summary$HCT1_log2.NRQ[summary$treatment=="T"]+10)) 
var.test(1/(summary$HCT1_log2.NRQ[summary$treatment=="P"]+10), 1/(summary$HCT1_log2.NRQ[summary$treatment=="PT"]+10)) 
var.test(1/(summary$HCT1_log2.NRQ[summary$treatment=="T"]+10), 1/(summary$HCT1_log2.NRQ[summary$treatment=="PT"]+10)) 

# model
model_HCT1 <- lmer(HCT1_log2.NRQ ~ factor2 * factor3 + (1|block), REML=FALSE, data=summary)

# get a diagnostic plot with estimated and actual random intercepts as well as residuals
qint <- ranef(model_HCT1)$block[["(Intercept)"]]
qres <- residuals(model_HCT1)
layout(matrix(c(1:2),1,2,byrow=T))
qqnorm(qint, ylab="Estimated random intercepts", main="Random intercepts")
qqline(qint)
qqnorm(qres, ylab="Estimated residuals", main="Residuals")
qqline(qres)

# get model results
summary(model_HCT1)
Anova(model_HCT1)
plot(model_HCT1)
cftest(model_HCT1)
lsmeans(model_HCT1, pairwise ~ factor2 + factor3, adjust="none")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# PG4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# assumption: normal distribution
# check qq-plots with different data transformations
qqp(summary$PG4_log2.NRQ, "norm") 
qqp(log(summary$PG4_log2.NRQ+10), "norm") 
qqp(log10(summary$PG4_log2.NRQ+10), "norm") 
qqp(sqrt(summary$PG4_log2.NRQ+10), "norm")  
qqp(1/sqrt(summary$PG4_log2.NRQ+10), "norm")
qqp(1/(summary$PG4_log2.NRQ+10), "norm")

shapiro.test((summary$PG4_log2.NRQ))
ks.test((summary$PG4_log2.NRQ), "pnorm", mean=mean((summary$PG4_log2.NRQ), na.rm=T), sd=sd((summary$PG4_log2.NRQ), na.rm=T))

# assumption: equal variances
var.test((summary$PG4_log2.NRQ[summary$treatment=="C"]), (summary$PG4_log2.NRQ[summary$treatment=="P"])) 
var.test((summary$PG4_log2.NRQ[summary$treatment=="C"]), (summary$PG4_log2.NRQ[summary$treatment=="T"])) 
var.test((summary$PG4_log2.NRQ[summary$treatment=="C"]), (summary$PG4_log2.NRQ[summary$treatment=="PT"])) 
var.test((summary$PG4_log2.NRQ[summary$treatment=="P"]), (summary$PG4_log2.NRQ[summary$treatment=="T"])) 
var.test((summary$PG4_log2.NRQ[summary$treatment=="P"]), (summary$PG4_log2.NRQ[summary$treatment=="PT"])) 
var.test((summary$PG4_log2.NRQ[summary$treatment=="T"]), (summary$PG4_log2.NRQ[summary$treatment=="PT"])) 

# model
model_PG4 <- lmer(PG4_log2.NRQ ~ factor2 * factor3 + (1|block), REML=FALSE, data=summary)

# get a diagnostic plot with estimated and actual random intercepts as well as residuals
qint <- ranef(model_PG4)$block[["(Intercept)"]]
qres <- residuals(model_PG4)
layout(matrix(c(1:2),1,2,byrow=T))
qqnorm(qint, ylab="Estimated random intercepts", main="Random intercepts")
qqline(qint)
qqnorm(qres, ylab="Estimated residuals", main="Residuals")
qqline(qres)

# get model results
summary(model_PG4)
Anova(model_PG4)
plot(model_PG4)
cftest(model_PG4)
lsmeans(model_PG4, pairwise ~ factor2 + factor3, adjust="none")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# PR1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# assumption: normal distribution
# check qq-plots with different data transformations
qqp(summary$PR1_log2.NRQ, "norm") 
qqp(log(summary$PR1_log2.NRQ+10), "norm") 
qqp(log10(summary$PR1_log2.NRQ+10), "norm") 
qqp(sqrt(summary$PR1_log2.NRQ+10), "norm")  
qqp(1/sqrt(summary$PR1_log2.NRQ+10), "norm")
qqp(1/(summary$PR1_log2.NRQ+10), "norm")

shapiro.test((summary$PR1_log2.NRQ))
ks.test((summary$PR1_log2.NRQ), "pnorm", mean=mean((summary$PR1_log2.NRQ), na.rm=T), sd=sd((summary$PR1_log2.NRQ), na.rm=T))

# assumption: equal variances
var.test((summary$PR1_log2.NRQ[summary$treatment=="C"]), (summary$PR1_log2.NRQ[summary$treatment=="P"])) 
var.test((summary$PR1_log2.NRQ[summary$treatment=="C"]), (summary$PR1_log2.NRQ[summary$treatment=="T"])) 
var.test((summary$PR1_log2.NRQ[summary$treatment=="C"]), (summary$PR1_log2.NRQ[summary$treatment=="PT"])) 
var.test((summary$PR1_log2.NRQ[summary$treatment=="P"]), (summary$PR1_log2.NRQ[summary$treatment=="T"])) 
var.test((summary$PR1_log2.NRQ[summary$treatment=="P"]), (summary$PR1_log2.NRQ[summary$treatment=="PT"])) 
var.test((summary$PR1_log2.NRQ[summary$treatment=="T"]), (summary$PR1_log2.NRQ[summary$treatment=="PT"])) 

# model
model_PR1 <- lmer(PR1_log2.NRQ ~ factor2 * factor3 + (1|block), REML=FALSE, data=summary)

# get a diagnostic plot with estimated and actual random intercepts as well as residuals
qint <- ranef(model_PR1)$block[["(Intercept)"]]
qres <- residuals(model_PR1)
layout(matrix(c(1:2),1,2,byrow=T))
qqnorm(qint, ylab="Estimated random intercepts", main="Random intercepts")
qqline(qint)
qqnorm(qres, ylab="Estimated residuals", main="Residuals")
qqline(qres)

# get model results
summary(model_PR1)
Anova(model_PR1)
plot(model_PR1)
cftest(model_PR1)
lsmeans(model_PR1, pairwise ~ factor2 + factor3, adjust="none")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRX2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# assumption: normal distribution
# check qq-plots with different data transformations
qqp(summary$PRX2_log2.NRQ, "norm") 
qqp(log(summary$PRX2_log2.NRQ+10), "norm") 
qqp(log10(summary$PRX2_log2.NRQ+10), "norm") 
qqp(sqrt(summary$PRX2_log2.NRQ+10), "norm")  
qqp(1/sqrt(summary$PRX2_log2.NRQ+10), "norm")
qqp(1/(summary$PRX2_log2.NRQ+10), "norm")

shapiro.test((summary$PRX2_log2.NRQ))
ks.test((summary$PRX2_log2.NRQ), "pnorm", mean=mean((summary$PRX2_log2.NRQ), na.rm=T), sd=sd((summary$PRX2_log2.NRQ), na.rm=T))

# assumption: equal variances
var.test((summary$PRX2_log2.NRQ[summary$treatment=="C"]), (summary$PRX2_log2.NRQ[summary$treatment=="P"])) 
var.test((summary$PRX2_log2.NRQ[summary$treatment=="C"]), (summary$PRX2_log2.NRQ[summary$treatment=="T"])) 
var.test((summary$PRX2_log2.NRQ[summary$treatment=="C"]), (summary$PRX2_log2.NRQ[summary$treatment=="PT"])) 
var.test((summary$PRX2_log2.NRQ[summary$treatment=="P"]), (summary$PRX2_log2.NRQ[summary$treatment=="T"])) 
var.test((summary$PRX2_log2.NRQ[summary$treatment=="P"]), (summary$PRX2_log2.NRQ[summary$treatment=="PT"])) 
var.test((summary$PRX2_log2.NRQ[summary$treatment=="T"]), (summary$PRX2_log2.NRQ[summary$treatment=="PT"])) 

# model
model_PRX2 <- lmer(PRX2_log2.NRQ ~ factor2 * factor3 + (1|block), REML=FALSE, data=summary)

# get a diagnostic plot with estimated and actual random intercepts as well as residuals
qint <- ranef(model_PRX2)$block[["(Intercept)"]]
qres <- residuals(model_PRX2)
layout(matrix(c(1:2),1,2,byrow=T))
qqnorm(qint, ylab="Estimated random intercepts", main="Random intercepts")
qqline(qint)
qqnorm(qres, ylab="Estimated residuals", main="Residuals")
qqline(qres)

# get model results
summary(model_PRX2)
Anova(model_PRX2)
plot(model_PRX2)
cftest(model_PRX2)
lsmeans(model_PRX2, pairwise ~ factor2 + factor3, adjust="none")

# ------------------------------------------------------
# VISUALIZATION: BOXPLOTS
# ------------------------------------------------------
# boxplots with mean ct values for each gene
# settings for A4 
layout(matrix(c(1:6),2,2,byrow =TRUE))
par(mar=c(5,4,4,4))
# ELF1
boxplot(summary$ELF1_mean.ct ~ summary$treatment,las=2,ylab="ct", main="ELF1", ylim=c(0,40))
stripchart(summary$ELF1_mean.ct ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)
# PR1
boxplot(summary$PR1_mean.ct ~ summary$treatment,las=2,ylab="ct", main="PR1", ylim=c(0,40))
stripchart(summary$PR1_mean.ct ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)
# PG4
boxplot(summary$PG4_mean.ct ~ summary$treatment,las=2,ylab="ct", main="PG4", ylim=c(0,40))
stripchart(summary$PG4_mean.ct ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)
# HCT1
boxplot(summary$HCT1_mean.ct ~ summary$treatment,las=2,ylab="ct", main="HCT1", ylim=c(0,40))
stripchart(summary$HCT1_mean.ct ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)
# PRX2
boxplot(summary$PRX2_mean.ct ~ summary$treatment,las=2,ylab="ct", main="PRX2", ylim=c(0,40))
stripchart(summary$PRX2_mean.ct ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)
# CDH
boxplot(summary$CDH_mean.ct ~ summary$treatment,las=2,ylab="ct", main="CDH", ylim=c(0,40))
stripchart(summary$CDH_mean.ct ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)

# boxplot with log2.NRQ for each gene
# PR1
boxplot(summary$PR1_log2.NRQ ~ summary$treatment,las=2,ylab="log2NRQ", main="PR1")
stripchart(summary$PR1_log2.NRQ ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)
# PG4
boxplot(summary$PG4_log2.NRQ ~ summary$treatment,las=2,ylab="log2NRQ", main="PG4")
stripchart(summary$PG4_log2.NRQ ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)
# HCT1
boxplot(summary$HCT1_log2.NRQ ~ summary$treatment,las=2,ylab="log2NRQ", main="HCT1")
stripchart(summary$HCT1_log2.NRQ ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)
# PRX2
boxplot(summary$PRX2_log2.NRQ ~ summary$treatment,las=2,ylab="log2NRQ", main="PRX2")
stripchart(summary$PRX2_log2.NRQ ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)
# CDH
boxplot(summary$CDH_log2.NRQ ~ summary$treatment,las=2,ylab="log2NRQ", main="CDH")
stripchart(summary$CDH_log2.NRQ ~ summary$treatment, col="orangered", vertical = TRUE, add=T, pch=16)

# ------------------------------------------------------
# VISUALIZATION: BARPLOTS
# ------------------------------------------------------
# preperations for barplots: get mean values and standard error values for each treatment
# PR1
mean.PR1 <- sapply(split(summary$PR1_log2.NRQ, summary$treatment), mean, na.rm=T)
se.PR1 <- sapply(split(summary$PR1_log2.NRQ, summary$treatment), std.err)
# PG4
mean.PG4 <- sapply(split(summary$PG4_log2.NRQ, summary$treatment), mean, na.rm=T)
se.PG4 <- sapply(split(summary$PG4_log2.NRQ, summary$treatment), std.err)
# HCT1
mean.HCT1 <- sapply(split(summary$HCT1_log2.NRQ, summary$treatment), mean, na.rm=T)
se.HCT1<-sapply(split(summary$HCT1_log2.NRQ,summary$treatment), std.err)
# PRX2
mean.PRX2 <- sapply(split(summary$PRX2_log2.NRQ, summary$treatment), mean, na.rm=T)
se.PRX2 <- sapply(split(summary$PRX2_log2.NRQ, summary$treatment), std.err)
# CDH
mean.CDH <- sapply(split(summary$CDH_log2.NRQ, summary$treatment), mean, na.rm=T)
se.CDH <- sapply(split(summary$CDH_log2.NRQ, summary$treatment),std.err)

# barplot: generate a barplot and save it directly in a svg file
# settings for the barplots
angle <- c(0.6, 1.7, 2.8, 3.9)
mycol <- c("gray100", "#ffaaaaff", "#333333ff", "#ff0000ff")
par(mar=c(0.25, 2.5, 0.25, 0.25), cex=1, mfrow=c(1,1))

# PRX2
svg(file="PRX2.svg", width=3, height=3)
barplot(mean.PRX2, 
        xlab="", ylab="", main="",
        ylim=c(-2,2),
        col=mycol,
        xaxt='n',
        las=2,
        space=0.1,
        lwd=0.5,
        tck=-0.05)
MP<-mean.PRX2 + se.PRX2
MM<-mean.PRX2 - se.PRX2
segments(angle, MP, angle, MM)
segments(angle-0.1, MP, angle+0.1, MP)
segments(angle-0.1, MM, angle+0.1, MM)
dev.off()

# PR1
svg(file="PR1.svg", width=3, height=3)
barplot(mean.PR1, 
        xlab="", ylab="",main="",
        ylim=c(-2,2),
        col=mycol,
        xaxt='n',
        las=2,
        space=0.1,
        lwd=0.5,
        tck=-0.05)
MP <- mean.PR1 + se.PR1
MM <- mean.PR1 - se.PR1
segments(angle, MP, angle, MM)
segments(angle-0.1, MP, angle+0.1, MP)
segments(angle-0.1, MM, angle+0.1, MM)
dev.off()

# HCT1
svg(file="HCT1.svg", width=3, height=3)
barplot(mean.HCT1, 
        xlab="", ylab="",main="",
        ylim=c(-2,2),
        col=mycol,
        xaxt='n',
        las=2,
        space=0.1,
        lwd=0.5,
        tck=-0.05)
MP <- mean.HCT1 + se.HCT1
MM <- mean.HCT1 - se.HCT1
segments(angle, MP, angle, MM)
segments(angle-0.1, MP, angle+0.1, MP)
segments(angle-0.1, MM, angle+0.1, MM)
dev.off()

# PG4
svg(file="PG4.svg", width=3, height=3)
barplot(mean.PG4, 
        xlab="", ylab="",main="",
        ylim=c(-2,2),
        col=mycol,
        xaxt='n',
        las=2,
        space=0.1,
        lwd=0.5,
        tck=-0.05)
MP <- mean.PG4 + se.PG4
MM <- mean.PG4 - se.PG4
segments(angle, MP, angle, MM)
segments(angle-0.1, MP, angle+0.1, MP)
segments(angle-0.1, MM, angle+0.1, MM)
dev.off()

# CDH
svg(file="CDH.svg", width=3, height=3)
barplot(mean.CDH, 
        xlab="", ylab="",main="",
        ylim=c(-2,4),
        col=mycol,
        xaxt='n',
        las=2,
        space=0.1,
        lwd=0.5,
        tck=-0.05)
MP<-mean.CDH+se.CDH
MM<-mean.CDH-se.CDH
segments(angle, MP, angle,MM)
segments(angle-0.1, MP, angle+0.1, MP)
segments(angle-0.1, MM, angle+0.1, MM)
dev.off()

# ------------------------------------------------------
# VISUALIZATION: HEATMAP
# ------------------------------------------------------
# generate a new empty dataframe for mean values 
mean <- matrix(data=NA, nrow=5, ncol=length(unique(summary$treatment)), byrow=FALSE, dimnames=NULL)

# fill the mean dataframe with mean values for each gene and each treatment 
for(z in (1:length(unique(summary$treatment)))) {
  x<-subset(summary,summary$treatment==unique(summary$treatment)[z])
  mean[1,z] <- mean(x$CDH_log2.NRQ, na.rm=T)
  mean[2,z] <- mean(x$HCT1_log2.NRQ, na.rm=T)
  mean[3,z] <- mean(x$PG4_log2.NRQ, na.rm=T)
  mean[4,z] <- mean(x$PR1_log2.NRQ, na.rm=T)
  mean[5,z] <- mean(x$PRX2_log2.NRQ, na.rm=T)
}

# get row and column names
rownames(mean) <- c("CDH","HCT1","PG4","PR1","PRX2")
colnames(mean) <- c("C","P","T","PT")

# make heatmap function
make.heatmap1<-function(x){Heatmap(x,
                                   cluster_columns=FALSE,
                                   cluster_rows=FALSE,
                                   col=colorRamp2(c(-2.5,0,4),
                                                  c("royalblue3", "white", "red4")),
                                   row_names_side="left", na_col="gray30",
                                   clustering_distance_rows="none",
                                   clustering_method_rows="none",
                                   km=0, show_heatmap_legend=TRUE, 
                                   use_raster=FALSE,
                                   column_title_side="top",
                                   column_title_gp=gpar(fontsize=5, 
                                                        fontface="bold"),
                                   row_names_gp=gpar(fontsize=10), 
                                   rect_gp=gpar(col="black"),
                                   heatmap_legend_param=list 
                                   (title="relative position"),
                                   row_names_max_width=unit(100,"mm"))
}

# generate the heatmap based on the dataframe with the mean values
svg(file="heatmap_L0.svg", 
    width=2.5, 
    height=2.5)
make.heatmap1(mean)
dev.off()
