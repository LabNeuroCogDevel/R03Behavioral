library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(lsmeans)
library(mgcv)
library(itsadug)
library(lme4)
library(lsmeans)
library(stats)
library(psych)
library(LNCDR)
library(FactoMineR)
library(corrplot)
library(mgcv)
############
compositecols<-function(cols,data){
  #return composite score (z scored) of data[,cols]
  #seperate function to increase flexibilit
  data_z<-scale(data[,cols])
  compositeout<-scale(base::rowSums(data_z,na.rm=TRUE))
  return(compositeout)
}

########
coglongdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190918.outlierremoved.csv")
####vars############
####################
allfactorvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","VGS_CRLat","SOC.Problems.solved.in.minimum.moves",
                 "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
                 "nfixbreak_fl","best_acc_m_exclude_fl","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.Percent.correct","SSP.Span.length","nfixbreak_fl","best_acc_m_exclude_fl")

latencyvars<-c("Anti_CRLat","Mix_CRLat","VGS_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency","first_lat_m_exclude")
##convert DMS to perc
coglongdata$DMS.PC<-coglongdata$DMS.Percent.correct/100
###covert acc measures to be positive=better####
coglongdata$best_acc_m_exclude_fl<-coglongdata$best_acc_m_exclude*-1
coglongdata$nfixbreak_fl<-coglongdata$nfixbreak*-1
####################
###composites#######
coglongdata$Latencycomposite<-compositecols(latencyvars,coglongdata)
coglongdata$Accuracycomposite<-compositecols(accvars,coglongdata)
###factors#######
latencyfa<-psych::fa(coglongdata[,latencyvars],nfactors=1,fm="ml",missing=TRUE,impute="median",scores="tenBerge")
coglongdata$Latencyfactorscore<-latencyfa$scores
accuracyfa<-psych::fa(coglongdata[,accvars],nfactors=1,fm="ml",missing=TRUE,impute="median",scores="tenBerge")
coglongdata$Accuracyfactorscores<-accuracyfa$scores
write.csv(coglongdata,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190918.outlierremoved.compositeacclat.csv")




