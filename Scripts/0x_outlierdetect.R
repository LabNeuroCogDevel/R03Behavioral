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
######functions#####
windzthrehold<-function(x,sdcutoff=5){
  x<-LNCDR::zscore(x)
  x[abs(x)>sdcutoff]<-NA
  return(x)
}
gamleveragedetect<-function(outcome,covar,idvar=NULL,sdcutoff=4){
  #this function does leverage statitiscs based on gam residuals
  outcome<-unlist(outcome)
  covar<-unlist(covar)
  idvar<-as.factor(idvar)
  model<-"outcome~s(covar)"
  if (!is.null(idvar)){
    model<-"outcome~s(covar)+s(idvar,bs='re')"
  }
  f<-as.formula(model)
  m<-mgcv::gam(f)
  zresid<-scale(m$residuals)
  outputoutcome<-outcome
  outputoutcome[abs(zresid)>sdcutoff]<-NA
  return(outputoutcome)
}
distributeNAsfrommeans<-function(meandf,longdf){
  substocensor<-meandf[!complete.cases(meandf),]
  for (i in 1:nrow(substocensor)){
    subi<-substocensor[i,]
    NAcols<-names(subi)[which(is.na(subi))]
    longdf[longdf$id==subi$id,NAcols]<-NA
    print(subi$id)
    print(length(which(longdf$id==subi$id)))
  }
  return(longdf)
}

############
coglongdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190313.csv")
####vars############
####################
allfactorvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","VGS_CRLat","SOC.Problems.solved.in.minimum.moves",
                 "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
                 "nfixbreak","best_acc_m_exclude","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.Percent.correct","SSP.Span.length","nfixbreak")
accvarsperc<-c("Anti_CRR","Mix_CRR","DMS.PC")
accvarscount<-c("SSP.Span.length","nfixbreak")
##convert DMS to perc
coglongdata$DMS.PC<-coglongdata$DMS.Percent.correct/100
###add visit column###
coglongdata<-coglongdata %>% group_by(id) %>% mutate(visit=rank(d8))
coglongdata<-coglongdata[coglongdata$id!=10408,]
####################
coglongwinds<-coglongdata
###visits#######
coglongwinds[,allfactorvars]<-lapply(coglongwinds[,allfactorvars],gamleveragedetect,covar=coglongwinds$Ageatvisit,idvar=coglongwinds$id)
###across all measures 30 visits included one measures that was removed from analysis 
mvouts<-scale(psych::outlier(coglongwinds[,allfactorvars]))
coglongmdcensor<-coglongwinds
coglongmdcensor$mvout<-0
coglongmdcensor$mvout[mvouts>4]<-1
coglongmdcensor<-coglongmdcensor[coglongmdcensor$mvout==0,]
###7 visits removed entierly#########
####invidiaul subjects########
coglongwinds_means<-coglongmdcensor %>% dplyr::select(c(allfactorvars,"id")) %>% group_by(id) %>% summarise_all(mean)
coglongwinds_means[,allfactorvars]<-abs(scale(coglongwinds_means[,allfactorvars]))
coglongwinds_means_vars<-coglongwinds_means[,allfactorvars]
coglongwinds_means_vars[coglongwinds_means_vars>4]<-NA
coglongwinds_means_varswitid<-coglongwinds_means_vars
coglongwinds_means_varswitid$id<-coglongwinds_means$id

###
coglongmdcensorwithsubcensor<-coglongmdcensor
coglongmdcensorwithsubcensor<-distributeNAsfrommeans(coglongwinds_means_varswitid,coglongmdcensorwithsubcensor)
###642 complete visits 670 total visits#####
##########write data########
#write.csv(coglongmdcensor,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190918.outlierremoved.csv")
write.csv(coglongmdcensorwithsubcensor,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190918.outlierremoved.csv")









