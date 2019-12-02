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
library(cowplot)
library(tidyr)
##########################
ICCforlongformat<-function(df,vars,subjectvar,visitvar,maxvisits=2){
  #when df is longformat (row=subject at visit) returns ICC values for all measures in vars
  ###only uses complete cases###
  #max visits (how many visits should be useds in ICC); default to use the first two visits
  ###to do...make default max visits the number of unique visits in df? (ASK WF for help on this)
  ICCbind<-NULL
  for (v in vars){
    ####set up data####
    ICCdflong<-df[,c(subjectvar,visitvar,v)]
    ICCdfwide<-tidyr::spread(ICCdflong,visitvar,v)
    ICCdfwidethreshold<-ICCdfwide
    ICCdfwidethreshold[,subjectvar]<-NULL
    ICCdfwidethreshold[,maxvisits+1:ncol(ICCdfwidethreshold)]<-NULL ####clip to max visits
    ICCdfwidethreshold<-ICCdfwidethreshold[complete.cases(ICCdfwidethreshold),]
    ########perform ICC
    ICCsum<-as.data.frame(psych::ICC(ICCdfwidethreshold,missing=TRUE)$results)
    ICCwide<-as.data.frame(t(ICCsum$ICC))
    names(ICCwide)<-ICCsum$type
    ###########return details####
    ICCinfo<-data.frame(var=v,observations_used=nrow(ICCdfwidethreshold),visits_used=maxvisits)
    ########bind info
    ICCout<-cbind(ICCwide,ICCinfo)
    #####loop bind
    ICCbind<-rbind(ICCbind,ICCout)
  }
  return(ICCbind)
}
#################read data######
coglongdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190918.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$id)
####vars#########
allfactorvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","VGS_CRLat","SOC.Problems.solved.in.minimum.moves",
                 "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
                 "nfixbreak_fl","best_acc_m_exclude_fl","first_lat_m_exclude","Latencycomposite","Accuracycomposite","Accuracyfactorscores","Latencyfactorscore")

accvars<-c("Anti_CRR","Mix_CRR","DMS.Percent.correct","SSP.Span.length","nfixbreak_fl","best_acc_m_exclude_fl")

latencyvars<-c("Anti_CRLat","Mix_CRLat","VGS_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency","first_lat_m_exclude","Latencycomposite")
latencyvarsplot<-c("Latencycomposite","Anti_CRLat","first_lat_m_exclude","Mix_CRLat","VGS_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency")
accvarsplot<-c("Accuracycomposite","Anti_CRR","best_acc_m_exclude_fl","Mix_CRR","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","SSP.Span.length")
####ICC########
ICCS<-ICCforlongformat(coglongdata,allfactorvars,"id","visit")
write.csv(ICCS,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/ICCs.csv")
###############################################variable correlations############
#visit 1 only
coglongdata_v1<-coglongdata[coglongdata$visit==1,]
#no age removal####
zeroordercors<-cor(coglongdata[,allfactorvars],use="pairwise.complete")
write.csv(zeroordercors,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/zeroordercors.csv")
####age removal#####



