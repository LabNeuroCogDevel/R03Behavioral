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
#######
gampredvisit<-function(df,outcomevars,visitvar,covar,idvar=NULL){
  visitpredout<-NULL
  for (ov in outcomevars){
    print(ov)
  outcome<-df[,ov]
  outcome<-scale(outcome)
  covard<-df[,covar]
  visitvard<-df[,visitvar]
  idvard<-df[,idvar]
  model<-"outcome~s(covard)+s(visitvard)"
  if (!is.null(idvard)){
    model<-"outcome~s(covard)+s(visitvard)+s(idvard,bs='re')"
  }
  f<-as.formula(model)
  m<-mgcv::gam(f,optimizer="efs")
  print(summary(m)$s.pv[2])
  interval_inc<-.1
  v <- m$model[, "visitvard"]
  cond_list <- list(seq(min(v), max(v), by=interval_inc))
  names(cond_list) <- "visitvard"
  visitpred <- itsadug::get_predictions(m, cond = cond_list,rm.ranef = TRUE)
  visitpred$var<-ov
  visitpredi<-visitpred[,c("visitvard","fit","CI","var")]
  visitpredout<-rbind(visitpredout,visitpredi)
  }
  return(visitpredout)
}
############
coglongdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190918.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$id)
####vars############
####################
allfactorvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","VGS_CRLat","SOC.Problems.solved.in.minimum.moves",
                 "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
                 "nfixbreak_fl","best_acc_m_exclude_fl","first_lat_m_exclude","Latencycomposite","Accuracycomposite","Accuracyfactorscores","Latencyfactorscore")

accvars<-c("Anti_CRR","Mix_CRR","DMS.Percent.correct","SSP.Span.length","nfixbreak_fl","best_acc_m_exclude_fl")

latencyvars<-c("Anti_CRLat","Mix_CRLat","VGS_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency","first_lat_m_exclude","Latencycomposite")
latencyvarsplot<-c("Latencycomposite","Anti_CRLat","first_lat_m_exclude","Mix_CRLat","VGS_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency")
accvarsplot<-c("Accuracycomposite","Anti_CRR","best_acc_m_exclude_fl","Mix_CRR","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","SSP.Span.length")
################
accvisitpredsz<-gampredvisit(coglongdata,accvarsplot,"visit","Ageatvisit","id")
accvisitplot <- ggplot(accvisitpredsz,aes(x = visitvard, y = fit,colour=var,fill=var)) + 
  geom_line(size = 2) + ylab("") + 
  xlab("")
accvisitplot_luna<-lunaize(accvisitplot)

latvisitpredsz<-gampredvisit(coglongdata,latencyvarsplot,"visit","Ageatvisit","id")
latvisitplot <- ggplot(latvisitpredsz,aes(x = visitvard, y = fit,colour=var,fill=var)) + 
  geom_line(size = 2) + ylab("") + xlab("") #+ geom_ribbon(aes(ymin = fit - CI, ymax = fit + CI),alpha = 0.3)
  
latvisitplot_luna<-lunaize(latvisitplot)








+geom_ribbon(aes(ymin = fit - CI, ymax = fit + CI),alpha = 0.3)
+
  
  geom_point(data = df, aes(y = outcome, x = pred), alpha = 0.2)+
  geom_line(data = df, aes(y = outcome, group = id), alpha = 0.2)+scale_color_manual(values = c("white","#F8766D","#00BFC4"))+
  scale_fill_manual(values=c("white","#F8766D","#00BFC4"))




