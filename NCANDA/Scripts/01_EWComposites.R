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
coglongdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20191115.outlierremoved.csv")
####vars############
####################
###Exclude crystalized intelligence measures from composites 

accvars<-c("cnp_cpf_ifac_tot","cnp_cpw_iwrd_tot","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr","cnp_pmat24a_pmat24_a_cr","cnp_cpfd_dfac_tot","cnp_cpwd_dwrd_tot",
           "cnp_shortvolt_svt","cnp_er40d_er40_cr","cnp_pcet_pcet_acc2","cnp_medf36_medf36_a","cnp_pvrt_pvrt_pc","cnp_svdelay_svt_ld","latentdd")   

latvars<-c("cnp_cpf_ifac_rtc","cnp_cpw_iwrd_rtc","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","cnp_pmat24a_pmat24_a_rtcr","cnp_cpfd_dfac_rtc",
           "cnp_cpwd_dwrd_rtc","cnp_shortvolt_svtcrt","cnp_er40d_er40_crt","cnp_pcet_pcetrtcr","cnp_medf36_medf36_t","cnp_pvrt_pvrtrtcr","cnp_svdelay_svtldrtc","stroop_total_mean","latentgroove")

allfactorvars<-c(accvars,latvars)
####################
###composites#######
coglongdata$Latencycomposite<-compositecols(latvars,coglongdata)
coglongdata$Accuracycomposite<-compositecols(accvars,coglongdata)
###factors#######
latencyfa<-psych::fa(coglongdata[,latvars],nfactors=1,fm="ml",missing=TRUE,impute="median",scores="tenBerge")
coglongdata$Latencyfactorscore<-latencyfa$scores
accuracyfa<-psych::fa(coglongdata[,accvars],nfactors=1,fm="ml",missing=TRUE,impute="median",scores="tenBerge")
coglongdata$Accuracyfactorscores<-accuracyfa$scores
write.csv(coglongdata,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20191115.outlierremoved.compositeacclat.csv")




