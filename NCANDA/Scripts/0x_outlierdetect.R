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
library(gamm4)
######functions#####
windzthrehold<-function(x,sdcutoff=5){
  x<-LNCDR::zscore(x)
  x[abs(x)>sdcutoff]<-NA
  return(x)
}
gamm4leveragedetect<-function(outcome,covar,idvar=NULL,sdcutoff=4){
  #this function does leverage statitiscs based on gam residuals
  outcome<-unlist(outcome)
  covar<-unlist(covar)
  idvar<-as.factor(idvar)
  if (!is.null(idvar)){
    m<-gamm4(outcome~s(covar,k=4),random=~(1|idvar))
  }else{
    m<-gamm4(outcome~s(covar,k=4))
  }
  zresid<-scale(m$gam$residuals)
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
NCANDAdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20191115.csv")
####vars############
####################
accvars<-c("cnp_cpf_ifac_tot","cnp_cpw_iwrd_tot","cnp_spcptnl_scpt_tp","cnp_sfnb2_sfnb_mcr","cnp_pmat24a_pmat24_a_cr","cnp_cpfd_dfac_tot","cnp_cpwd_dwrd_tot",
           "cnp_shortvolt_svt","cnp_er40d_er40_cr","cnp_pcet_pcet_acc2","cnp_medf36_medf36_a","cnp_pvoc_pvoccr","cnp_pvrt_pvrt_pc","cnp_svdelay_svt_ld","np_wais4_rawscore_computed","latentdd")   

latvars<-c("cnp_cpf_ifac_rtc","cnp_cpw_iwrd_rtc","cnp_spcptnl_scpt_tprt","cnp_sfnb2_sfnb_mrtc","cnp_pmat24a_pmat24_a_rtcr","cnp_cpfd_dfac_rtc",
          "cnp_cpwd_dwrd_rtc","cnp_shortvolt_svtcrt","cnp_er40d_er40_crt","cnp_pcet_pcetrtcr","cnp_medf36_medf36_t","cnp_pvoc_pvocrtcr","cnp_pvrt_pvrtrtcr","cnp_svdelay_svtldrtc","stroop_total_mean","latentgroove")

allfactorvars<-c(accvars,latvars)
####################
coglongwinds<-NCANDAdata
coglongwinds<-coglongwinds[!is.na(coglongwinds$cnp_age),]
###visits#######
coglongwinds[,allfactorvars]<-lapply(coglongwinds[,allfactorvars],gamm4leveragedetect,covar=coglongwinds$cnp_age,idvar=coglongwinds$subject)
###determine missing visits 20191115#####
mvouts<-scale(psych::outlier(coglongwinds[,allfactorvars]))
coglongmdcensor<-coglongwinds
coglongmdcensor$mvout<-0
coglongmdcensor$mvout[mvouts>4]<-1
coglongmdcensor<-coglongmdcensor[coglongmdcensor$mvout==0,]
###35 visits removed entierly#########
####invidiaul subjects########
coglongwinds_means<-coglongmdcensor %>% dplyr::select(c(allfactorvars,"subject")) %>% group_by(subject) %>% summarise_all(mean)
coglongwinds_means[,allfactorvars]<-abs(scale(coglongwinds_means[,allfactorvars]))
coglongwinds_means_vars<-coglongwinds_means[,allfactorvars]
coglongwinds_means_vars[coglongwinds_means_vars>4]<-NA
coglongwinds_means_varswitid<-coglongwinds_means_vars
coglongwinds_means_varswitid$id<-coglongwinds_means$id

###
coglongmdcensorwithsubcensor<-coglongmdcensor
coglongmdcensorwithsubcensor<-distributeNAsfrommeans(coglongwinds_means_varswitid,coglongmdcensorwithsubcensor)
###2591 complete visits 3399 total visits#####
##########write data########
write.csv(coglongmdcensorwithsubcensor,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20191115.outlierremoved.csv")









