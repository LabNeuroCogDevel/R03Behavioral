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
#####0X NCDANDA cognitive data#######
###get cognitive data (SAVE) and exploratory factor#######
###WAIS #######
WAIS<-read.csv("/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/wais4.csv")
WAIS_merge<-WAIS[,c("subject","visit","np_wais4_rawscore_computed")]
#######stroop####
stroop<-read.csv("/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/stroop.csv")
stroop_merge<-stroop[,c("subject","visit","stroop_total_mean")]
#####grooved peg board######
groovedpegboard<-read.csv("/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/grooved_pegboard.csv")
groovedpegboard_merge<-groovedpegboard[,c("subject","visit","np_gpeg_dh_time","np_gpeg_ndh_time")]
fgroove<-fa(groovedpegboard_merge[,3:4],nfactors=1,rotate ="bifactor",fm="ml",scores="tenBerge")
groovedpegboard_merge$latentgroove<-fgroove$scores[,1]
####delay discounting#####
delaydiscounting_100<-read.csv("/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/dd100.csv")
delaydiscounting_1k<-read.csv("/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/dd1000.csv")
allddlogk<-merge(delaydiscounting_100,delaydiscounting_1k,by=c("subject","visit"))
allddlogk$ddage<-allddlogk$dd100_age
allddlogk_d<-allddlogk[,c("dd100_logk_1d","dd100_logk_7d","dd100_logk_1mo","dd100_logk_6mo","dd1000_logk_1d","dd1000_logk_7d","dd1000_logk_1mo","dd1000_logk_6mo")]

fdd<-fa(allddlogk_d,nfactors=1,rotate ="bifactor",fm="ml",scores="tenBerge")
#plot(allddlogk$dd100_age,fdd$scores[,1])
allddlogk$latentdd<-fdd$scores[,1]
ddmerge_merge<-allddlogk[,c("subject","visit","latentdd")]
#########CNP##########
cnpdatadict<-read.csv("/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/datadict/redcap/cnp_datadict.csv")
cpndata<-read.csv("/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/NCANDA_RELEASE_4Y_REDCAP_MEASUREMENTS_V02/summaries/redcap/cnp.csv")
cpndata<-merge(cpndata,WAIS_merge,by=c("subject","visit")) %>% merge(stroop_merge,by=c("subject","visit")) %>% merge(groovedpegboard_merge,by=c("subject","visit")) %>%
  merge(ddmerge_merge,by=c("subject","visit"))

efficiency_function<-function(df,vars){
  ####have to flip the sign!
  eff<-scale(df[,vars[1]],center=TRUE,scale=TRUE)+scale(df[,vars[2]],center=TRUE,scale=TRUE)
  print(cor(df[,vars[1]],df[,vars[2]],use="complete"))
  return(eff)
} 
#cpndata<-cpndata[cpndata$visit=="baseline",]
cpf<-c("cnp_cpf_ifac_tot","cnp_cpf_ifac_rtc")
cpndata$cpf_eff<-efficiency_function(cpndata,cpf)
cpw<-c("cnp_cpw_iwrd_tot","cnp_cpw_iwrd_rtc")
cpndata$cpw_eff<-efficiency_function(cpndata,cpw)
spctnl<-c("cnp_spcptnl_scpt_tp","cnp_spcptnl_scpt_tprt")
cpndata$spcptnl_eff<-efficiency_function(cpndata,spctnl)
sfnb<-c("cnp_sfnb2_sfnb_mcr","cnp_sfnb2_sfnb_mrtc")
cpndata$sfnb_eff<-efficiency_function(cpndata,sfnb)
pmat<-c("cnp_pmat24a_pmat24_a_cr","cnp_pmat24a_pmat24_a_rtcr")
cpndata$pmat_eff<-efficiency_function(cpndata,pmat)
cpfd<-c("cnp_cpfd_dfac_tot","cnp_cpfd_dfac_rtc")
cpndata$cpfd_eff<-efficiency_function(cpndata,cpfd)
cpwd<-c("cnp_cpwd_dwrd_tot","cnp_cpwd_dwrd_rtc")
cpndata$cpwd_eff<-efficiency_function(cpndata,cpwd)
shortvolt<-c("cnp_shortvolt_svt","cnp_shortvolt_svtcrt")
cpndata$shortvolt_eff<-efficiency_function(cpndata,shortvolt)
er40d<-c("cnp_er40d_er40_cr","cnp_er40d_er40_crt")
cpndata$er40d_eff<-efficiency_function(cpndata,er40d)
pcet<-c("cnp_pcet_pcet_acc2","cnp_pcet_pcetrtcr")
cpndata$pcet_eff<-efficiency_function(cpndata,pcet)
medf36<-c("cnp_medf36_medf36_a","cnp_medf36_medf36_t")
cpndata$medf36_eff<-efficiency_function(cpndata,medf36)
pvoc<-c("cnp_pvoc_pvoccr","cnp_pvoc_pvocrtcr")
cpndata$pvoc_eff<-efficiency_function(cpndata,pvoc)
pvrt<-c("cnp_pvrt_pvrt_pc","cnp_pvrt_pvrtrtcr")
cpndata$pvrt_eff<-efficiency_function(cpndata,pvrt)
svdelay<-c("cnp_svdelay_svt_ld","cnp_svdelay_svtldrtc")
cpndata$svdelay_eff<-efficiency_function(cpndata,svdelay)
allneuropsychvars<-c(cpf,cpw,spctnl,sfnb,pmat,cpfd,cpwd,
                     shortvolt,er40d,pcet,medf36,pvoc,pvrt,svdelay)
allneuropsychvarpairs<-c("cpf","cpw","spctnl","sfnb","pmat","cpfd","cpwd",
                         "shortvolt","er40d","pcet","medf36","pvoc","pvrt","svdelay")
allneuropsychvarpairs_eff<-c("cpf_eff","cpw_eff","spcptnl_eff","sfnb_eff","pmat_eff","cpfd_eff",
                             "cpwd_eff","shortvolt_eff","er40d_eff","pcet_eff","medf36_eff","pvoc_eff",
                             "pvrt_eff","svdelay_eff")
allnueropsychaccvars<-c(cpf[1],cpw[1],spctnl[1],sfnb[1],pmat[1],cpfd[1],cpwd[1],shortvolt[1],er40d[1],pcet[1],medf36[1],pvoc[1],pvrt[1],svdelay[1])
allnueropsychrtvars<-c(cpf[2],cpw[2],spctnl[2],sfnb[2],pmat[2],cpfd[2],cpwd[2],shortvolt[2],er40d[2],pcet[2],medf36[2],pvoc[2],pvrt[2],svdelay[2])

nocnpvars<-c("np_wais4_rawscore_computed","stroop_total_mean","latentgroove","latentdd")

neuropsychvars<-c(allneuropsychvars,nocnpvars)

cpndatatosave<-cpndata
cpndatatosave<-cpndatatosave[,c("subject","visit","cnp_age",neuropsychvars)]
write.csv(cpndatatosave,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/NCANDA/Data/btc_NCANDAscoredmeasures_20191115.csv")










####outlier detect
outiersettoNA<-function(x,sd=4){
  xz<-scale(x,center=TRUE,scale=TRUE)
  xout<-ifelse(abs(xz)>sd,NA,x)
  return(xout)
}
####outlierremove
cpndata<-mutate_at(cpndata,neuropsychvars,outiersettoNA)
agemerge<-cpndata[,c("subject","visit","cnp_age")]
#####factor######
library(psych)
library(pan)
library(GPArotation)
library(nFactors)
#######################
cpndata_reduce<-cpndata[,c("subject","visit",as.character(neuropsychvars))]
cpndata_reduce_baseline<-cpndata_reduce[cpndata_reduce$visit=="baseline",]
cpndata_reduce_baseline_complete<-cpndata_reduce_baseline[complete.cases(cpndata_reduce_baseline),]
cpndata_reduce_baseline_complete_d<-cpndata_reduce_baseline_complete[,neuropsychvars]
cpndata_reduce_baseline_complete_d_scale<-as.data.frame(scale(cpndata_reduce_baseline_complete_d,center=TRUE,scale=TRUE))
#ncnpfactors<-psych::fa.parallel(cpndata_reduce_baseline_complete_d_scale,n.iter=1000,fm="ml",fa="fa")
f<-fa(cpndata_reduce_baseline_complete_d_scale,nfactors=8,rotate ="bifactor",fm="ml",scores="tenBerge")
fscores_df<-as.data.frame(f$scores)
fscores_df$sub<-cpndata_reduce_baseline_complete$subject
write.table(fscores_df,"/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/CNP data/cnpfactorscores_baseline.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
########all visits at one time######
cpndata_reduce_full_complete<-cpndata_reduce[complete.cases(cpndata_reduce),]
cpndata_reduce_full_complete_d<-cpndata_reduce_full_complete[,neuropsychvars]
cpndata_reduce_full_complete_d_scale<-as.data.frame(scale(cpndata_reduce_full_complete_d,center=TRUE,scale=TRUE))
#ncnpfactors<-psych::fa.parallel(cpndata_reduce_full_complete_d_scale,n.iter=1000,fm="ml",fa="fa")
ffull<-fa(cpndata_reduce_full_complete_d_scale,nfactors=8,rotate ="bifactor",fm="ml",scores="tenBerge")
fscoresfull_df<-as.data.frame(ffull$scores)
fscoresfull_df$subject<-cpndata_reduce_full_complete$subject
fscoresfull_df$visit<-cpndata_reduce_full_complete$visit
fscoresfull_dfwithage<-merge(fscoresfull_df,agemerge,by=c("subject","visit"))
names(fscoresfull_df)[1]<-paste0("general_rt",names(fscoresfull_df)[1])
names(fscoresfull_df)[2]<-paste0("general_accuracy",names(fscoresfull_df)[2])
write.csv(fscoresfull_df,"/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/CNP data/cnpfactorscores_allvisits.csv")
#####################################
cpndata_reduce_followup1y<-cpndata_reduce[cpndata_reduce$visit=="followup_1y",]
cpndata_reduce_followup1y_complete<-cpndata_reduce_followup1y[complete.cases(cpndata_reduce_followup1y),]
cpndata_reduce_followup1y_complete_d<-cpndata_reduce_followup1y_complete[,neuropsychvars]
cpndata_reduce_followup1y_complete_d_scale<-as.data.frame(scale(cpndata_reduce_followup1y_complete_d,center=TRUE,scale=TRUE))

#ncnpfactors<-psych::fa.parallel(cpndata_reduce_followup1y_complete_d_scale,n.iter=1000,fm="ml",fa="fa")
fffollowup1y<-fa(cpndata_reduce_followup1y_complete_d_scale,nfactors=8,rotate ="bifactor",fm="ml",scores="tenBerge")
fscoresfollowup1y_df<-as.data.frame(fffollowup1y$scores)
fscoresfollowup1y_dff$sub<-cpndata_reduce_followup1y_complete$subject
#write.table(fscoresfollowup1y_dff,"/Users/brendenclemmens/Desktop/Projects/NCANDA_Behavioral/Data/CNP data/cnpfactorscores_baseline.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

cor(f$loadings,fffollowup1y$loadings)
cor(ffull$loadings,f$loadings)
cor(ffull$loadings,fffollowup1y$loadings)
###m Ml 2 and 4 questionable#####

########################split by age################
cpndata_reduce_withage<-cpndata[,c("subject","visit","cnp_age",as.character(neuropsychvars))]
cpndata_reduce_withage_complete<-cpndata_reduce_withage[complete.cases(cpndata_reduce_withage),]
########below 18
cpndata_reduce_withage_complete_lessthan18<-cpndata_reduce_withage_complete[cpndata_reduce_withage_complete$cnp_age<=18,]

cpndata_reduce_withage_complete_lessthan18_d<-cpndata_reduce_withage_complete_lessthan18[,neuropsychvars]
cpndata_reduce_withage_complete_lessthan18_d_scale<-as.data.frame(scale(cpndata_reduce_withage_complete_lessthan18_d,center=TRUE,scale=TRUE))
#ncnpfactors_lessthan18<-psych::fa.parallel(cpndata_reduce_withage_complete_lessthan18_d_scale,n.iter=1000,fm="ml",fa="fa")
falessthan18<-fa(cpndata_reduce_withage_complete_lessthan18_d_scale,nfactors=8,rotate ="bifactor",fm="ml",scores="tenBerge")

cor(ffull$loadings,falessthan18$loadings)
######above 18
cpndata_reduce_withage_complete_greaterthan18<-cpndata_reduce_withage_complete[cpndata_reduce_withage_complete$cnp_age>18,]
cpndata_reduce_withage_complete_greaterthan18_d<-cpndata_reduce_withage_complete_greaterthan18[,neuropsychvars]
cpndata_reduce_withage_complete_greaterthan18_d_scale<-as.data.frame(scale(cpndata_reduce_withage_complete_greaterthan18_d,center=TRUE,scale=TRUE))
#ncnpfactors_greaterthan18<-psych::fa.parallel(cpndata_reduce_withage_complete_greaterthan18_d_scale,n.iter=1000,fm="ml",fa="fa")
fagreaterthan18<-fa(cpndata_reduce_withage_complete_greaterthan18_d_scale,nfactors=8,rotate ="bifactor",fm="ml",scores="tenBerge")

cor(ffull$loadings,fagreaterthan18$loadings)
cor(falessthan18$loadings,fagreaterthan18$loadings)





allneuropsychdata_dn_scale<-

allneuropsychdata_dn_scale_complete<-allneuropsychdata_dn_scale[complete.cases(allneuropsychdata_dn_scale),]
completecasesubjects<-allneuropsychdata$subjectkey[complete.cases(allneuropsychdata_dn_scale)]
nfactorspsych<-psych::fa.parallel(allneuropsychdata_dn_scale_complete,n.iter=1000,fm="ml",fa="fa")



