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
library(lavaan)
#####################
#########function
variancepartitiongam<-function(df,additivevars,basevars,agevar){
  pairs<-as.data.frame(expand.grid(additivevars,basevars))
  names(pairs)<-c("addvar","basevar")
  varout<-NULL
  for (p in 1:nrow(pairs)){
    print(as.character(pairs$addvar[p]))
    df$addvar<-unlist(df[,as.character(pairs$addvar[p])])
    df$basevar<-unlist(df[,as.character(pairs$basevar[p])])
    df$agevar<-df[,agevar]
    modeid<-mgcv::gam(agevar~s(id,bs="re"), data=df,optimizer = "efs")
    modb <- mgcv::gam(agevar~s(basevar,k=4,fx = T)+s(id,bs="re"), data=df,optimizer = "efs")
    moda <- mgcv::gam(agevar~s(basevar,k=4,fx = T)+s(addvar,k=4,fx = T)+s(id,bs="re"), data=df,optimizer = "efs")
    idvar<-summary(modeid)$dev.expl
    totalagevar<-summary(moda)$dev.expl-idvar 
    basevar<-summary(modb)$dev.expl-idvar  
    addvariance<-totalagevar-basevar
    
    basevarout<-c(as.character(pairs$addvar[p]),"common",totalagevar,basevar)
    addvarout<-c(as.character(pairs$addvar[p]),"specific",totalagevar,addvariance)
    varouttemp<-rbind(basevarout,addvarout)
    varout<-rbind(varout,varouttemp)
  }
  return(varout)
}

##########
allcoglongdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/alldatawithwithinbetween.20190906.csv")
allcoglongdata$id<-as.factor(allcoglongdata$id)
###################
allfactorvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","VGS_CRLat","SOC.Problems.solved.in.minimum.moves",
                 "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
                 "nfixbreak_fl","best_acc_m_exclude_fl","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.Percent.correct","SSP.Span.length","nfixbreak")
accvarsperc<-c("Anti_CRR","Mix_CRR","DMS.PC")
accvarscount<-c("SSP.Span.length","nfixbreak")

allfactorvarsorder<-c("AntCRR","MixCRR","nfixbrk","bstc","SOCP","DMSPr","SSPSp",
                      "AntCRL","MxCRLt","VGSCRL","frst","DMSM")
#######################
DGAwithinvarpart<-variancepartitiongam(allcoglongdata,allfactorvarsorder,"DGAwithin","Agetvst")
DGAwithinvarpartdf<-data.frame(DGAwithinvarpart)
names(DGAwithinvarpartdf)<-c("var","commsp","total","relvar")
DGAwithinvarpartdf$relvar<-as.numeric(as.character(DGAwithinvarpartdf$relvar))
DGAwithinvarpartdf$total<-as.numeric(as.character(DGAwithinvarpartdf$total))
DGAwithinvarpartdf$percentvar<-DGAwithinvarpartdf$relvar/DGAwithinvarpartdf$total
DGAwithinvarpartdf$commsp2<-factor(DGAwithinvarpartdf$commsp,levels=c("specific","common"))
gp<-ggplot(DGAwithinvarpartdf,aes(x=var,y=percentvar,fill=commsp))+geom_bar(position="stack", stat="identity")+
  ylab("Age-Related Variance\n")+xlab("")+theme(text=element_text(size=28))+theme(legend.title = element_blank())

ggplot(data, aes(fill=condition, y=percentvar, x=specie)) + 
  
