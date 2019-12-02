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
#####CANTAB############################
CANTAB<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/CogCANTAB_20190128.csv")
CANTAB<-CANTAB[CANTAB$Age>7,]
CANTAB<-CANTAB[CANTAB$Age<35,]
CANTAB_nowarn<-CANTAB[,!(names(CANTAB) %in% grep("Warning",names(CANTAB),value=TRUE))]
CANTAB_nowarn$SOC.Overallmeaninitialthinkingtime<- apply(CANTAB_nowarn[c("SOC.Mean.initial.thinking.time..2.moves.","SOC.Mean.initial.thinking.time..3.moves.","SOC.Mean.initial.thinking.time..4.moves.","SOC.Mean.initial.thinking.time..5.moves.")], 1, median)
CANTAB_nowarn$SOC.Overallmeansubsequentthinkingtime<-apply(CANTAB_nowarn[c("SOC.Mean.subsequent.thinking.time..2.moves.","SOC.Mean.subsequent.thinking.time..3.moves.","SOC.Mean.subsequent.thinking.time..4.moves.","SOC.Mean.subsequent.thinking.time..5.moves.")],1, median)
CANTAB_nowarn$d8<-as.numeric(format(mdy_hms(CANTAB_nowarn$Session.start.time),"%Y%m%d"))
CANTAB_nowarn$id<-gsub("_.","",gsub("[xX]([0-9])","",CANTAB_nowarn$Subject.ID))
###gender write from CANTAB (Double check this 20190825)
GENDER_R<-CANTAB_nowarn[,c("id","Gender")]
write.csv(GENDER_R,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/gender_20190825.csv")
#####
#MOT Median Latency MOT Mean Error###
###SOC.Problems.solved.in.minimum.moves" #
###MEan SOC.Mean.initial.thinking.time..3.moves###
#"DMS.Median.correct.latency"    
#DMS.Percent.correct
#SSP.Span.length

###Measures
CANTAB_Tmeasures<-CANTAB_nowarn[,c("id","Age","d8","MOT.Median.latency","MOT.Mean.error","SOC.Problems.solved.in.minimum.moves","SOC.Overallmeaninitialthinkingtime","SOC.Overallmeansubsequentthinkingtime","DMS.Percent.correct","DMS.Median.correct.latency","SSP.Span.length","SSP.Mean.time.to.first.response","SSP.Mean.time.to.last.response")]
cor(CANTAB_Tmeasures,use="complete.obs",method ="kendall")
ggplot(CANTAB_Tmeasures, aes(Age, SSP.Span.length)) + 
  geom_point() +
  geom_smooth(method = "loess", se = FALSE)

CANTABdemo<-CANTAB_nowarn[,c("Age","SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","SSP.Span.length")]
CANTABdemo_d<-CANTABdemo[,c("SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","SSP.Span.length")]
#ncnpfactors<-psych::fa.parallel(CANTABdemo_d,n.iter=1000,fm="ml",fa="fa")
f<-fa(CANTABdemo_d,nfactors=1,rotate ="bifactor",fm="ml",scores="tenBerge")
CANTABdemo$scores<-f$scores

ggplot(CANTABdemo, aes(Age, scores)) + 
  geom_point() +
  geom_smooth(method = "loess", se = FALSE)

###########################
####MGS scored#####
MGS<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/mgs_btc_summaries_20190313.csv")
MGS$id<-MGS$subj
MGS$d8<-MGS$date
###Eye Tracking ############
library(dplyr)
library(tidyr)
Eyedata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/all_et_score.csv")
Eyedata<-Eyedata[Eyedata$lat>60,]
Eyedata<-Eyedata[!is.na(Eyedata$id),]
Eyedata$scoringtype<-"auto"
Eyedata$scoringtype[grep("manual",Eyedata$fn)]<-"manual"
Eyedata_wide <- 
  Eyedata %>%
  filter(Count!=-1) %>% # remove drop
  # summarise
  group_by(id,d8,vtype,task,scoringtype) %>%
  summarise(CRR=length(which(Count==1))/n(),
            CRLat=median(lat[Count==1]),
            Ageatvisit=mean(age),
            nvisits=length(unique(age))) %>% 
  gather(metric,value,  CRR, CRLat) %>%
  unite(task_metric,task,metric) %>%
  spread(task_metric, value )

#Eyedata_wide<-Eyedata_wide[Eyedata_wide$Ageatvisit,]
Eyedata_wide<-Eyedata_wide[Eyedata_wide$Ageatvisit<45,]
Eyedata_wide$antiminusvgslat<-Eyedata_wide$Anti_CRLat-Eyedata_wide$VGS_CRLat
Eyedata_wide$antiminusvgslat_normvgslat<-(Eyedata_wide$Anti_CRLat-Eyedata_wide$VGS_CRLat)/Eyedata_wide$VGS_CRLat
Eyedata_wide$antiminusvgslat_normantiplusvgs<-(Eyedata_wide$Anti_CRLat-Eyedata_wide$VGS_CRLat)/(Eyedata_wide$Anti_CRLat+Eyedata_wide$VGS_CRLat)

Eyedata_wide$iddate<-paste(Eyedata_wide$id,Eyedata_wide$d8,sep="_")
#######check manual and autoscoring#######
a<-reshape2::melt(Eyedata_wide, id.vars=c("iddate","vtype","scoringtype"))
m <-
  a %>% filter(scoringtype=='auto') %>%
  merge(a %>% filter(scoringtype=='manual'),
        by=c("iddate","vtype","variable"),
        suffixes=c('.auto','.man') )
m %>% group_by(variable) %>% summarise(cor(value.auto,value.man, use='pairwise.complete.obs',method="kendall"))

ggplot(m %>% filter(!variable %in% c("id","d8","Ageatvisit","nvisits")),aes(x=value.man,y=value.auto,colour=variable))+geom_point()+facet_wrap(~variable,scales="free")

Eyedata_wide_auto<-m[m$scoringtype=="auto",]
Eyedata_wide_manual<-m[m$scoringtype=="manual",]
eyedata_wide_automanualtest<-merge(m,m,by=c("iddate","vtype"),suffixes=c(".auto",".man"))
#########

# bandedids<-Eyedata_wide$iddate[Eyedata_wide$Anti_CRLat==200]
# bandediddf<-Eyedata_wide[Eyedata_wide$Anti_CRLat==200,]
# write.table(bandedids,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/idstocheck.txt")
# eyedata_bandedids<-Eyedata[Eyedata$id %in% bandediddf$id,]
# eyedata_bandedids<-eyedata_bandedids[eyedata_bandedids$scoringtype=="manual",]

# library(ggplot2)
# ggplot(eyedata_wide_automanualtest, aes(Ageatvisit, Anti_CRR,colour=scoringtype)) + 
#   geom_point() +
#   geom_smooth(method = "loess", se = FALSE)
# 
# ggplot(Eyedata_wide, aes(Ageatvisit, Anti_CRLat)) + 
#   geom_point() +
#   geom_smooth(method = "loess", se = FALSE)
# 
# ggplot(Eyedata_wide, aes(Ageatvisit, VGS_CRLat,colour=scoringtype)) + 
#   geom_point() +
#   geom_smooth(method = "loess", se = FALSE)
# 
# ggplot(Eyedata_wide, aes(Ageatvisit, antiminusvgslat_normvgslat)) + 
#   geom_point() +
#   geom_smooth(method = "loess", se = FALSE)
# 
# ggplot(Eyedata_wide, aes(Ageatvisit, antiminusvgslat_normantiplusvgs)) + 
#   geom_point() +
#   geom_smooth(method = "loess", se = FALSE)
# 
# ggplot(Eyedata_wide, aes(Ageatvisit, antiminusvgslat_normantiplusvgs)) + 
#   geom_point() +
#   geom_smooth(method = "loess", se = FALSE)
# 
# ggplot(Eyedata_wide, aes(Ageatvisit, antiminusvgslat_normantiplusvgs)) + 
#   geom_point() +
#   geom_smooth(method = "loess", se = FALSE)
# 
# Eyedata_wide$d8_f<-">200507"
# Eyedata_wide$d8_f[Eyedata_wide$d8<20050705]<-"<200507"
# 
# ggplot(Eyedata_wide, aes(d8,colour=d8_f)) + 
#   geom_histogram() +facet_wrap(~scoringtype)
# 
# min(Eyedata_wide$d8[Eyedata_wide$scoringtype=="manual"])
#########################merge all data####
Eyedata_wide_autoonly<-Eyedata_wide[Eyedata_wide$scoringtype=="auto",]
Eyedata_wide_autoonlywithcantab<-merge(Eyedata_wide_autoonly,CANTAB_Tmeasures,by=c("id","d8"),all=TRUE)
allcoglongmeasures<-merge(Eyedata_wide_autoonlywithcantab,MGS,by=c("id","d8"),all = TRUE)

allcoglongmeasures_formodels<-allcoglongmeasures[,c("id","d8","Ageatvisit","Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR",
                                                    "VGS_CRLat","SOC.Problems.solved.in.minimum.moves","SOC.Overallmeaninitialthinkingtime","SOC.Overallmeansubsequentthinkingtime",
                                                    "DMS.Percent.correct","DMS.Median.correct.latency","SSP.Span.length","nfixbreak","best_acc_m_exclude","first_lat_m_exclude")]
allcoglongmeasures_formodels_complete<-allcoglongmeasures_formodels[complete.cases(allcoglongmeasures_formodels),]
#######factor##########
allcoglongmeasures_formodels_complete_d<-allcoglongmeasures_formodels_complete[c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","VGS_CRLat","SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","DMS.Median.correct.latency","SSP.Span.length","nfixbreak","best_acc_m_exclude","first_lat_m_exclude")]
cor(allcoglongmeasures_formodels_complete_d,method="kendall")#####
#psych::fa.parallel(allcoglongmeasures_formodels_complete_d,n.iter=1000,fm="ml",fa="fa")
faggregate<-fa(allcoglongmeasures_formodels_complete_d,nfactors=3,rotate ="bifactor",fm="ml",scores="tenBerge")
######between only###############
allcoglongmeasures_formodels_complete_between<-allcoglongmeasures_formodels_complete %>% group_by(id) %>% summarise_all(funs(mean))

allcoglongmeasures_formodels_complete_between_d<-allcoglongmeasures_formodels_complete_between[,c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","VGS_CRLat","SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","DMS.Median.correct.latency","SSP.Span.length","nfixbreak","best_acc_m_exclude","first_lat_m_exclude")]
#psych::fa.parallel(allcoglongmeasures_formodels_complete_between_d,n.iter=1000,fm="ml",fa="fa")
fabetween<-fa(allcoglongmeasures_formodels_complete_between_d,nfactors=3,rotate ="bifactor",fm="ml",scores="tenBerge")
fabetweenfactanal<-factanal(allcoglongmeasures_formodels_complete_between_d,factors=3,rotation = "bifactor",scores="Bartlett")
allcoglongmeasures_formodels_complete_between$f1score<-fabetween$scores[,1]
#######
score_new_data <- function(fit, new_data, original_data) {
  means <- sapply(original_data[,row.names(fit$correlation)], mean)
  sds <- sapply(original_data[,row.names(fit$correlation)], sd)
  z <- as.matrix(scale(new_data[,row.names(fit$correlation)], 
                       center = means,
                       scale = sds))
  z %*% solve(fit$correlation, fit$loadings)
}

allscoresfrombetweenfactors<-score_new_data(fabetweenfactanal,allcoglongmeasures_formodels_complete,allcoglongmeasures_formodels_complete_between)

allcoglongmeasures_formodels_complete$f1score<-allscoresfrombetweenfactors[,1]
allcoglongmeasures_formodels_complete$f2score<-allscoresfrombetweenfactors[,2]
allcoglongmeasures_formodels_complete$f3score<-allscoresfrombetweenfactors[,3]
###############
gp<-ggplot(allcoglongmeasures_formodels_complete, aes(Ageatvisit, f1score)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)+ylab("Cognitive Control")+xlab("Age")

gpfinal<-LNCDR::lunaize(gp)
gpfinal


visitsummary<-allcoglongmeasures_formodels_complete %>% group_by(id) %>% summarise(nvisits=n())
write.csv(allcoglongmeasures_formodels_complete,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190313.csv")
write.csv(allcoglongmeasures,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03cleaneddata_20190626.csv")
##################
med_df<-allcoglongmeasures_formodels_complete
med_df$invage<-1/med_df$Ageatvisit
testmediation<-names(allcoglongmeasures_formodels_complete_between_d)
model.mf1<-lmer(as.formula("f1score~invage+(1|id)"),data=med_df)
model.mf2<-lmer(as.formula("f2score~invage+(1|id)"),data=med_df)
library(mediation)
for (var in testmediation){
  testagemodel<-lmer(as.formula(paste0(var,"~invage+(1|id)")),data=med_df)
  pvaltest<-LNCDR::lmer_extract(testagemodel,"invage")[3]
  print(var)
  print(as.numeric(pvaltest))
  if (pvaltest < .05){
    print("running f1 mediation")
    modely.yf1<-lmer(paste0(var,"~invage+f1score+(1|id)"),data=med_df)
    modely.yf2<-lmer(paste0(var,"~invage+f2score+(1|id)"),data=med_df)
    medf1 <- mediate(model.mf1, modely.yf1, treat = "invage", mediator = "f1score", boot=FALSE,sims=1000) 
  print(medf1$d0.p)
  medf2 <- mediate(model.mf2, modely.yf2, treat = "invage", mediator = "f2score", boot=FALSE,sims=1000)
  print(medf2$d0.p)
}
}


##################loadings figure##########
loadings<-as.data.frame(fabetweenfactanal$loadings[,1:3])
loadings$varname<-row.names(loadings)
loadingslong<-gather(loadings,key="factor_name","value",-varname)

loadingslong$loading_f<-as.character(round(loadingslong$value, digits = 3))
loadingslong$loading_f<-gsub("0\\.",".",loadingslong$loading_f)

tile<-ggplot(loadingslong,aes(x = factor_name,y = varname, fill = value,colour=I("black"))) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, space = "Lab",breaks=c(1,0,-1),limits=c(-.85,.85))+
  geom_text(aes(label=loading_f),size=8)

tile2<-LNCDR::lunaize(tile)+theme(axis.text.y = element_blank(),axis.text.x=element_blank())+ylab("")+xlab("")
ggsave("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/loadings.pdf",tile2,height=6,width=6)
ggsave("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/loadings.2.pdf",tile,height=6,width=6)
#########
nvisits<-allcoglongmeasures_formodels_complete %>% group_by(id) %>% summarise(nvisits=n())
allcoglongmeasures_formodels_complete<-merge(allcoglongmeasures_formodels_complete,nvisits,by="id")
#allcoglongmeasures_formodels_complete<-allcoglongmeasures_formodels_complete[allcoglongmeasures_formodels_complete$nvisits>4,]

mod<-gam(f3score~s(Ageatvisit)+s(id, bs="re"),data=allcoglongmeasures_formodels_complete)
agepred<-get_predictions(mod,cond=list(Ageatvisit=seq(min(allcoglongmeasures_formodels_complete$Ageatvisit),max(allcoglongmeasures_formodels_complete$Ageatvisit)),rm.ranef=TRUE,print.summary=FALSE))
agepred$invage<-1/agepred$Ageatvisit

allcoglongmeasures_formodels_complete$invage<-1/allcoglongmeasures_formodels_complete$Ageatvisit
allcoglongmeasures_formodels_complete$agesq<-allcoglongmeasures_formodels_complete$Ageatvisit^2

lmermod<-lmer(f3score~invage+(1|id),data=allcoglongmeasures_formodels_complete)
car::Anova(lmermod)
lmerpred<-as.data.frame(lsmeans(lmermod,~invage,at=list(invage=unique(agepred$invage))))
lmerpred$Ageatvisit<-1/lmerpred$invage

agepred$lsmean<-agepred$fit
alldatawithbothfitsgp<-ggplot(lmerpred,aes(x=Ageatvisit,y=lsmean))+geom_line(colour="red",size=2)+geom_point(data=allcoglongmeasures_formodels_complete,aes(x=Ageatvisit,y=f3score),alpha=.2)+
  geom_line(data=allcoglongmeasures_formodels_complete,aes(x=Ageatvisit,y=f1score,group=id),alpha=.2)+geom_line(data=agepred,aes(x=Ageatvisit,y=lsmean),colour="blue",size=2)

alldatawithbothfitsgp2<-LNCDR::lunaize(alldatawithbothfitsgp)+ylab("Cognitive Control\n")+xlab("\nAge")
ggsave("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/agebycognitivecontrol1.alldata.GAMplusinverse.pdf",alldatawithbothfitsgp2,height=8,width=10)


alldatawithinversegp<-ggplot(lmerpred,aes(x=Ageatvisit,y=lsmean))+geom_line(colour="red",size=2)+geom_point(data=allcoglongmeasures_formodels_complete,aes(x=Ageatvisit,y=f3score),alpha=.2)+
  geom_line(data=allcoglongmeasures_formodels_complete,aes(x=Ageatvisit,y=f1score,group=id),alpha=.2)

alldatawithinversegp2<-LNCDR::lunaize(alldatawithinversegp)+ylab("Cognitive Control\n")+xlab("\nAge")


ggsave("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/agebycognitivecontrol1.alldata.inverse.pdf",alldatawithinversegp2,height=8,width=10)

######Factor 2###############










ggplot(lmerpred,aes(x=Ageatvisit,y=lsmean))+geom_line(colour="red",size=2)+geom_point(data=allcoglongmeasures_formodels_complete,aes(x=Ageatvisit,y=f1score,alpha=.4))+geom_line(aes(group="id"))




