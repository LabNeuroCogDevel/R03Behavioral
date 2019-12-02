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
########################
######################
coglongdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190918.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$id)
coglongdata$age<-coglongdata$Ageatvisit
gender<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/gender_20190825.csv")
gender$X<-NULL
#10132=MALE (Single Visit it was logged incorrectly)
#10292=FEMALE (Single Visit it was logged incorrectly)
gender$Gender[gender$id==10132]<-"Male"
gender$Gender[gender$id==10292]<-"Female"
nodupgender<-gender[!duplicated(gender),]
nodupgenderm<-nodupgender[nodupgender$id %in% unique(coglongdata$id),]
dupes<-nodupgenderm[duplicated(nodupgenderm$id) | duplicated(nodupgenderm$id, fromLast = TRUE),]
dupes_allvisits<-gender[gender$id %in% dupes$id,]

coglongdatag<-merge(coglongdata,nodupgenderm,by="id",all.x=TRUE,all.y=FALSE)
#####################################
age_ranked <- waterfall_group(coglongdatag)

p <- ggplot(age_ranked) + aes(x = age, y = age_id, group = age_id) + geom_line() + geom_point()
p <- lunaize(p) + ylab("Subjects\n") + xlab("\n Age") + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())




p <- ggplot(age_ranked) + aes(x = age, y = age_id, group = age_id,colour=Gender) + geom_line() + geom_point()
p <- lunaize(p) + ylab("") + xlab("Age") + theme(axis.title.y = element_blank(), 
                                                 axis.ticks.y = element_blank(), axis.text.y = element_blank())+
  scale_colour_manual(values=c("white","#F8766D","#00BFC4"))+theme(legend.title = element_blank(),legend.position="top") 

###hey


