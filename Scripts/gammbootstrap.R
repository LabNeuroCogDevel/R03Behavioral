#!/usr/bin/env Rscript
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
library(parallel)
library(lme4) # bootMer
library(MASS) ## for mvrnorm
library(cowplot)
lunaize_geomraster<-function(x){
  x+theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(axis.ticks.y=element_line(colour="white"))+theme(axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())+
    theme(legend.position = "none")#+theme(axis.title.x=element_blank(),axis.text.x=element_blank())
  }
#library(gbp) # bin bpacking problem solver
######################
#' @output - this thing
coglongdf<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/data/btc_R03scoredmeasures_20190313.csv")
# btc_simgam <- function() {
#  coglongdf<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/data/btc_R03scoredmeasures_20190313.csv")
#  fml <- f1score~s(Ageatvisit)+s(id, bs="re")
#  m <- gam(fml, data=coglongdf)
#  ci <- ci_from_simgam(m, 'Ageatvisit')
# }

ci_from_simgam <- function(m, agevar, n=10000) {
  simdiff <- sim_diff1_from_gam(m, agevar, n=10000)
  ci <- ci_from_simdiff1(simdiff$pred, simdiff$ages)
}

sim_diff1_from_gam <- function(m, agevar, id='id', n.iterations=10000) {

   ## From 
   # https://stats.stackexchange.com/questions/190348/can-i-use-bootstrapping-to-estimate-the-uncertainty-in-a-maximum-value-of-a-gam

   v <- m$model[, agevar]
   cond_list <- list(seq(min(v), max(v), by=.1))
   pp <- data.frame(a=cond_list[[1]], b=Inf)
   # names should match what went into the model
   names(pp) <- c(agevar, id)
    

   Xp <- predict(m, pp, type="lpmatrix")
   pp2<-pp
   pp2$visit=1
   Xp2<-predict(m, pp2, type="lpmatrix")

   mu_beta <- coef(m)
   sigma_Vb <- vcov(m) #  variance-covariance matrix of the main parameters  fitted model
   # used as: a positive-definite symmetric matrix specifying the covariance matrix of the variables.

   set.seed(10)
   mrand <- mvrnorm(n.iterations, mu_beta, sigma_Vb) 

   ilink <- family(m)$linkinv

   # assumes last column is the id/raneff 
   Xp_noid <- Xp[,1:ncol(Xp)-1]
   mrand_noid <- mrand[,1:ncol(mrand)-1]

   # generate a whole bunch of plausable values, get the diff
   pred <- lapply(seq_len(n.iterations), function(i)  {
                  pred <- ilink(Xp_noid %*% mrand_noid[i, ])
                  dff <- c(NA,diff(pred))
                  #ci <- quantile(dff, c(.025,.975)) ## get 95% CI
                 })

   return(list(pred=pred, ages=pp[,1]))
}

ci_from_simdiff1 <- function(pred, ages) {

   names(pred) <- 1:length(pred)
   mm <- t(bind_rows(pred))

   # this is the ouptut !
   mean_dff <- apply(mm, 2, mean)
   ci <- apply(mm, 2, quantile, c(.025,.975), na.rm=T)
   colnames(ci) <- ages
   return(list(ci=ci,mean_dff=mean_dff))

   # this is for fun
   ages[which.min(ci[1,])]
   ages[which.min(ci[2,])]

   plot(ages, means_dff)
   for(i in 1:10) lines(ages, pred[[i]])
}
####ploting thresholded derivs####
plotgammfactorwithderiv<-function(df,model,derivs,agevar,yvar,idvar,xplotname="Age",yplotname="fit",savename="growthrate"){
  ci<-data.frame(derivs$ci)
  derivages<-as.numeric(gsub("X","",names(ci)))
  names(ci)<-derivages
  
  cit<-as.data.frame(t(ci))
  names(cit)<-c("low","high")
  cit$age<-row.names(cit)
  
  meanderivdf<-as.data.frame(derivs$mean_dff)
  names(meanderivdf)<-"deriv"
  meanderivdf$age<-derivages
  sigages<-merge(meanderivdf,cit,by="age")
  sigages$derivthresh<-sigages$deriv
  sigages$derivthresh[sign(sigages$low)!=sign(sigages$high)]<-0
  maturationpoint<-data.frame(age=min(sigages$age[sigages$derivthresh==0],na.rm=TRUE))
  print(maturationpoint$age)
  maturationpoint$height<-1
  
  tile<-ggplot(sigages,aes(x=age,y=1,fill=derivthresh))+geom_raster(interpolate=TRUE)+scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0, space = "Lab",breaks=c(max(meanderivdf$deriv,na.rm=TRUE),0,min(meanderivdf$deriv,na.rm=TRUE)),limits=c(min(meanderivdf$deriv,na.rm=TRUE),max(meanderivdf$deriv,na.rm=TRUE)))
  tile2<-tile+geom_segment(aes(x=maturationpoint$age,xend=maturationpoint$age,y=.5,yend=1.5),linetype=2,colour="black")+xlab("\nAge")
  tile3<-lunaize_geomraster(tile2)+theme(text = element_text(size=36))
  #tile4<-tile3+theme(panel.border = element_rect(colour = "black", fill=NA, size=1),plot.margin = margin(-1, -1, -1, -1, "cm"))
  
  modeldata<-data.frame(ydata=model$y,agevar=model$model[,agevar])
  agepred<-get_predictions(model,cond=list(Ageatvisit=derivages))
  
  ageplot<-ggplot(agepred,aes_string(x=agepred[,agevar],y='fit'))+geom_line(colour="black",size=2)+geom_point(data=modeldata,aes(x=agevar,ydata),alpha=.2)+
    geom_line(data=df,aes_string(x=agevar,y=yvar,group=idvar),alpha=.2)+ylab(yplotname)+xlab(xplotname)
  ageplot2<-ageplot#+#geom_vline(xintercept = maturationpoint$age,linetype=2)
  ageplot3<-LNCDR::lunaize(ageplot2)+theme(text = element_text(size=36))+theme(axis.title.x=element_blank(),axis.text.x=element_blank())

  library(grid)
  library(gridExtra)
  plotsavename<-paste0("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/",paste0(savename,".pdf"))
  tilegrob<-ggplotGrob(tile3)
  agegrob<-ggplotGrob(ageplot3)
  
  g<-rbind(agegrob,tilegrob,size="first")
  panels <- g$layout$t[grep("panel", g$layout$name)]
  g$heights[panels] <- unit(c(1,.1),"null")
  pdf(plotsavename, height = 9, width = 12)
  grid.draw(g)
  dev.off()

}
#######
coglongdf_agerank<-coglongdf %>% group_by(id) %>% mutate(minage=min(Ageatvisit)+id/10000) %>% ungroup %>% mutate(agerank=rank(minage))
coglongdf_agerankwithvisits<-coglongdf_agerank %>% group_by(id) %>% mutate(nvisits=n())
coglongdf_agerankwithvisits$fiveormore<-"1-4 Visits"
coglongdf_agerankwithvisits$fiveormore[coglongdf_agerankwithvisits$nvisits>4]<-"5 or More Visits"

mainageplot<-ggplot(coglongdf_agerankwithvisits,aes(x=Ageatvisit,y=agerank,group=id,colour=as.factor(fiveormore)))+geom_line()+geom_point()
mainageplot<-mainageplot+scale_colour_manual(values=c("black","black"))+ylab("Subjects")+xlab("Age")
mainageplotgp<-LNCDR::lunaize(mainageplot)+theme(text = element_text(size=36))+theme(legend.position = "none")
mainageplotgp<-mainageplotgp+theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())
ggsave("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/ageplot.pdf",mainageplotgp,height=8,width=10)
#######all data
#####factor 1############
m1_factormodel<-gam(f1score~s(Ageatvisit)+s(id, bs="re"),data=coglongdf)
derivsfactor1<-ci_from_simgam(m1_factormodel, 'Ageatvisit')
derivsfactor1<-LNCDR::gam_growthrate(derivsfactor1, 'Ageatvisit')
LNCDR::gam_growthrate_plot(coglongdf,m1_factormodel,derivsfactor1,'Ageatvisit','f1score','id')
plotgammfactorwithderiv(coglongdf,m1_factormodel,derivsfactor1,'Ageatvisit','f1score','id',yplotname="Cognitive Control",savename="factor1")
####factor 2##############
m2_factormodel<-gam(f2score~s(Ageatvisit)+s(id, bs="re"),data=coglongdf)
derivsfactor2<-ci_from_simgam(m2_factormodel, 'Ageatvisit')
plotgammfactorwithderiv(coglongdf,m2_factormodel,derivsfactor2,'Ageatvisit','f2score','id',yplotname="Pure Latency",savename="factor2")
derivsfactor2<-LNCDR::gam_growthrate(m2_factormodel, 'Ageatvisit')
gam_growthrate_plot(coglongdf,m2_factormodel,derivsfactor2,'Ageatvisit','f2score','id')



###factor 1, 5 or more visits#######
coglongdfwithvisits<-coglongdf %>% group_by(id) %>% mutate(nvisits=n())
coglongdfwithvisits_fiveormore<-coglongdfwithvisits[coglongdfwithvisits$nvisits>4,]

##############
m1_factormodel<-gam(f1score~s(Ageatvisit)+s(id, bs="re"),data=coglongdfwithvisits_fiveormore)
derivsfactor1<-ci_from_simgam(m1_factormodel, 'Ageatvisit')
plotgammfactorwithderiv(coglongdfwithvisits_fiveormore,m1_factormodel,derivsfactor1,'Ageatvisit','f1score','id',yplotname="Cognitive Control",savename="factor1.fiveormorevisits")
#########
m2_factormodel<-gam(f2score~s(Ageatvisit)+s(id, bs="re"),data=coglongdfwithvisits_fiveormore)
derivsfactor2<-ci_from_simgam(m2_factormodel, 'Ageatvisit')
plotgammfactorwithderiv(coglongdfwithvisits_fiveormore,m2_factormodel,derivsfactor2,'Ageatvisit','f2score','id',yplotname="Pure Latency",savename="factor2.fiveormorevisits")
###############
##################control for practice effects######
coglongdfwithvisitno<-coglongdf %>% group_by(id) %>% mutate(visit=rank(d8))

m1_factormodel<-gam(f1score~s(Ageatvisit)+s(visit)+s(id, bs="re"),data=coglongdfwithvisitno)
derivsfactor1<-LNCDR::gam_growthrate(m1_factormodel, 'Ageatvisit')

gam_growthrate_plot(coglongdf,m1_factormodel,derivsfactor1,'Ageatvisit','f1score','id')


plotgammfactorwithderiv(coglongdfwithvisitno,m1_factormodel,derivsfactor1,'Ageatvisit','f1score','id',yplotname="Cognitive Control",savename="factor1.visitcovariation")





######
















