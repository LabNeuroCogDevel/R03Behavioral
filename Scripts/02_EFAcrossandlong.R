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
####functions####
factorloadingcolourplot<-function(loadingdf,varcolumn,loadingcolumns,order=NULL){
  loadingdflong<-gather(loadingdf, factor, loading, loadingcolumns, factor_key=TRUE)
  loadingdflong$varso<-loadingdflong[,varcolumn]
  if (!is.null(order)){
    loadingdflong$varso<- factor(loadingdflong$varso, levels = order)
  }
  loading_range<-range(loadingdflong$loading)
  loadinggp<-ggplot(loadingdflong,aes(x=factor,y=varso,fill=loading))+
    geom_tile(colour="black")+scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, loading_range)), limits = loading_range)
  loadinggp
  return(loadinggp)
}
lunaize_geomraster<-function(x){
  x+
    theme_bw()+
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      axis.title.x     = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.x      = element_blank(),
      legend.position  = "none")
}
lunaize_topplot<-function(x){
  x+
    theme_bw()+theme(text = element_text(size = 36))+
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      panel.border = element_blank(),
      legend.position  = "none")
}
too_small <- function(x) abs(x) < 10^-15
clip_on_sig <- function(ci){
  # if confidence interval includes zero
  # signs of x and y will be different, -x * +y  < 0
  # or if both high and low are extremly close to zero
  not_sig <- ci$ci_low * ci$ci_high < 0 |
    (too_small(ci$ci_low) & too_small(ci$ci_high))
  ci$mean_dff_clip <- ci$mean_dff
  ci$mean_dff_clip[not_sig] <- 0
  return(ci)
}
lunaize_topplot<-function(x){
  x+
    theme_bw()+theme(text = element_text(size = 36))+
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      panel.border = element_blank(),
      legend.position  = "top")
}
############
coglongdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190918.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$id)
####vars############
####################
allfactorvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","VGS_CRLat","SOC.Problems.solved.in.minimum.moves",
                 "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
                 "nfixbreak_fl","best_acc_m_exclude_fl","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.Percent.correct","SSP.Span.length","nfixbreak")
accvarsperc<-c("Anti_CRR","Mix_CRR","DMS.PC")
accvarscount<-c("SSP.Span.length","nfixbreak")

allfactorvarsorder<-c("Anti_CRR","Mix_CRR","nfixbreak_fl","best_acc_m_exclude_fl","SOC.Problems.solved.in.minimum.moves","DMS.PC","SSP.Span.length",
                      "Anti_CRLat","Mix_CRLat","VGS_CRLat","first_lat_m_exclude","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency")
##convert DMS to perc
coglongdata$DMS.PC<-coglongdata$DMS.Percent.correct/100
###add visit column###
coglongdata<-coglongdata %>% group_by(id) %>% mutate(visit=rank(d8))
coglongdata<-coglongdata[coglongdata$id!=10408,]
######EFA###############
####################
###Cross-sectional Means##########################################
coglongdata_did<-coglongdata[,c("id",allfactorvars)]
coglongdata_means<-coglongdata_did %>% group_by(id) %>% summarise_all(funs(mean))
#All continuous########
#psych::fa.parallel(coglongdata_means[,allfactorvars],n.iter=5000,fm="ml",fa="fa")
####5 factors by parallel analysis ###3 factors by eigen value of 1 rule#
fameans<-fa(coglongdata_means[,allfactorvars],nfactors=5,rotate ="bifactor",fm="ml",scores="tenBerge")
#Mixed Using factorminer#####
# categoricalvars<-c("SSP.Span.length","nfixbreak")
# coglongdata_means_mixed<-coglongdata_means %>% mutate_at(categoricalvars,ordered)
# fameansmixed<-FAMD(coglongdata_means_mixed[,allfactorvars],ncp=5)
####Cross-sectional First Visit####################################
coglongdata_v1<-coglongdata[coglongdata$visit==1,]
###scale#####
coglongdata_v1[,allfactorvars]<-lapply(coglongdata_v1[,allfactorvars], scale)
#psych::fa.parallel(coglongdata_v1[,allfactorvars],n.iter=5000,fm="ml",fa="fa")
 ##4 factors by parallel analysis 20190826###
favisit1<-fa(coglongdata_v1[,allfactorvars],nfactors=4,rotate ="geominQ",fm="ml",scores="tenBerge")
visit1temp<-coglongdata_v1[,allfactorvars]
names(visit1temp)<-abbreviate(names(visit1temp),minlength = 7,named=FALSE)
#corrplot(cor(visit1temp), method="number")
##loading vis##
loadingdf<-as.data.frame(unclass(favisit1$loadings))
loadingdf$vars<-row.names(loadingdf)
loadingcolumns<-names(loadingdf)[names(loadingdf)!="vars"]
factorloadingcolourplot(loadingdf,"vars",loadingcolumns,order=rev(allfactorvarsorder))
####Cross-sectional aggregate######################################
#psych::fa.parallel(coglongdata[,allfactorvars],n.iter=5000,fm="ml",fa="fa")
##4 factors by parallel analysis 20190826####
favaggregate<-fa(coglongdata[,allfactorvars],nfactors=4,rotate ="geominQ",fm="ml",scores="tenBerge")
######longitudinal##################################################
##################################################################
lagfun<-function(x){
  ds<-x-(dplyr::lag(x))
}
lagvars<-c(allfactorvars,"visit")
coglongdata$initialvisit<-coglongdata$visit
coglongdatadiffs<-coglongdata %>% group_by(id) %>% filter(n() > 1) %>% mutate_at(allfactorvars,lagfun)
coglongdatadiffsnona<-coglongdatadiffs[complete.cases(coglongdatadiffs),]
###longitudinal first visit###################################
#############################################################
###All subs##########
coglongdatadiffsnonav2_v1<-coglongdatadiffsnona[coglongdatadiffsnona$visit==2,]
coglongdatadiffsnonav2_v1[,allfactorvars]<-scale(coglongdatadiffsnonav2_v1[,allfactorvars])
#psych::fa.parallel(coglongdatadiffsnonav2_v1[,allfactorvars],n.iter=5000,fm="ml",fa="fa")
###PA two factors#####
fadiffv1v2<-fa(coglongdatadiffsnonav2_v1[,allfactorvars],nfactors=2,rotate ="oblimin",fm="ml",scores="tenBerge")
difftemp<-coglongdatadiffsnonav2_v1[,allfactorvars]
names(difftemp)<-abbreviate(names(difftemp),minlength = 7,named=FALSE)
##loading vis##
loadingdf<-as.data.frame(unclass(fadiffv1v2$loadings))
loadingdf$vars<-row.names(loadingdf)
loadingcolumns<-names(loadingdf)[names(loadingdf)!="vars"]
factorloadingcolourplot(loadingdf,"vars",loadingcolumns,order=rev(allfactorvarsorder))


###Under 18##########
coglongdatadiffsnonav2_v1_20<-coglongdatadiffsnona[coglongdatadiffsnona$Ageatvisit<18,]
coglongdatadiffsnonav2_v1_20z<-coglongdatadiffsnonav2_v1_20
coglongdatadiffsnonav2_v1_20z[,allfactorvars]<-scale(coglongdatadiffsnonav2_v1_20z[,allfactorvars])
#psych::fa.parallel(coglongdatadiffsnonav2_v1_20z[,allfactorvars],n.iter=5000,fm="ml",fa="fa")
fadiffv1v2U18<-fa(coglongdatadiffsnonav2_v1_20z[,allfactorvars],nfactors=4,rotate ="geomin",fm="ml",scores="tenBerge")
#####################
######################
###Multilevel SEM#####
####http://lavaan.ugent.be/tutorial/multilevel.html
coglongdata_ML<-coglongdata
names(coglongdata_ML)<-abbreviate(names(coglongdata_ML),minlength = 7,named=FALSE)
names(coglongdata_ML)<-gsub("\\.","",names(coglongdata_ML))
names(coglongdata_ML)<-gsub("_","",names(coglongdata_ML))
MLvars<-c("AntCRL","AntCRR","MxCRLt","MixCRR","VGSCRL","SOCP","SOCOvrllmnn","DMSPr","DMSM","SSPSp",
          "nfixbrk","bstc","frst")
coglongdata_ML[,MLvars]<-scale(coglongdata_ML[,MLvars])
####Two by two########
twobytwo <- '
        level: 1
          DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc
          LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM+frst
        level: 2
          DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc
          LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM+frst
'
fittwobytwo <- sem(model = twobytwo, data = coglongdata_ML, cluster = "id")
summary(fittwobytwo,standardized=TRUE)
fitMeasures(fittwobytwo)
mi<-modindices(fittwobytwo)
mi2<-mi[rev(order(mi$mi)),]
####Onebyone######
onebyone <- '
level: 1
DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst
level: 2
DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst
'
fitonebyone <- sem(model = onebyone, data = coglongdata_ML, cluster = "id")
summary(fitonebyone,standardized=TRUE)
fitMeasures(fitonebyone)
mi<-modindices(fitonebyone)
mi2<-mi[rev(order(mi$mi)),]
#####EFA Results##########
EFA <- '
        level: 1
          DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst+SOCOvrllmnn
          LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM+frst+SOCOvrllmnn
          DGA ~~ 0*LAT
        level: 2
          ACCeye =~ AntCRR+MixCRR+nfixbrk+bstc
          Lateye =~ AntCRL+MxCRLt+VGSCRL+frst
          ACCccan =~ SOCP+DMSPr+SSPSp+DMSM+SOCOvrllmnn
'
fitEFA <- sem(model = EFA, data = coglongdata_ML, cluster = "id")
summary(fitEFA,standardized=TRUE)
fitMeasures(fitEFA)
mi<-modindices(fitEFA)
mi2<-mi[rev(order(mi$mi)),]
####EFA bifactor both#######
EFA2bif <- '
        level: 1
        DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst+SOCOvrllmnn
        LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM+frst+SOCOvrllmnn
        AntCRL ~~ AntCRL
        AntCRR ~~ MixCRR
        MixCRR ~~ MxCRLt
        SOCP ~~ SOCOvrllmnn
        DGA ~~ 0*LAT
      level: 2
        DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst+SOCOvrllmnn
        LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM+frst+SOCOvrllmnn
        AntCRL ~~ AntCRL
        AntCRR ~~ MixCRR
        SOCP ~~ SOCOvrllmnn
        MixCRR ~~ MxCRLt
        DGA ~~ 0*LAT
'
EFA2bif <- '
        level: 1
DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst
LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM+frst
AntCRL ~~ AntCRL
AntCRR ~~ MixCRR
MixCRR ~~ MxCRLt
MixCRR ~~ nfixbrk
DGA ~~ 0*LAT
level: 2
DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst
LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM+frst
AntCRL ~~ AntCRL
AntCRR ~~ MixCRR
MixCRR ~~ MxCRLt
MixCRR ~~ nfixbrk
DGA ~~ 0*LAT
'

vars<-c("AntCRR","MixCRR","SOCP","DMSPr","SSPSp","nfixbrk","bstc","AntCRL","MxCRLt","VGSCRL","DMSM","frst")
coglongdata_MLdatacomplete<-coglongdata_ML[,c("id","visit","Agetvst",vars)]
coglongdata_MLdatacomplete<-coglongdata_MLdatacomplete[complete.cases(coglongdata_MLdatacomplete),]

fitEFA2bif <- sem(model = EFA2bif, data = coglongdata_MLdatacomplete, cluster = "id",missing = "ML")
summary(fitEFA2bif,standardized=TRUE)
fitMeasures(fitEFA2bif)
mi<-modindices(fitEFA2bif)
mi2<-mi[rev(order(mi$mi)),]
predlevel1<-as.data.frame(lavPredict(fitEFA2bif,level=1))
coglongdata_MLdatacomplete$DGAwithin<-predlevel1$DGA
coglongdata_MLdatacomplete$LATwithin<-predlevel1$LAT

predlevel2<-as.data.frame(lavPredict(fitEFA2bif,level=2,newdata = coglongdata_MLdatacomplete))
predlevel2$id<-unique(coglongdata_MLdatacomplete$id)
names(predlevel2)<-c("DGAbetween","LATbetween","id")

coglongdata_MLdatacompletewithinbetween<-merge(coglongdata_MLdatacomplete,predlevel2,by="id")
betweenmeans<- coglongdata_MLdatacompletewithinbetween %>% group_by(id) %>% summarise(medianage=median(Agetvst),minage=min(Agetvst),DGAbetween=median(DGAbetween),LATbetween=median(LATbetween))
write.csv(coglongdata_MLdatacompletewithinbetween,"/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/alldatawithwithinbetween.20190906.csv")

#####within figures######
##DGAwithin###
modw <- mgcv::gam(DGAwithin~s(Agetvst)+s(id,bs="re"), data=coglongdata_MLdatacompletewithinbetween,optimizer = "efs")
m<-modw
interval_inc<-.1
v <- m$model[, "Agetvst"]
cond_list <- list(seq(min(v), max(v), by=interval_inc))
names(cond_list) <- "Agetvst"
agepredDGA <- itsadug::get_predictions(m, cond = cond_list,rm.ranef = TRUE)
#genderpred <- itsadug::get_predictions(genderinttest, cond = condlist,rm.ranef = TRUE)
agepredDGA$var<-"within"
agepredDGA$rm.ranef<-NULL
agepredDGA$id<-NULL
##DGAbetween###
modb <- mgcv::gam(DGAbetween~s(medianage), data=betweenmeans,optimizer = "efs")
m<-modb
interval_inc<-.1
v <- m$model[, "medianage"]
cond_list <- list(seq(min(v), max(v), by=interval_inc))
names(cond_list) <- "medianage"
agepredDGAb <- itsadug::get_predictions(m, cond = cond_list)
#genderpred <- itsadug::get_predictions(genderinttest, cond = condlist,rm.ranef = TRUE)
agepredDGAb$var<-"between"
names(agepredDGAb)[names(agepredDGAb)=="medianage"]<-"Agetvst"
######top plot###
agepredwithinbetween<-rbind(agepredDGA,agepredDGAb)
ageplot <- ggplot(agepredwithinbetween,aes(x = Agetvst, y = fit,colour=var,fill=var)) + 
  geom_line(size = 2) + ylab("Cognitive Control")+
  xlab("Age")+geom_ribbon(aes(ymin = fit - CI, ymax = fit + CI),alpha = 0.3)+scale_colour_manual(values = c("#00BFC4","#F8766D"))+scale_fill_manual(values = c("#00BFC4","#F8766D"))

ageplot<-lunaize_topplot(ageplot)+theme(legend.title = element_blank())+theme(legend.position="none")#+theme(legend.text=element_text(size=28))
###TILES####
##within###
ciw <- LNCDR::gam_growthrate(modw, 'Agetvst', "id",n = 10000, qnt = c(0.025, 0.975))
ciw<-clip_on_sig(ciw)
ciw$se<-(ciw$ci_high-ciw$mean_dff)/1.96
ciw$sd<-ciw$se*sqrt(10000) #n interations from model 
ciw$tstat<-(ciw$mean_dff-0)/(ciw$sd/sqrt(10000))
ciw$tstat[ciw$mean_dff_clip==0]<-0
###between###
cib <- LNCDR::gam_growthrate(modb,'medianage',idvar=NULL,n = 10000, qnt = c(0.025, 0.975))
cib<-clip_on_sig(cib)
cib$se<-(cib$ci_high-cib$mean_dff)/1.96
cib$sd<-cib$se*sqrt(10000) #n interations from model 
cib$tstat<-(cib$mean_dff-0)/(cib$sd/sqrt(10000))
cib$tstat[cib$mean_dff_clip==0]<-0

deriv_range<-range(c(ciw$tstat,ciw$tstat),na.rm=TRUE)
tilewithin <- ggplot(ciw) + aes_string(x = "ages", y = 1, fill = "tstat") + geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "grey", mid = "white", high ="#F8766D" , midpoint = 0, space = "Lab", na.value="white",breaks = sort(c(0, deriv_range)), limits = deriv_range)
tile_lunawithin<-lunaize_geomraster(tilewithin)

deriv_range<-range(c(cib$tstat,cib$tstat),na.rm=TRUE)
tilebetween<- ggplot(cib) + aes_string(x = "ages", y = 1, fill = "tstat") + geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "#00BFC4", midpoint = 0, space = "Lab", na.value="white",breaks = sort(c(0, deriv_range)), limits = deriv_range)
tile_lunabetween<-lunaize_geomraster(tilebetween)

allplots<-list(tile_lunabetween,tile_lunawithin)
tilegrid<-cowplot::plot_grid(plotlist=allplots,ncol = 1)
######all plots######
scatter<-list(ageplot)
allplotstp2<-c(scatter,allplots)
topplottilegrid<-cowplot::plot_grid(plotlist=allplotstp2,ncol = 1,rel_heights=c(10,rep(1,length(allplots))))
cowplot::save_plot("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/Multivariate/SEM/DGAwithinbetween.pdf", topplottilegrid, ncol = 1, base_height=11,base_width=12)

legend <- cowplot::get_legend(tilewithin)
dev.off()
grid.newpage()
plotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/within.DGA.pdf")
pdf(plotname,width=2,height=2)
grid.draw(legend)
dev.off()

legend <- cowplot::get_legend(tilebetween)
dev.off()
grid.newpage()
plotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/between.DGA.pdf")
pdf(plotname,width=2,height=2)
grid.draw(legend)
dev.off()

####BETWEEN FIGURES####
#####within figures######
##LATwithin###
modw <- mgcv::gam(LATwithin~s(Agetvst)+s(id,bs="re"), data=coglongdata_MLdatacompletewithinbetween,optimizer = "efs")
m<-modw
interval_inc<-.1
v <- m$model[, "Agetvst"]
cond_list <- list(seq(min(v), max(v), by=interval_inc))
names(cond_list) <- "Agetvst"
agepredDGA <- itsadug::get_predictions(m, cond = cond_list,rm.ranef = TRUE)
#genderpred <- itsadug::get_predictions(genderinttest, cond = condlist,rm.ranef = TRUE)
agepredDGA$var<-"within"
agepredDGA$rm.ranef<-NULL
agepredDGA$id<-NULL
##LATbetween###
modb <- mgcv::gam(LATbetween~s(medianage), data=betweenmeans,optimizer = "efs")
m<-modb
interval_inc<-.1
v <- m$model[, "medianage"]
cond_list <- list(seq(min(v), max(v), by=interval_inc))
names(cond_list) <- "medianage"
agepredDGAb <- itsadug::get_predictions(m, cond = cond_list)
#genderpred <- itsadug::get_predictions(genderinttest, cond = condlist,rm.ranef = TRUE)
agepredDGAb$var<-"between"
names(agepredDGAb)[names(agepredDGAb)=="medianage"]<-"Agetvst"
######top plot###
agepredwithinbetween<-rbind(agepredDGA,agepredDGAb)
ageplot <- ggplot(agepredwithinbetween,aes(x = Agetvst, y = fit,colour=var,fill=var)) + 
  geom_line(size = 2) + ylab("Cognitive Control")+
  xlab("Age")+geom_ribbon(aes(ymin = fit - CI, ymax = fit + CI),alpha = 0.3)+scale_colour_manual(values = c("#00BFC4","#F8766D"))+scale_fill_manual(values = c("#00BFC4","#F8766D"))

ageplot<-lunaize_topplot(ageplot)+theme(legend.title = element_blank())+theme(legend.position="none")#+theme(legend.text=element_text(size=28))
###TILES####
##within###
ciw <- LNCDR::gam_growthrate(modw, 'Agetvst', "id",n = 10000, qnt = c(0.025, 0.975))
ciw<-clip_on_sig(ciw)
ciw$se<-(ciw$ci_high-ciw$mean_dff)/1.96
ciw$sd<-ciw$se*sqrt(10000) #n interations from model 
ciw$tstat<-(ciw$mean_dff-0)/(ciw$sd/sqrt(10000))
ciw$tstat[ciw$mean_dff_clip==0]<-0
###between###
cib <- LNCDR::gam_growthrate(modb,'medianage',idvar=NULL,n = 10000, qnt = c(0.025, 0.975))
cib<-clip_on_sig(cib)
cib$se<-(cib$ci_high-cib$mean_dff)/1.96
cib$sd<-cib$se*sqrt(10000) #n interations from model 
cib$tstat<-(cib$mean_dff-0)/(cib$sd/sqrt(10000))
cib$tstat[cib$mean_dff_clip==0]<-0

deriv_range<-range(ciw$tstat,na.rm=TRUE)
tilewithin <- ggplot(ciw) + aes_string(x = "ages", y = 1, fill = "tstat") + geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "#F8766D", mid = "white", high ="grey" , midpoint = 0, space = "Lab", na.value="white",breaks = sort(c(0, deriv_range)), limits = deriv_range)
tile_lunawithin<-lunaize_geomraster(tilewithin)

deriv_range<-range(cib$tstat)
tilebetween<- ggplot(cib) + aes_string(x = "ages", y = 1, fill = "tstat") + geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "#00BFC4", midpoint = 0, space = "Lab", na.value="white",breaks = sort(c(0, deriv_range)), limits = deriv_range)
tile_lunabetween<-lunaize_geomraster(tilebetween)

allplots<-list(tile_lunabetween,tile_lunawithin)
tilegrid<-cowplot::plot_grid(plotlist=allplots,ncol = 1)
######all plots######
scatter<-list(ageplot)
allplotstp2<-c(scatter,allplots)
topplottilegrid<-cowplot::plot_grid(plotlist=allplotstp2,ncol = 1,rel_heights=c(10,rep(1,length(allplots))))
cowplot::save_plot("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/Multivariate/SEM/LATwithinbetween.pdf", topplottilegrid, ncol = 1, base_height=11,base_width=12)
#####

legend <- cowplot::get_legend(tilewithin)
dev.off()
grid.newpage()
plotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/within.LAT.pdf")
pdf(plotname,width=2,height=2)
grid.draw(legend)
dev.off()

# legend <- cowplot::get_legend(tilebetween)
# dev.off()
# grid.newpage()
# plotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/between.DGA.pdf")
# pdf(plotname,width=2,height=2)
# grid.draw(legend)
# dev.off()






ageplot <- ggplot(agepredwithinbetween,aes(x = Agetvst, y = fit,colour=var,fill=var)) + 
  geom_line(size = 2) + ylab("Cognitive Control")+theme(legend.title=element_blank())+
  xlab("Age")+geom_ribbon(aes(ymin = fit - CI, ymax = fit + CI),alpha = 0.3)+scale_colour_manual(values = c("#00BFC4","#F8766D"))+scale_fill_manual(values = c("#00BFC4","#F8766D"))

legend <- cowplot::get_legend(ageplot)
dev.off()
grid.newpage()
plotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/Multivariate/SEM/withinbetweenleg.pdf")
pdf(plotname,width=2,height=2)
grid.draw(legend)
dev.off()


mod <- mgcv::gam(DGAwithin~s(Agetvst)+s(id,bs="re"), data=coglongdata_MLdatacompletewithinbetween,optimizer = "efs")
ci <- LNCDR::gam_growthrate(mod, 'Agetvst', "id",n = 10000, qnt = c(0.025, 0.975))
plist <- gam_growthrate_plot(coglongdata_MLdatacompletewithinbetween, mod, ci, 'Agetvst', idvar="id",xplotname='Age',draw_points = T)

mod <- mgcv::gam(LATwithin~s(Agetvst)+s(id,bs="re"), data=coglongdata_MLdatacompletewithinbetween,optimizer = "efs")
ci <- LNCDR::gam_growthrate(mod, 'Agetvst', "id",n = 10000, qnt = c(0.025, 0.975))
plist <- gam_growthrate_plot(coglongdata_MLdatacompletewithinbetween, mod, ci, 'Agetvst', idvar="id",xplotname='Age',draw_points = F)

mod <- mgcv::gam(DGAbetween~s(medianage), data=betweenmeans,optimizer = "efs")
ci <- LNCDR::gam_growthrate(mod, 'medianage',n = 10000, qnt = c(0.025, 0.975))
plist <- gam_growthrate_plot(betweenmeans, mod, ci, 'medianage', xplotname='Age',draw_maturation = F)

mod <- mgcv::gam(LATbetween~s(medianage), data=betweenmeans,optimizer = "efs")
ci <- LNCDR::gam_growthrate(mod, 'medianage',n = 10000, qnt = c(0.025, 0.975))
plist <- gam_growthrate_plot(betweenmeans, mod, ci, 'medianage', xplotname='Age',draw_maturation = F)



idage<-coglongdata_MLdatacomplete[,c("id","Agetvst")]
idage<-idage[!duplicated(idage),]





########
EFA2bifsecondlevelspec <- '
        level: 1
        DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst+SOCOvrllmnn
        LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM+frst+SOCOvrllmnn
        Eye
        DGA ~~ 0*LAT
        AntCRL ~~ AntCRL
        AntCRR ~~ MixCRR
        SOCP ~~ SOCOvrllmnn
        level: 2
        DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst+SOCOvrllmnn
        Lateye =~ AntCRL+MxCRLt+VGSCRL+frst
        ACCccan =~ SOCP+DMSPr+SSPSp+DMSM+SOCOvrllmnn
        DGA ~~ 0*Lateye
        DGA ~~ 0*ACCccan
        Lateye ~~ 0*ACCccan
        AntCRL ~~ AntCRL
        AntCRR ~~ MixCRR
        SOCP ~~ SOCOvrllmnn
'
fitEFA2bifsecond <- sem(model = EFA2bif, data = coglongdata_ML, cluster = "id")
summary(fitEFA2bifsecond,standardized=TRUE)
fitMeasures(fitEFA2bifsecond)
mi<-modindices(fitEFA2bifsecond)
mi2<-mi[rev(order(mi$mi)),]
#####################
###################





######
twobytwoBEFA <- '
        level: 1
          IN =~ AntCRR+MixCRR+nfixbrk
          LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM
          WM =~  SSPSp+bstc+frst+SOCP+DMSPr
          AntCRR ~~ MixCRR
        level: 2
          DGA =~ AntCRR+MixCRR+SOCP+DMSPr+SSPSp+nfixbrk+bstc+AntCRL+MxCRLt+VGSCRL+DMSM+frst
          LAT =~ AntCRL+MxCRLt+VGSCRL+DMSM+frst
          DGA ~~ 0*LAT
          AntCRR ~~ MixCRR
          AntCRR ~~ bstc
          MixCRR ~~ bstc
'










####transform and windsorize data#######
log10transformplusconstant<-function(x,c=(1-(min(x)))){
  #x=to-be-transformed data
  #c=constant
  xtran<-log10((x+c))
}
log10transformminusconstant<-function(x,c=(1+(max(x)))){
  #x=to-be-transformed data
  #c=constant
  xtran<-log10(c-x)
}
windzthrehold<-function(x,sdcutoff){
  x[abs(x)>sdcutoff]<-NA
}
outlierdetect<-function(df,sdcutoff=3){
  dfz<-zscorecols(df)
}
coglongdata_t<-coglongdata %>% mutate_at(accvarsperc,asin) %>% mutate_at(accvarscount,)


coglongdata_t_fd<-coglongdata_t[,allfactorvars]
coglongdata_fd<-coglongdata[,allfactorvars]






