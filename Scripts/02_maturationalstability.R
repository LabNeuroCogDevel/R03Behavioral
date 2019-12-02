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
######functions#####
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
multi_gam_growthrate_plotseperate<-function(df,outcomevars,predvars='Ageatvisit',model,plotspecifier,idvar){
  ###wrapper script for gam_growthrate######
  pairs<-as.data.frame(expand.grid(outcomevars,predvars))
  names(pairs)<-c("outcome","pred")
  for (p in 1:nrow(pairs)){
  df$outcome<-unlist(df[,as.character(pairs$outcome[p])])
  df$pred<-unlist(df[,as.character(pairs$pred[p])])
  m<-mgcv::gam(model,data=df)
  ggr<-gam_growthrate(m,agevar="pred",idvar=idvar)
  plotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/%s.%s.pdf",plotspecifier,as.character(pairs$outcome[p]))
  gam_growthrate_plot(df,m,ggr,agevar='pred',yvar="outcome",idvar=idvar,xplotname="age",yplotname=as.character(pairs$outcome[p]),plotsavename = plotname)
  ggr$var<-as.character(pairs$outcome[p])
  }
}
multi_gam_growthrate<-function(df,outcomevars,predvars='Ageatvisit',model,idvar){
  ###wrapper script for gam_growthrate######
  pairs<-as.data.frame(expand.grid(outcomevars,predvars))
  names(pairs)<-c("outcome","pred")
  ggrout<-NULL
  for (p in 1:nrow(pairs)){
    print(as.character(pairs$outcome[p]))
    df$outcome<-unlist(df[,as.character(pairs$outcome[p])])
    df$pred<-unlist(df[,as.character(pairs$pred[p])])
    m<-mgcv::gam(model,data=df)
    ggr<-gam_growthrate(m,agevar="pred",idvar=idvar)
    ggr$var<-as.character(pairs$outcome[p])
    ggrout<-rbind(ggrout,ggr)
  }
  return(ggrout)
}
multi_gam_growthrate_multiplot<-function(df,outcomevars,predvars='Ageatvisit',model,plotspecifier,idvar,topscatter=NULL,topmodel=NULL,totaltiles=NULL){
  ###wrapper script for gam_growthrate######
  pairs<-as.data.frame(expand.grid(outcomevars,predvars))
  names(pairs)<-c("outcome","pred")
  ggrout<-NULL
  for (p in 1:nrow(pairs)){
    df$outcome<-unlist(df[,as.character(pairs$outcome[p])])
    df$pred<-unlist(df[,as.character(pairs$pred[p])])
    m<-mgcv::gam(model,data=df,optimizer="efs")
    ggr<-gam_growthrate(m,agevar="pred",idvar=idvar)
    ggr$var<-as.character(pairs$outcome[p])
    ggrout<-rbind(ggrout,ggr)
  }
  ggrouttoplot<-ggrout[!is.na(ggrout$mean_dff),]
  ggrouttoplotclip<-clip_on_sig(ggrouttoplot)
  ggrouttoplotclip$se<-(ggrouttoplotclip$ci_high-ggrouttoplotclip$mean_dff)/1.96
  ggrouttoplotclip$sd<-ggrouttoplotclip$se*sqrt(10000) #n interations from model 
  ggrouttoplotclip$tstat<-(ggrouttoplotclip$mean_dff-0)/(ggrouttoplotclip$sd/sqrt(10000))
  
  ggrouttoplotclip$percentdiff<-ggrouttoplotclip$mean_dff/ggrouttoplotclip$fit
  ggrouttoplotclip$tstat[ggrouttoplotclip$mean_dff_clip==0]<-0
  
  ggrouttoplotclip %>% group_by(var) %>% summarise(r=range(tstat)[1],r2=range(tstat)[2])
  deriv_range<-range(ggrouttoplotclip$tstat)
  
  allplots<-list()
  for (vi in 1:length(unique(ggrouttoplotclip$var))){
    print(vi)
    v<-unique(ggrouttoplotclip$var)[vi]
    ci<-ggrouttoplotclip[ggrouttoplotclip$var==v,]
    tile <- ggplot(ci) + aes_string(x = "ages", y = 1, fill = "tstat") + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)
    #+ggtitle(v)
    if (vi !=length(unique(ggrouttoplotclip$var))){
    tile_luna <- lunaize_geomraster(tile) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())#,axis.text.x=element_text(color="white"),axis.ticks.x=element_line(color="white"))
    }else{
      tile_luna <- lunaize_geomraster(tile) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    }
    allplots[[vi]]<-tile_luna
  }
  if (!is.null(totaltiles)){
    print(totaltiles)
    addtiles<-totaltiles-length(allplots)
    cifill<-ci
    cifill$tstat<-0
    tilefill <- ggplot(cifill) + aes_string(x = "ages", y = 1, fill = "tstat") + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)
    tilefillluna <- lunaize_geomraster(tilefill) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    for (ai in seq(1,addtiles)){
      print(ai)
      ati<-length(allplots)+ai
      allplots[[ati]]<-tilefillluna
    }
  }
  legend <- cowplot::get_legend(tile)
  plotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/%s.pdf",plotspecifier)
  tilegrid<-cowplot::plot_grid(plotlist=allplots,ncol = 1)
  save_plot(plotname, tilegrid, ncol = 1, base_height=8,base_asp = 1.1)
  dev.off()
  grid.newpage()
  plotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/%s.legend.pdf",plotspecifier)
  pdf(plotname,width=2,height=2)
  grid.draw(legend)
  dev.off()
 
  if (!is.null(topscatter) && !is.null(topmodel)){
    print("adding top scatter based on top model and predvars[1]")
    topmodelm<-gam(topmodel,data=df,optimizer="efs")
    ci<-gam_growthrate(topmodelm,predvars[1],idvar="id")
    scatter<-LNCDR::gam_growthrate_plot(df,topmodelm,ci,agevar=predvars[1],idvar=idvar)
    scatter<-list(lunaize_topplot(scatter$ageplot))
    allplotstp<-allplots
    allplotstp2<-c(scatter,allplotstp)
    topplottilegrid<-cowplot::plot_grid(plotlist=allplotstp2,ncol = 1,rel_heights=c(10,rep(1,length(allplots))))
    topplotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/%s.withtopplot.pdf",plotspecifier)
    save_plot(topplotname, topplottilegrid, ncol = 1, base_height=12,base_width =10)
    }
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
latencyvarsplot<-c("Latencycomposite","Anti_CRLat","first_lat_m_exclude","Mix_CRLat","Mix_CRLat","VGS_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency")
accvarsplot<-c("Accuracycomposite","Anti_CRR","best_acc_m_exclude_fl","Mix_CRR","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","SSP.Span.length")

####################plot seperate########
model<-outcome ~ s(pred) + s(id, bs = "re")
#multi_gam_growthrate_plotseperate(coglongdata,allfactorvars,predvars='Ageatvisit',model,plotspecifier="/Univariate/NoVisit/no.visit",idvar="id")

coglongdata$visitf<-as.factor(coglongdata$visit)
model2<-outcome ~ s(pred)+s(visit)+s(id, bs = "re")
#multi_gam_growthrate_plotseperate(coglongdata,allfactorvars,predvars='Ageatvisit',model2,plotspecifier="/Univariate/VisitCov/with.visit",idvar="id")
###########plot together##################
# multi_gam_growthrate_multiplot(coglongdata,latencyvars,predvars='Ageatvisit',model,plotspecifier="/Univariate/NoVisit/grid.novisit",idvar="id")
# multi_gam_growthrate_multiplot(coglongdata,latencyvars,predvars='Ageatvisit',model2,plotspecifier="/Univariate/VisitCov/grid.visit",idvar="id")

##Latency
topmodel<-Latencycomposite~s(Ageatvisit)+s(id,bs="re")
multi_gam_growthrate_multiplot(coglongdata,latencyvarsplot,predvars='Ageatvisit',model2,plotspecifier="/Univariate/VisitCov/grid.visit",idvar="id",topscatter=TRUE,topmodel=topmodel,totaltiles = length(accvarsplot))
##########
#######
###accuracy
topmodel<-Accuracycomposite~s(Ageatvisit)+s(visit)+s(id,bs="re")
multi_gam_growthrate_multiplot(coglongdata,accvarsplot,predvars='Ageatvisit',model2,plotspecifier="/Univariate/VisitCov/acc.grid.visit",idvar="id",topscatter=TRUE,topmodel=topmodel)










###visit1########
mggacc<-multi_gam_growthrate(coglongdata,accvars,predvars='Ageatvisit',model,idvar="id")
coglongwindsv1<-coglongwinds[coglongwinds$visit==1,]
model<-outcome ~ s(pred)
multi_gam_growthrate(coglongwindsv1,allfactorvars,predvars='Ageatvisit',model,plotspecifier="no.visit.visit1",idvar=NULL)








###  f<-paste0(as.character(pairs$outcome[p]),sprintf("~s(%s)",pairs$pred[p]))
# if (!is.null(covars)){
# }
# if (!is.null(idvar)){
#   f<-paste0(f,sprintf("+s(%s,bs=re)",idvar))
# }



