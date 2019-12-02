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
##############################
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
      legend.position  = "top")
}
##############################
coglongdata<-read.csv("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Data/btc_R03scoredmeasures_20190918.outlierremoved.compositeacclat.csv")
coglongdata$id<-as.factor(coglongdata$id)
#####gender file#####
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
coglongdatag$X.2<-NULL
coglongdatag$X.1<-NULL
coglongdatag$X<-NULL
####vars#########
allfactorvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","VGS_CRLat","SOC.Problems.solved.in.minimum.moves",
                 "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
                 "nfixbreak_fl","best_acc_m_exclude_fl","first_lat_m_exclude","Latencycomposite","Accuracycomposite","Accuracyfactorscores","Latencyfactorscore")

accvars<-c("Anti_CRR","Mix_CRR","DMS.Percent.correct","SSP.Span.length","nfixbreak_fl","best_acc_m_exclude_fl")

latencyvars<-c("Anti_CRLat","Mix_CRLat","VGS_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency","first_lat_m_exclude","Latencycomposite")
latencyvarsplot<-c("Latencycomposite","Anti_CRLat","first_lat_m_exclude","Mix_CRLat","VGS_CRLat","SOC.Overallmeaninitialthinkingtime","DMS.Median.correct.latency")
accvarsplot<-c("Accuracycomposite","Anti_CRR","best_acc_m_exclude_fl","Mix_CRR","nfixbreak_fl","SOC.Problems.solved.in.minimum.moves","DMS.Percent.correct","SSP.Span.length")
#######################################
gendermodelint<-outcome ~ s(pred,by=Gender) + s(id, bs = "re")
gendermodelme<-outcome ~ s(pred) + Gender + s(id, bs = "re")
coglongdatag$outcome<-coglongdatag$Accuracycomposite
coglongdatag$pred<-coglongdatag$Ageatvisit

genderinttest<-mgcv::gam(gendermodelint,data=coglongdatag,optimizer="efs")
gendermetest<-mgcv::gam(gendermodelme,data=coglongdatag,optimizer="efs")

AIC(genderinttest,gendermetest)
m<-genderinttest
interval_inc<-.1
v <- m$model[, "pred"]
cond_list <- list(seq(min(v), max(v), by=interval_inc))
names(cond_list) <- "pred"
cond_list$Gender<-c("Male","Female")
agepred <- itsadug::get_predictions(m, cond = cond_list,rm.ranef = TRUE)
#genderpred <- itsadug::get_predictions(genderinttest, cond = condlist,rm.ranef = TRUE)

ageplot <- ggplot(agepred,aes(x = pred, y = fit,colour=Gender,fill=Gender)) + 
  geom_line(size = 2) + ylab("") + 
  xlab("")+geom_ribbon(aes(ymin = fit - CI, ymax = fit + CI),alpha = 0.3)+
  geom_point(data = coglongdatag, aes(y = outcome, x = pred), alpha = 0.2)+
  geom_line(data = coglongdatag, aes(y = outcome, group = id), alpha = 0.2)+scale_color_manual(values = c("white","#F8766D","#00BFC4"))+
  scale_fill_manual(values=c("white","#F8766D","#00BFC4"))

ageplot <- ggplot(agepred,aes(x = pred, y = fit,colour=Gender,fill=Gender)) + 
  geom_line(size = 2) + ylab("") + 
  xlab("Age")+geom_ribbon(aes(ymin = fit - CI, ymax = fit + CI),alpha = 0.3)+scale_color_manual(values = c("#00BFC4","#F8766D"))+
  scale_fill_manual(values=c("#00BFC4","#F8766D"))

ageplot<-lunaize_topplot(ageplot)+theme(legend.title = element_blank())+theme(legend.position="none")
####approximate cohen's d?
forscaleing<-coglongdatag %>% group_by(Gender) %>% summarise(sdoutcome=sd(outcome))
psd<-sqrt(mean(forscaleing$sdoutcome,na.rm=TRUE))

agepredf<-agepred[agepred$Gender=="Female",c("Gender","fit","pred","CI")]
names(agepredf)<-c("Genderf","fitf","pred","CIf")

agepredm<-agepred[agepred$Gender=="Male",c("Gender","fit","pred","CI")]
names(agepredm)<-c("Genderm","fitm","pred","CIm")

agepredmf<-merge(agepredf,agepredm,by=c("pred"))
agepredmf$fitdiff<-agepredmf$fitf-agepredmf$fitm
agepredmf$estd<-agepredmf$fitdiff/psd
agepredmf$estdclip<-agepredmf$estd
agepredmf$fitfhi<-agepredmf$fitf+agepredmf$CIf
agepredmf$fitflo<-agepredmf$fitf-agepredmf$CIf
agepredmf$estdclip[agepredmf$fitm < agepredmf$fitfhi & agepredmf$fitm > agepredmf$fitflo]<-0

deriv_range<-range(agepredmf$estdclip)
tile <- ggplot(agepredmf) + aes_string(x = "pred", y = 1, fill = "estdclip") + geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "#00BFC4", mid = "white", high = "#F8766D", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)
tile_luna<-lunaize_geomraster(tile)
#############
######all plots######
scatter<-list(ageplot)
tilebot<-list(tile_luna)
allplotstp2<-c(scatter,tilebot)
topplottilegrid<-cowplot::plot_grid(plotlist=allplotstp2,ncol = 1,rel_heights=c(10,rep(1,length(tilebot))))
cowplot::save_plot("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/Genderdiff/Acccompositebygender.pdf", topplottilegrid, ncol = 1, base_height=10,base_width=12)

legend <- cowplot::get_legend(tile)
dev.off()
grid.newpage()
plotname<-sprintf("/Users/brendenclemmens/Desktop/Projects/R03_behavioral/Figures/Genderdiff/gendercohensd.pdf")
pdf(plotname,width=2,height=2)
grid.draw(legend)
dev.off()


####look more into this#######
####https://www.fromthebottomoftheheap.net/2017/10/10/difference-splines-i/

library(tidymv)

plot_difference(
  genderinttest,
  series = pred,
  difference = list(Gender = c("Female", "Male"))
)




smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
}

diffs<-smooth_diff(genderinttest,coglongdatag,"Female","Male","Gender")



library(scales)
show_col(hue_pal()(2))


