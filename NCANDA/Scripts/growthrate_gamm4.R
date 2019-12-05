####convienence functions
ttor<-function(ts,nobs){
  rfromt<-sqrt((t^2)/(t^2+nobs))
}
find_covars_gam <- function(fml, ...) {
  ind <- as.character(fml)[3]
  # s(x1) + x2 +s(x3,"re") -> x1, x2
  vars <- unlist(strsplit(ind, "\\+"))  # formula split by +
  vars <- gsub(" ", "", vars) # remove spaces
  vars <- gsub("\\w+\\((.*?)\\)", "\\1", vars) # remove surrounding s()
  # remove random effect
  no_re <- grep("re[\"']", vars, value=T, invert=T)
  no_re <- gsub(",.*", "", no_re) # nothing after comma
  # remove anything else (likely passed in agevar)
  if (length(list(...)) != 0L) {
    no_re <- no_re[ ! no_re  %in% c(...)]
  }
  return(no_re)
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
###bootstrap and jacknife wrappers#####
jackknife_gammgrowthrate<-function(modelformula=as.formula("outcome~s(pred,k=4,fx=T)"),df,jacknifevar,idvar=NULL,mc.cores){
  mclapply(unique(df[,jacknifevar]),mc.cores=mc.cores,function(j){
    str(j)
    jdf<-df[df[,jacknifevar]!=j,]
    jm<-gamm4(modelformula,data=jdf,random=~(1|idvect))
    jf_ci<-gamm4_growthrate_maturationpoint(jm,agevar="pred",idvar="idvect")
    jf_ci$jf_id<-as.character(j)
    return(jf_ci)
  })
}



####main wrapper functions##########
gamm4_growthrate<-function(m, agevar, idvar = NULL, n.iterations = 10000, qnt = c(0.025, 
                                                                                 0.975)) 
{
  simdiff <- sim_diff1_from_gamm4(m, agevar, idvar, n.iterations = n.iterations)
  ci <- ci_from_simdiff1(simdiff$pred, simdiff$ages, qnt = qnt)
  ci$fit <- simdiff$fit
  return(ci)
}

gamm4_growthrate_maturationpoint<-function(m, agevar, idvar = NULL, n.iterations = 10000, qnt = c(0.025, 
                                                                                  0.975)) 
{
  simdiff <- sim_diff1_from_gamm4(m, agevar, idvar, n.iterations = n.iterations)
  ci <- ci_from_simdiff1(simdiff$pred, simdiff$ages, qnt = qnt)
  ci$fit <- simdiff$fit
  ci$matpoint <-gam_maturation_pointmaxsig(ci)
  ciclip<-clip_on_sig(ci)
  return(ciclip)
}

gam_maturation_pointmaxsig <- function(ci) {
  
  # when ci bounds include 0 (different sign), no longer signficant
  # clip out insignificant derivitive
  if (is.na(ci$ci_low[1])) ci <- ci[-1, ]
  
  # get mean_df_clip column
  if (! "mean_dff_clip" %in% names(ci)) ci <- clip_on_sig(ci)
  
  # find maturation point after the first signficant age
  onset_sig <- ci$ages[ci$mean_dff_clip != 0]
  maturation_pnt <- NA
  if (length(onset_sig)>0L && !all(is.na(onset_sig))) {
    #mat_points_idx <- ci$mean_dff_clip==0 & ci$ages > onset_sig[1]
    mat_points_idx<- which((ci$ages > max(ci$ages[ci$mean_dff_clip!=0]))) ###Ages greater than maximum significant (better defintion of maturation 20191130)
    if (length(mat_points_idx) > 0L && any(mat_points_idx))
      maturation_pnt <- min(ci$ages[mat_points_idx], na.rm=T) 
    if(length(onset_sig)>0L && !any(mat_points_idx))
      maturation_pnt<-max(ci$ages)
  }
  return(maturation_pnt)
}

sim_diff1_from_gamm4 <- function(m, agevar, idvar=NULL,
                               n.iterations=10000, interval_inc=.1) {
  m<-m$gam
  v <- m$model[, agevar]
  cond_list <- list(seq(min(v), max(v), by=interval_inc))
  pp <- data.frame(a=cond_list[[1]], b=Inf)
  # names should match what went into the model
  names(pp) <- c(agevar, idvar)
  
  # what if idvar is factor (Inf wont work)
  if (is.null(idvar)) {
    # do nothing. no idvar
  } else if (is.factor(m$model[, idvar])){
    # select idvar with the middle most random effect
    # random effects are coefficents like s(idvar).xxxxx
    # where xxxx is the index of the specific idvar factor name
    idvarpatt <- sprintf("s\\(%s\\)", idvar)
    idvarpatt. <- sprintf("s\\(%s\\).", idvar)
    randeff <- m$coefficients[ grep(idvarpatt, names(m$coefficients)) ]
    medval <- sort(randeff)[floor(length(randeff)/2)]
    med_re_name <- names(which(randeff == medval))
    median_idx <- gsub(idvarpatt., "", med_re_name)
    median_subj <- levels(m$model[, idvar])[as.numeric(median_idx)]
    warning("gam w/factor idvar, ",
            "setting the middle most random effect subject: ",
            median_subj)
    pp[, 2] <- median_subj
    
    # alternatively, select the first
    # pp[, 2] <- m$model[1, idvar]
  } else {
    warning("predition with continous (non-factor) idvar will give 'Inf' fit")
    # maybe pick middle value instead?
    # pp[, 2] <- mean(m$model[, idvar], na.rm=T)
  }
  
  # for all covars, pick out the mean
  for (cv in find_covars_gam(m$formula, agevar)) {
    x <- m$model[, cv]
    if (is.character(x) || is.factor(x) ){
      warning("gam w/factor covar, setting all sim to the first!")
      y <- x[1]
      # TODO: maybe pracma::Mode ?
    } else {
      y <- mean(x, na.rm=T)
    }
    pp[, cv] <- y
  }
  
  Xp <- predict(m, pp, type="lpmatrix")
  
  mu_beta <- coef(m)
  sigma_Vb <- vcov(m)
  # variance-covariance matrix of the main parameters  fitted model
  # used as: a positive-definite symmetric matrix specifying
  #  the covariance matrix of the variables.
  
  # set.seed(10)
  mrand <- MASS::mvrnorm(n.iterations, mu_beta, sigma_Vb)
  
  # ilink <- family(m)$linkinv
  # ilink<-m$family$linkinv()
  # only want inetercept and agevar
  keep_cols <- grep(paste0("Intercept|", agevar), dimnames(Xp)[[2]], value=T)
  Xp_agevar <- Xp[, keep_cols]
  mrand_agevar <- mrand[, keep_cols]
  
  # generate a whole bunch of plausable values, get the diff
  diffs <- lapply(seq_len(n.iterations), function(i)  {
    fit <- m$family$linkinv((Xp_agevar %*% mrand_agevar[i, ]))
    dff <- c(NA, diff(fit))
    return(dff)
  })
  
  return(list(pred=diffs, ages=pp[, 1], fit=predict(m, pp)))
}

ci_from_simdiff1 <- function(pred, ages, qnt=c(.025, .975)) {
  
  names(pred) <- 1:length(pred)
  mm <- t(dplyr::bind_rows(pred))
  
  # this is the ouptut !
  mean_dff <- apply(mm, 2, mean)
  ci <- apply(mm, 2, quantile, qnt, na.rm=T)
  colnames(ci) <- ages
  out <- data.frame(mean_dff=mean_dff, ages=ages)
  ci_out <- t(ci)
  dimnames(ci_out)[[2]] <- c("ci_low", "ci_high")
  return(cbind(out, ci_out))
  
  # NEVER REACHED -- left as bad documentation
  # old: return just ci and mean_dff
  return(list(ci=ci, mean_dff=mean_dff))
  
  # this is for fun
  ages[which.min(ci[1, ])]
  ages[which.min(ci[2, ])]
  
  plot(ages, mean_dff)
  for (i in 1:10) lines(ages, pred[[i]])
}


gamm_growthrate_plot <-
  function(d, model, ci, agevar, covar,idvar=NULL,
           yvar=as.character(model$formula[2]),
           plotsavename=NA, xplotname="Age", yplotname=yvar,
           draw_maturation=T, draw_points=T, show_all_fill=F,
           ci_plot=T,adjustresids=NULL){
    
    require(ggplot2)
    require(itsadug)
    modellme<-model$mer
    model<-model$gam
    # TODO:
    # remove or replace first row mean_dff
    #   NA draws weird first color on spectrum
    
    # make sure we have what we say we want
    if (! "gam" %in% class(model) ) stop("model given must be a gam model!")
    if (! "data.frame" %in% class(d) ) stop("d given must be a data.frame!")
    if (! "data.frame" %in% class(ci) ) stop("ci is not growthrate_gam() output")
    if (! yvar %in% names(model$model) ) stop(yvar, "not in model dataframe!")
    
    ci$mean_dff_clip <- ci$mean_dff
    # when ci bounds include 0 (different sign), no longer signficant
    ci <- clip_on_sig(ci)
    maturation_pnt <- gam_maturation_point(ci)
    
    # warn about no matruation point
    if (is.na(maturation_pnt) && draw_maturation) {
      warning("No maturation point!")
      draw_maturation <- F
    }
    
    # show even unsignficant change in raster if show_all_fill
    fill_column <- ifelse(show_all_fill, "mean_dff", "mean_dff_clip")
    
    ## setup derivitive raster plot
    deriv_range <- range(ci$mean_dff, na.rm=T)
    tile <-
      ggplot(ci[-1, ]) + # don't plot first row (is NA)
      aes_string(x="ages", y=1, fill=fill_column) +
      geom_raster(interpolate=TRUE) +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0,
        space = "Lab",
        breaks=sort(c(0, deriv_range)), # assumes range covers 0
        limits=deriv_range
      ) +
      xlab(sprintf("\n%s", xplotname))
    
    # draw dotted line where maturation point is
    if (draw_maturation)
      tile <- tile +
      geom_segment(
        linetype=2, colour="black",
        aes(x=maturation_pnt, xend=maturation_pnt, y=.5, yend=1.5))
    
    # lunaize the figure
    tile_luna <- lunaize_geomrasterxkeep(tile) +
      theme(text = element_text(size=36))
    
    
    # predictions
    modeldata<-data.frame(ydata=model$y, agevar=model$model[, agevar])
    condlist <- list(a=ci$ages)
    names(condlist) <- agevar
    for (cv in find_covars_gam(model$formula, agevar)) {
      x <- model$model[, cv]
      if (is.character(x) || is.factor(x) ){
        warning("gam w/factor covar, setting all sim to the first!")
        y <- x[1]
        # TODO: maybe pracma::Mode ?
      } else {
        y <- mean(x, na.rm=T)
      }
      pp[, cv] <- y
    }
    preddata<-data.frame(var=ci$ages,covar=y)
    names(preddata)<-c(agevar,covar)
    
    yhats <- predict(model,preddata,se.fit=TRUE)
    agepred<-cbind(preddata,yhats$fit,yhats$se.fit)
    names(agepred)<-c(agevar,covar,"fit","se")
    agepred$CI<-1.96*agepred$se
    ageplot<-
      ggplot(agepred) +
      aes_string(x=agevar, y="fit") +
      # solid bold line for fitted model
      geom_line(colour="black", size=2) +
      # label plot
      ylab(yplotname) +
      xlab(xplotname)
    
    if (ci_plot) {
      ageplot <- ageplot +
        geom_ribbon(aes(ymin=fit - CI, ymax=fit + CI), alpha=.3)
    }
    
    
    modeldata$resid<-adjustresids
    modeldata$id<-model$model[,idvar]
    covarname<-find_covars_gam(model$formula, agevar)
    modeldata$covar<-model$model[,covarname]
    names(modeldata)[names(modeldata)=="covar"]<-covarname
    modeldata$pred<-model$model[, agevar]
    names(modeldata)[names(modeldata=="pred")]<-agevar
    
    modelranef<-as.data.frame(ranef(modellme)$id)
    modelranef$id<-gsub("1/","",row.names(modelranef))
    names(modelranef)[names(modelranef)=="(Intercept)"]<-"randint"
    
    modeldatawithranef<-merge(modeldata,modelranef,by=c("id"))
    modeldatawithranef$yhats<- predict(model,modeldatawithranef)
    modeldatawithranef$adjustoutcome<-modeldatawithranef$yhats+modeldatawithranef$randint+modeldatawithranef$resid
    
    
    # individual points for actual data
    if (draw_points) ageplot <- ageplot +
      geom_point(data=modeldatawithranef, aes(y=adjustoutcome, x=agevar), alpha=.2)
    
    # add connecting lines if we have an idvar
    if (!is.null(idvar) && draw_points)
      ageplot <- ageplot +
      geom_line(data=modeldatawithranef, aes(y=adjustoutcome, group=id), alpha=.2)
    
    # lunaize main plot
    ageplot_luna<-LNCDR::lunaize(ageplot)+
      theme(text = element_text(size=36),
            axis.title.x=element_blank(),
            axis.text.x=element_blank())
    
    # save to file if we have plotsavename
    g <- gam_growthrate_plot_combine(ageplot_luna, tile_luna, plotsavename)
    
    list_of_plots <- list(tile=tile_luna, ageplot=ageplot_luna, both=g)
    # give back everything we created
    return(list_of_plots)
  }


#' combine age plot and tile slop heatmap into one figure (w/grob and grid)
#'
#' @description save two figures (only use if you need to mess with titles)
#' @export
#' @param ageplot_luna     ggplot plot of subject coef by age (top part of figure)
#' @param tile_luna        tile heatmap of slope  (bottom part of figure)
#' @param PDFout           PDF name to save output into, NA no saved, NULL not plotted
#' @examples
#'  data <- data.frame(age=1:100,fd_mean=1:100,subj=as.factor(letters[1:25]), conn_ahpc_vmpfc=randu[1:100,1])
#'  mod<-mgcv::gam(conn_ahpc_vmpfc~s(age)+s(fd_mean)+s(subj, bs="re"), data=data)
#'  ci<-LNCDR::gam_growthrate(mod, 'age', n = 10000, qnt = c(0.025, 0.975), idvar='subj')
#'  plist <- gam_growthrate_plot(data, mod, ci, 'age', idvar='subj')
#'  plist$tile <- plist$tile + xlab('AGE')
#'  g <- gam_growthrate_plot_combine(plist$ageplot, plist$tile, 'gammod.pdf')
gam_growthrate_plot_combine <- function(ageplot_luna, tile_luna, PDFout=NA) {
  require(grid)
  require(gridExtra)
  
  tilegrob<- ggplotGrob(tile_luna)
  agegrob <- ggplotGrob(ageplot_luna)
  
  
  g<-rbind(agegrob, tilegrob, size="first")
  panels <- g$layout$t[grep("panel", g$layout$name)]
  g$heights[panels] <- unit(c(1, .1), "null")
  
  # NULL is no draw
  # NA is draw to screen
  # filename is save to pdf
  if (is.null(PDFout)){
    return(g)
  } else if (is.na(PDFout))  {
    grid.draw(g)
  } else {
    
    # check we are saving pdf
    ext <- rev(strsplit(PDFout, "\\.")[[1]])[1]
    if (ext != "pdf") stop(PDFout, " must end in .pdf!")
    
    # draw into pdf
    pdf(PDFout, height = 9, width = 12)
    grid.draw(g)
    dev.off()
  }
  return(g)
}

lunaize_geomrasterxkeep<-function(x){
  x+
    theme_bw()+
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      legend.position  = "none")
}

# sim_absmax_from_gamm4<-function(m, agevar, idvar=NULL,
#                              n.iterations=10000, interval_inc=.1) {
#   m<-m$gam
#   v <- m$model[, agevar]
#   cond_list <- list(seq(min(v), max(v), by=interval_inc))
#   pp <- data.frame(a=cond_list[[1]], b=Inf)
#   # names should match what went into the model
#   names(pp) <- c(agevar, idvar)
#   
#   # what if idvar is factor (Inf wont work)
#   if (is.null(idvar)) {
#     # do nothing. no idvar
#   } else if (is.factor(m$model[, idvar])){
#     # select idvar with the middle most random effect
#     # random effects are coefficents like s(idvar).xxxxx
#     # where xxxx is the index of the specific idvar factor name
#     idvarpatt <- sprintf("s\\(%s\\)", idvar)
#     idvarpatt. <- sprintf("s\\(%s\\).", idvar)
#     randeff <- m$coefficients[ grep(idvarpatt, names(m$coefficients)) ]
#     medval <- sort(randeff)[floor(length(randeff)/2)]
#     med_re_name <- names(which(randeff == medval))
#     median_idx <- gsub(idvarpatt., "", med_re_name)
#     median_subj <- levels(m$model[, idvar])[as.numeric(median_idx)]
#     warning("gam w/factor idvar, ",
#             "setting the middle most random effect subject: ",
#             median_subj)
#     pp[, 2] <- median_subj
#     
#     # alternatively, select the first
#     # pp[, 2] <- m$model[1, idvar]
#   } else {
#     warning("predition with continous (non-factor) idvar will give 'Inf' fit")
#     # maybe pick middle value instead?
#     # pp[, 2] <- mean(m$model[, idvar], na.rm=T)
#   }
#   
#   # for all covars, pick out the mean
#   for (cv in find_covars_gam(m$formula, agevar)) {
#     x <- m$model[, cv]
#     if (is.character(x) || is.factor(x) ){
#       warning("gam w/factor covar, setting all sim to the first!")
#       y <- x[1]
#       # TODO: maybe pracma::Mode ?
#     } else {
#       y <- mean(x, na.rm=T)
#     }
#     pp[, cv] <- y
#   }
#   
#   Xp <- predict(m, pp, type="lpmatrix")
#   
#   mu_beta <- coef(m)
#   sigma_Vb <- vcov(m)
#   # variance-covariance matrix of the main parameters  fitted model
#   # used as: a positive-definite symmetric matrix specifying
#   #  the covariance matrix of the variables.
#   
#   # set.seed(10)
#   mrand <- MASS::mvrnorm(n.iterations, mu_beta, sigma_Vb)
#   
#   # ilink <- family(m)$linkinv
#   # ilink<-m$family$linkinv()
#   # only want inetercept and agevar
#   keep_cols <- grep(paste0("Intercept|", agevar), dimnames(Xp)[[2]], value=T)
#   Xp_agevar <- Xp[, keep_cols]
#   mrand_agevar <- mrand[, keep_cols]
#   
#   # generate a whole bunch of plausable values, get the diff
#   maxes <- lapply(seq_len(n.iterations), function(i)  {
#     fit <- m$family$linkinv((Xp_agevar %*% mrand_agevar[i, ]))
#     max <- max(abs(fit))
#     return(c(max))
#   })
#   
#   return(maxes)
# }
