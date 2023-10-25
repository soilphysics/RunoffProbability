require(svMisc)
require(digest)

ThresholdGLS <- function(m.th,df.th,log.trans=FALSE){
  m.var <- attr(m.th$terms,"term.labels")[2]
  #factor.levels<-as.character(levels(df.sub$Site))
  factor.levels<-rownames(m.th$contrasts$Site)
  results<-data.frame(threshold=rep(0,length(factor.levels)),row.names = factor.levels)
  results$Treatment<-rownames(results)
  results$nbth <- 0
  results$nath <- 0
  for(i in 1:length(factor.levels)){
    ## first we determine if the factor level is the ==============================
    ## first factor level in the list. In that case, we only need to consider 
    ## the (intercept) and the precipitation to calculate the threshold
    if(i==1){
      th <- as.numeric(-m.th$coefficients[1]/m.th$coefficients[2])
    }else{
      ### extract the variable-factor level combination ------------------------------
      m.vf <- paste0(m.var,factor.levels[i])
      m.vf2 <- paste0(attr(m.th$terms,"term.labels")[1],":",m.vf)
      ### determine the intercept ------------------------------
      m.int <-m.th$coefficients[1]+m.th$coefficients[names(m.th$coefficients)==m.vf]
      ### determine the slope ------------------------------
      m.slope<- m.th$coefficients[2]+m.th$coefficients[names(m.th$coefficients)==m.vf2]
      th <- -m.int/m.slope
    }
    results[rownames(results)==factor.levels[i],]$threshold <- th
    # add below/above th points
    if(log.trans==TRUE){
      th.2<-10^th
    }else{
      th.2<-th
    }
    df.sub.th <- subset(df.th,Site==factor.levels[i])
    # calculate number of obs. that have runoff and are below th
    nbth <- sum(df.sub.th$Precipitation<th.2 & df.sub.th$Runoff>0)
    # calculate number of obs. that have no runoff and are above th
    nath <- sum(df.sub.th$Precipitation>th.2 & df.sub.th$Runoff==0)
    results[rownames(results)==factor.levels[i],]$nbth <- nbth
    results[rownames(results)==factor.levels[i],]$nath <- nath
    }
  return(results)
}


BootstrappingGLS <- function(m.bs, n.bootstrap, df.th,log.trans=FALSE) {
  # create an ID for the data set
  md5.string <-
    paste0(df.th$Source[1], n.bootstrap, log.trans,m.bs$call, collapse = '')
  # create a file name and check if it already exists
  fn <- paste0("bootstrapping/", digest(md5.string), ".rda")
  if (!file.exists(fn)) {
    # set seed to obtain reproducible results
    set.seed(123)
    # extract the number of variables
    n.v <- length(coef(m.bs))
    # create empty coefficient matrix
    coefmat <- matrix (NA, n.bootstrap, n.v)
    # extract residuals
    resids = residuals(m.bs, type = "normalized")
    # extract fitted values
    preds <- fitted(m.bs)
    # perform resampling
    for (i in 1:n.bootstrap) {
      ## create new y values by adding randomly adding
      ## residuals to the fitted model values
      
      bs.y <- preds + sample(resids, replace = T)
      ## add bs.y to df.sub so the update function 'finds' it
      df.sub$bs.y <- bs.y
      ## Then perform a new model fit
      try({
        m.bs.hat <- update(m.bs, bs.y ~ .)
        ## write coeffients into the coeffient matrix
        coefmat[i, ] <- coef(m.bs.hat)
      })
    }
    cat("\n")
    # add column names to coefficient matrix
    coefmat <- data.frame(coefmat)
    names(coefmat) <- names(coef(m.bs))
    # extract factor levels for the variable Sites
    factor.levels <- rownames(m.bs$contrasts$Site)
    # create output data frame
    results <- data.frame(
      CI_low = rep(0, length(factor.levels)),
      CI_high = rep(0, length(factor.levels)),
      row.names = factor.levels
    )
    # calculate 2.5%-97.5% confidence intervalls
    for (i in 1:length(factor.levels)) {
      ## first we determine if the factor level is the ==============================
      ## first factor level in the list. In that case, we only need to consider
      ## the (intercept) and the precipitation to calculate the threshold
      m.var <- attr(m.bs$terms, "term.labels")[2]
      if (i == 1) {
        th <- as.numeric(-coefmat[, 1] / coefmat[, 2])
      } else{
        ### extract the variable-factor level combination ------------------------------
        m.vf <- paste0(m.var, factor.levels[i])
        m.vf2 <-
          paste0(attr(m.bs$terms, "term.labels")[1], ":", m.vf)
        ### determine the intercept ------------------------------
        m.int <-
          coefmat[, 1] + coefmat[, names(m.bs$coefficients) == m.vf]
        ### determine the slope ------------------------------
        m.slope <-
          coefmat[, 2] + coefmat[, names(m.bs$coefficients) == m.vf2]
        th <- as.numeric(-m.int / m.slope)
      }
      CI <- quantile(th, c(0.025, 0.975), na.rm = T)
      results[rownames(results) == factor.levels[i], ]$CI_low <-
        CI[1]
      results[rownames(results) == factor.levels[i], ]$CI_high <-
        CI[2]
    }
    results$Treatment <- rownames(results)
    r.th <- ThresholdGLS(m.bs, df.th,log.trans)
    results.m <- merge(results, r.th)
    # write file
    save(results.m, file=fn, compress = "xz")
  } else{
    load(fn)
  }
  return(results.m)
}

FormatBS <- function(threshold,CI_low, CI_high){
  test <- paste0(format(round(threshold,1),digits=1,nsmall=1)," [",
                 format(round(CI_low,1),digits=1,nsmall=1),";",
                 format(round(CI_high,1),digits=1,nsmall=1),"]")
  return(test)
}
PerformBS <- function(gls.model){
  # Extract residuals and fitted values
  bp.df<- data.frame(resids = residuals(gls.model,type = "normalized"),
                     fitted = fitted(gls.model))
  
  # Perform Breusch-Pagan test
  test <- bptest(resids ~ fitted,data=bp.df)
  return(test$p.value)
}

CreateResultsTable <- function(author.name,fl){
  df<-data.frame(Article=rep(author.name,length(fl)),
             Treatment = fl,
             lm =rep(0,length(fl)),
             lm_nbth =rep(0,length(fl)),
             lm_nath =rep(0,length(fl)),
             gls =rep(0,length(fl)),
             gls_nbth =rep(0,length(fl)),
             gls_nath =rep(0,length(fl)),
             log10.lm =rep(0,length(fl)),
             log10.lm_nbth =rep(0,length(fl)),
             log10.lm_nath =rep(0,length(fl)),
             log10.gls =rep(0,length(fl)),
             log10.gls_nbth =rep(0,length(fl)),
             log10.gls_nath =rep(0,length(fl)))
  return(df)
}
