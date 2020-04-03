ShannonP<-function(rasterm, w, na.tolerance, debugging){
  #
  ## Reshape values
  #
  values <- as.numeric( as.factor(rasterm) )
  rasterm_1 <- matrix(data = values, nrow = nrow(rasterm), ncol = ncol(rasterm))
  #
  ## Add additional columns and rows to match moving window
  #
  hor <- matrix(NA, ncol = ncol(rasterm), nrow = w)
  ver <- matrix(NA, ncol = w, nrow = nrow(rasterm)+ w * 2)
  trasterm <- cbind(ver, rbind(hor,rasterm_1,hor), ver)
  window = 2*w+1
  rm(hor, ver, rasterm_1, values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = ncol(rasterm), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  ShannonOP <- foreach(cl=(1+w):(ncol(rasterm)+w),.options.parallel = opts,.verbose = F) %dopar% {
    ShannonOut <- sapply((1+w):(nrow(rasterm)+w), function(rw) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        vv<-NA
        return(vv)
      } 
      else {
        tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }
        if( debugging ) {
          message("Shannon - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        p <- tw_values/sum(tw_values)
        vv <- (-(sum(p*log(p))))
        return(vv)
      }
    })
    return(ShannonOut)
  }
  return(ShannonOP)
}