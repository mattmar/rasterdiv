RenyiP <- function(rasterm, w, alpha, base, na.tolerance, debugging){
  # Some initial housekeeping
  window = 2*w+1
  message("\n\nProcessing alpha ",alpha, " Window ", window)
  # Set a progress bar
  pb <- txtProgressBar(title = "Iterative training", min = w, max = dim(rasterm)[2]+w, style = 3)
  #
  ## Reshape values
  #
  values <- as.numeric(as.factor(rasterm))
  rasterm_1 <- matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor <- matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver <- matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm <- cbind(ver,rbind(hor,rasterm_1,hor),ver)
  rm(hor,ver,rasterm_1,values); gc()
  #
  ## Start the parallelized loop over iter
  #
  RenyiOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
    setTxtProgressBar(pb, cl)
    if(debugging) {
      cat(paste(cl))
    }
    RenyiOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        vv <- NA
        return(vv)
      } 
      else {
        tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }
        if( debugging ) {
          message("Renyi - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        p <- tw_values/sum(tw_values)
        vv <- 1/(1-alpha) * drop(log(sum(p^alpha),base))
        return(vv)
      }
    })
    return(RenyiOut)
  } # End Renyi - parallelised
  message(("\n\n Parallel calculation of Renyi's index complete.\n"))
  return(RenyiOP)
}