HillP<-function(rasterm, w, alpha, na.tolerance,debugging){
  #
  ## Reshape values
  #
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  window = 2*w+1
  rm(hor,ver,rasterm_1,values); gc()
  #
  ## Progression bar
  #
  pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #
  ## Start the parallelized loop over iter
  #
  HillOP <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.parallel = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    HillOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
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
          message("Hill - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        p <- tw_values/sum(tw_values)
        vv <- drop(sum(p^alpha))^(1/(1-alpha))
        return(vv)
      }
    })
    return(HillOut)
  }
  return(HillOP)
}