RenyiS <- function(rasterm, w, alpha, base, na.tolerance, debugging){

  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Reshape values
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Add additional columns and rows for moving window
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  window = 2*w+1
  # Loop over each pixel
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
     if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
      out[rw-w,cl-w]<-NA
    } else {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw<-tw[-length(tw)]
      }
      if(debugging) {
        message("Renyi\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- tw_values/sum(tw_values)
      out[rw-w,cl-w]<-1/(1-alpha) * drop(log(sum(p^alpha),base))
    }
  }
  svMisc::progress(value=cl/(ncol(trasterm)-1)*100, max.value=100, progress.bar = F,init=T)
}
return(out)
}