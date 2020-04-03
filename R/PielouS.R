PielouS <- function(rasterm, w, na.tolerance, debugging){
  #message("\nStarting Pielou's index calculation:\n")
  # Reshape values
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows for moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  window = 2*w+1
  #
  ## Loop over all the pixels
  #
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        out[rw-w,cl-w]<-NA
      } else {
        tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }

        if(debugging) {
          message("\nPielou\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
        }

        tw_values <- as.vector(tw)
        maxS <- log(length(tw))
        p <- tw_values/sum(tw_values)
        p_log <- log(p)
        out[rw-w,cl-w] <- (-(sum(p*p_log)))/maxS

        if( debugging ) {
          message("\ncat: ",paste(tw,collapse=" ")," log S: ",maxS," Pielou: ",out[rw-w,cl-w])
        }

      }
    }
    svMisc::progress(value=cl/(ncol(trasterm)-1)*100, max.value=100, progress.bar = F,init=T)
  } 
  
  return(out)
}