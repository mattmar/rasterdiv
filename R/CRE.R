CRE <- function(x, window=9, mode="classic", rescale=FALSE, na.tolerance=0.0, simplify=3, np=1, cluster.type="SOCK", debugging=FALSE)
{
#
## Define function to check if a number is an integer
#
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
#
## Initial checks
#
    if( !(is(x,"matrix") | is(x,"SpatialGridDataFrame") | is(x,"RasterLayer") | is(x,"list")) ) {
        stop("\nNot a valid x object.")
    }
    if (mode=="multidimension" & length(x)<2){
        stop("x must be a list of length > 1. Exiting...")
    }

    if( is(x,"SpatialGridDataFrame") ) {
        x <- raster(x) # Change x matrix/ces names
    }
    else if( is(x,"matrix") | is(x,"RasterLayer")) {
        rasterm<-x
    } 
    else if( is(x,"list") ) {
        rasterm<-x[[1]]
    }
    if(na.tolerance>1.0 | na.tolerance<0.0){
        stop("na.tolerance must be in the [0-1] interval. Exiting...")
    }

# Deal with matrices and RasterLayer in a different way
# If data are raster layers
    if( is(x[[1]],"RasterLayer") ) {
        if( mode=="classic" ){
            isfloat<-FALSE # If data are float numbers, transform them in integer, this may allow for a shorter computation time on big datasets.
            if( !is.wholenumber(rasterm@data@min) | !is.wholenumber(rasterm@data@max) | is.infinite(rasterm@data@min) | !is.wholenumber(median(getValues(rasterm),na.rm=T)) ){
                message("Converting x data in an integer matrix...")
                isfloat<-TRUE
                mfactor<-100^simplify
                rasterm<-getValues(rasterm)*mfactor
                rasterm<-as.integer(rasterm)
                rasterm<-matrix(rasterm,nrow(x),ncol(x),byrow=TRUE)
                gc()
            }else{
                rasterm<-matrix(getValues(rasterm),ncol=ncol(x),nrow=nrow(x),byrow=TRUE)
            }
        }
        #Print user messages
        if( mode=="classic" ){
            message("Matrix check OK: \nCumulative Residual Entropy output matrix will be returned")
        }else if( mode=="multidimension" ){
            message(("Matrix check OK: \nA matrix with multimension Cumulative Residual Entropy will be returned"))
        }else{
            stop("Matrix check failed: \nNot a valid x | method | distance, please check all these options...")
        }
# If data are a matrix or a list
    }else if( is(x,"matrix") | is(x,"list") ) {
        if( mode=="classic" ){ 
            isfloat<-FALSE # If data are float numbers, transform them in integer
            if( !is.integer(rasterm) ){
                message("Converting x data in an integer matrix...")
                isfloat<-TRUE
                mfactor<-100^simplify
                rasterm<-as.integer(rasterm*mfactor)
                rasterm<-matrix(rasterm,nrow(x),ncol(x),byrow=TRUE)
                gc()
            }else{
                rasterm<-as.matrix(rasterm)
            }
        }
        if( mode=="classic" ){
            message("Matrix check OK: \nCumulative Residual Entropy output matrix will be returned")
        }else if( mode=="multidimension" ){
            message(("Matrix check OK: \nA matrix with multimension Cumulative Residual Entropy will be returned"))
        }else{
            stop("Matrix check failed: \nNot a valid x | method | distance, please check all these options")
        }
    }

    if(np>1) {
        if(mode=="multidimension"){
            message(
                "Multi-core is not supported for multidimensional Rao, proceeding with 1 core...")
            np=1
        }else{
            message("
##################### Starting parallel calculation #######################")
        }
    }
#
## Derive operational moving window
#
    if( window%%2==1 ){
        w <- (window-1)/2
    } else {
        stop("The size of moving window must be an odd number. Exiting...")
    }
#
## Preparation of output matrices
#
    if(np==1) {
        raoqe<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
    }
#
## If mode is classic Rao
#
    if(mode=="classic") {
#
# If classic Cumulative Residual Entropy is parallelized
#
        if(np>1) {
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
            rm(hor,ver,rasterm_1,values); gc()
            if(debugging){cat("#check: Cumulative Residual Entropy parallel function.")}
#       
## Create cluster object with given number of slaves
#
            if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
                cls <- makeCluster(np,type=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
            } else if( cluster.type=="MPI" ) {
                cls <- makeCluster(np,outfile="",useXDR=FALSE,methods=FALSE,output="")
            }
            doParallel::registerDoParallel(cls)
            on.exit(stopCluster(cls)) # Close the clusters on exit
            gc()
#
## Start the parallelized loop over iter
#
            pb <- txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
            raoqe <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.parallel = opts,.verbose = F) %dopar% {
                if(debugging) {
                    cat(paste(cl))
                }
                raout <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
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
                            message("Working on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window)
                        }
                        tw_labels <- names(tw)
                        tw_values <- as.vector(tw)
                        #if clause to exclude windows with only 1 category
                        if( length(tw_values) <2 ) {
                            vv<-NA
                            return(vv)
                        }
                        else {
                            vv <- .CRE_(tw_values)
                            return(vv)
                        }
                    }
                })
                return(raout)
            } # End classic RaoQ - parallelized
            message(("\n\nCalculation of Cumulative Residual Entropy complete.\n"))
#
## If classic RaoQ is sequential
#
        } else if(np==1) {
# Reshape values
            values<-as.numeric(as.factor(rasterm))
            rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
# Add additional columns and rows for moving window
            hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
            ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
            trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
# Derive distance matrix
            classes<-levels(as.factor(rasterm))
# Loop over each pixel
            for (cl in (1+w):(dim(rasterm)[2]+w)) {
                for(rw in (1+w):(dim(rasterm)[1]+w)) {
                    if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
                        raoqe[rw-w,cl-w]<-NA
                    } else {
                        tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
                        if( "NA's"%in%names(tw) ) {
                            tw<-tw[-length(tw)]
                        }
                        if(debugging) {
                            message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window)
                        }
                        tw_labels <- names(tw)
                        tw_values <- as.vector(tw)
                        #if clause to exclude windows with only 1 category
                        if(length(tw_values) < 2) {
                            raoqe[rw-w,cl-w]<-NA
                        } else {
                            if(isfloat) {
                                raoqe[rw-w,cl-w]<-.CRE_(tw_values/mfactor)
                            } else {
                                raoqe[rw-w,cl-w]<-.CRE_(tw_values)
                            }
                        }
                    } 
                    svMisc::progress(value=cl/(ncol(trasterm)-1)*100, max.value=100, progress.bar = F,init=T)
                } 
            } # End of for loop 
            message(("\nCalculation of Cumulative Residual Entropy complete.\n"))
        }
    }  # End classic RaoQ - sequential
    else if( mode=="multidimension" ){
        if(debugging) {
            message("#check: Into multidimensional clause.")
        }
#----------------------------------------------------#
#
## If multimensional Cumulative Residual Entropy
#
# Check if there are NAs in the matrices
        if ( is(rasterm,"RasterLayer") ){
            if(any(sapply(lapply(unlist(x),length),is.na)==TRUE))
                message("\n Warning: One or more RasterLayers contain NA's which will be treated as 0")
        } else if ( is(rasterm,"matrix") ){
            if(any(sapply(x, is.na)==TRUE) ) {
                message("\n Warning: One or more matrices contain NA's which will be treated as 0")
            }
        }
#
## Reshape values
#
        vls<-lapply(x, function(x) {raster::as.matrix(x)})
#
## Rescale and add additional columns and rows for moving w
#
        hor<-matrix(NA,ncol=dim(vls[[1]])[2],nrow=w)
        ver<-matrix(NA,ncol=w,nrow=dim(vls[[1]])[1]+w*2)
        if(rescale) {
            trastersm<-lapply(vls, function(x) {
                t1 <- raster::scale(raster(cbind(ver,rbind(hor,x,hor),ver)))
                t2 <- raster::as.matrix(t1)
                return(t2)
            })
        } else {
            trastersm<-lapply(vls, function(x) {
                cbind(ver,rbind(hor,x,hor),ver)
            })
        }
        if(debugging) {
            message("#check: After rescaling in multimensional clause.")
        }
#
## Loop over all the pixels in the matrices
#
        if( (ncol(vls[[1]])*nrow(vls[[1]]))> 10000) {
            message("\n Warning: ",ncol(vls[[1]])*nrow(vls[[1]])*length(vls), " cells to be processed, it may take some time... \n")
        }
        for (cl in (1+w):(dim(vls[[1]])[2]+w)) {
            for(rw in (1+w):(dim(vls[[1]])[1]+w)) {
                if( length(!which(!trastersm[[1]][c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
                    raoqe[rw-w,cl-w] <- NA
                } else {
                    tw<-lapply(trastersm, function(x) { x[(rw-w):(rw+w),(cl-w):(cl+w)]
                })
#
##Vectorize the matrices in the list and calculate
#Among matrix pairwase distances
                    lv <- lapply(tw, function(x) {as.vector(t(x))})
                    raoqe[rw-w,cl-w] <- .CRE_(data.frame(lv))
                }
            }
            progress(value=cl, max.value=dim(rasterm)[2]+w, progress.bar = FALSE)
        }
        if(exists("pb")) {
         close(pb) 
         message("\nCalculation of Multidimensional Cumulative Residual Entropy index complete.\n")
     }
 } else{
    message("Something went wrong when trying to calculate Rao's index.")
}  # end of multimensional RaoQ

#----------------------------------------------------#

#
## Return multiple outputs
#
if(debugging){
    message( "#check: return function." )
}
if( mode=="classic" ) {
   if( isfloat & np>1 ) {
    return(do.call(cbind,raoqe)/mfactor)
    if(debugging){
        message("#check: return function - classic.")
    }
} else if( !isfloat & np>1 ) {
    return(do.call(cbind,raoqe))
} else { return(raoqe) }
} else if( mode=="multidimension" ) {
    outl <- list(raoqe)
    names(outl) <- c("Multidimension_CRE")
    return(outl)
}
}

#----------------------------------------------------#

#
##Supporting function for CRE
#
.CRE_<-function(B,base=exp(1))
{
  #Cumulative Residual Entropy
  P=.Prob(B)
  Pcre=.CumRes(P)
  -sum(Pcre*log(Pcre,base)*.Deltas(P), na.rm=TRUE)
}

.Prob<-function(C)
{
  #Point probability
  if( is.null(dim(C))){L=length(C)
} else {L=dim(C)[1]}
table(C)/L
}

.Deltas<-function(P, first=0)
{
  #Difference among values of a table.
  #For multidimensional table the product is given 
  if ((length(dim(P)))==1){
      delta=c(first,diff(as.numeric(names(P))))} else {
          delta=1
          for (dim in dimnames(P)){
            delta=outer(delta,c(first,diff(as.numeric(dim))))
        }
    }
    drop(delta)
}

.CumRes<-function(a)
{
  #Calculate Cumulative Residual Probability
  D=dim(a)
  if (length(dim(a))==1){
    return( rev( cumsum(rev(a)) )) }
    atemp=a
    for(i in 1:length(D))
    {
        aa=.Rev(.Cumsum(.Rev(atemp,i),i),i)
        atemp=aa
    }
    atemp
}

.Reorder<-function(a,ax)
{
  #Reordering dimension required after the use of apply over more than one dimension  
  D=dim(a)
  maxdim=length(D)
  reorder=2:maxdim
  end=reorder[(ax):(maxdim-1)]
  if (maxdim==ax){end=reorder[0]}
  reorder=c(reorder[0:(ax-1)],1,end)
  aperm(a,reorder)
  
}
.Cumsum<-function(a,ax=1)
{
  #Cumulative sum along a specific dimension
  D=dim(a)
  dimen=1:length(D)
  .Reorder(apply(a, dimen[-ax],cumsum),ax)
}
.Rev<-function(a,ax)
{
  #Reverse order along a specific dimension
  D=dim(a)
  dimen=1:length(D)
  .Reorder(apply(a, dimen[-ax], rev),ax)
}