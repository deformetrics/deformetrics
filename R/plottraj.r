#'plottraj
#'
#' This function plot trajectory of data in two or three dimensions in the order they appear in the input matrix
#' @param matrix a matrix kx2 or kx3
#' @param factor character: variable factor that affiliates each row to a trajectory
#' @param label character vector: optional vector of names
#' @param sizelabel numeric: size of the label text 
#' @param col numeric: col argument in par()
#' @param alpha=NULL transparency value 
#' @param pch numeric: pch argument in par()
#' @param lty numeric: lty argument in par()
#' @param cex numeric: cex argument in par()
#' @param axlab=NULL
#' @param xlim numeric vector: xlim argument in par() (default=range(matrix[,1]))
#' @param ylim numeric vector: ylim argument in par() (default=range(matrix[,2]))
#' @param asp numeric: giving the aspect ratio y/x (default=abs((xlim[1]-xlim[2])/(ylim[1]-ylim[2])))
#' @param reorder logical: if TRUE group is reordered according the the order of levels appearance
#' @param title character: optional name for the title (default=NULL)
#' @param add logical: if TRUE the plot is added to existing plot
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' data(case1)
#' group<-factor(rep(1:5,each=21))
#' plottraj(procSym(case1,CSinit=T)$PCscores[,1:2],group,asp=1)######## a classic approach does not recover the cycle
#' dtcase1<-newdt(case1,group,CR=procSym(case1[,,firstsfac(group)],CSinit=F,reflect=F,scale=F,pcAlign=F)$mshape,locs=case1[,,firstsfac(group)],tolrot=10)
#' plottraj(procSym(dtcase1,CSinit=T)$PCscores[,1:2],group,asp=1) ####  after the parallel transport it is correctly recovered
#' }
#' @export
plottraj<-function(matrix,factor,label=NULL,sizelabel=0,col=NULL,alpha=NULL,pch=NULL,lty=NULL,cex=NULL,axlab=NULL,xlim=range(matrix[,1]),ylim=range(matrix[,2]),asp=abs((xlim[1]-xlim[2])/(ylim[1]-ylim[2])),reorder=T,title=NULL,add=F){
  require(calibrate)
  require(compositions)
  positions<-NULL
  for(k in 1:max(as.numeric(factor))){
    factor<-factor(factor,levels=unique(factor))
    positioni<-which(as.numeric(factor)==k)
    positions<-c(positions,list(positioni))
  }
  names(positions)<-levels(factor)
  sequ<-NULL
  for(i in 1:max(length(positions))){
    seqi<-c(1:length(positions[[i]]))
    sequ<-c(sequ,list(seqi))
  }
  
  pure<-rep(NA,length(unlist(positions)))
  for(j in 1:max(length(positions))){
    pure<-replace(pure,positions[[j]],c(sequ[[j]]))
  }
  
  
  if(!is.null(label)){
    label<-label}else{
      label<-paste(pure,as.character(factor),sep="_")
    }
  
  if(is.null(col)){col<-as.numeric(factor)}else{col<-col}
  
  if(is.null(pch)){pch<-rep(19,length(factor))}else{pch<-pch}
  
  if(is.null(alpha)){alpha<-rep(1,length(factor))}else{alpha<-alpha}
  
  if(is.null(cex)){cex<-rep(0.7,length(factor))}else{cex<-cex}
  
  if(is.null(lty)){lty<-rep(1,nlevels(factor))}else{lty<-lty}
  
  
  
  if(!is.null(axlab)){
    xlabel<-axlab[1]
    ylabel<-axlab[2]
  }else{
    xlabel<-c("X")
    ylabel<-c("Y")
  }
  
  
  if(!is.null(axlab)&dim(matrix)[2]>2){
    vlabs=c(axlab[1],axlab[2],axlab[3])
  }else{vlabs = c("x", "y", "z")}
  
  #if(asp==F){asp=1}else{asp=asp}
  
  if(dim(matrix)[2]<3){
    
    for (i in 1 :nlevels(factor)){
      plot(matrix[,1][as.numeric(factor)==i],matrix[,2][as.numeric(factor)==i],col=col[as.numeric(factor)==i],xlim=xlim,ylim=ylim,xlab=axlab[1],ylab=axlab[2],pch=pch[as.numeric(factor)==i],cex=cex[as.numeric(factor)==i],asp=asp,main=title)
      lines(matrix[,1][as.numeric(factor)==i],matrix[,2][as.numeric(factor)==i],col=col[as.numeric(factor)==i],lty=lty[i])
      par(new=T)
    }
    if(!is.null(label)&sizelabel>0){
      for(i in 1:nlevels(factor)){
        textxy(matrix[,1][as.numeric(factor)==i],matrix[,2][as.numeric(factor)==i],cex=sizelabel,label[as.numeric(factor)==i],col=col[as.numeric(factor)==i])
        par(new=T)
      }}
    par(new=F)
  }else{
    
    
    plot3D(cbind(matrix[,1],matrix[,2],matrix[,3]),col=col,1,axes=F,bbox=F,cex=cex,add=add)
    
    for (i in 1 :max(as.numeric(factor))){
      lines3d(cbind(matrix[,1][as.numeric(factor)==i],matrix[,2][as.numeric(factor)==i],matrix[,3][as.numeric(factor)==i]),col=col[as.numeric(factor)==i],add=T)
      
      
      
      if(!is.null(label)&sizelabel>0){
        rgl:::text3d(cbind(matrix[,1][as.numeric(factor)==i],matrix[,2][as.numeric(factor)==i],matrix[,3][as.numeric(factor) ==i]),col=col[as.numeric(factor)==i],text=label[as.numeric(factor)==i],cex=rep(sizelabel,length(matrix[,1][as.numeric(factor)==i])))
      }}
    
  }
  
}

