#'plotmyarrays
#'
#' This function plot one or more elements of an array of three diensional or two dimensional shapes
#' @param x numeric: a matrix kxm or an array kxmxn of landmark coordinates
#' @param l numeric: is a vector indicating the number of variables of 1st array dimension to plot (default=c(1:dim(x)[1]))
#' @param v numeric: is a vector indicating the number of variables of 2nd array dimension to plot (default=c(1:dim(x)[2]))
#' @param ind numeric: is a vector indicating which elements of 3rd array dimension should be plotted (default=c(1:dim(x)[3]))
#' @param group character vector: optional grouping variable; shapes will be colored according to group affiliation
#' @param links numeric list: an optional list of vectors (default = NULL) indicating the links between landmarks
#' @param xlim numeric vector: xlim argument in par()
#' @param ylim numeric vector: ylim argument in par()
#' @param col numeric: col argument in par()
#' @param txt logical: if TRUE the text of ordinal landmarks number is plotted
#' @param lwd numeric: lwd argument in par() for links
#' @param pch numeric: pch argument in par() for points in 2D
#' @param cex numeric: cex argument in par() for points
#' @param asp numeric: giving the aspect ratio y/x (default=1)
#' @param cextext numeric: size of the text 
#' @param xlab character: xlab argument in par() (default="")
#' @param ylab character: ylab argument in par() (default="")
#' @param zlab character: zlab argument in par() (default="")
#' @param xaxt character: xaxt argument in par() (default="s")
#' @param yaxt character: yaxt argument in par() (default="s")
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' data(my2d)
#' linksdors<-list(c(1,2),c(37,7),c(12,4),c(27,28),c(25,21),c(38,40),c(9,10),c(2,3),c(3,4),c(1,7),c(1,6),c(3,5),c(6,40),c(5,9),c(40,8),c(8,9),c(1,7),c(7,6),c(3,4),c(4,5),c(39,38),c(38,35),c(35,37),c(37,39),c(35,34),c(34,33),c(33,32),c(32,31),c(31,30),c(30,29),c(29,37),c(37,36),c(36,29),c(28,31),c(28,30),c(13,10),c(10,11),c(11,12),c(12,13),c(13,14),c(14,16),c(16,17),c(17,20),c(20,19),c(19,18),c(18,12),c(18,15),c(15,12),c(21,19),c(21,20),c(24,25),c(25,26),c(26,27),c(27,24),c(26,24),c(24,23),c(23,22),c(22,8),c(8,2))                                                                                                                                                           
#' my2d<-procSym(my2d,CSinit=T,scale=F)
#' pc1seq<-scoresequencearray(my2d$mshape,min(my2d$PCscores[,1]),max(my2d$PCscores[,1]),1,my2d$PCs)
#' par(mfrow=c(4,4))
#' for(i in 1:dim(pc1seq)[3]){
#'   plotmyarrays(pc1seq[,,i],links=linksdors,cex=0,txt=F)
#' }
#' @export 
plotmyarrays<-function(x,l=c(1:dim(x)[1]),v=c(1:dim(x)[2]),ind=c(1:dim(x)[3]),group=NULL,links=NULL,xlim=range(x),ylim=range(x),col=NULL,txt=T,lwd=1,pch=NULL,cex=NULL,asp=1,cextext=1,xlab="",ylab="",zlab="",xaxt="s",yaxt="s"){
  #### x is  a two- or three-dimensional array.
  #### l is a vector indicating the number of variables of 1th array dimension to plot
  #### v is a vector indicating the number of variables of 2th array dimension to plot
  #### ind is a vector indicating which elements of 3th array dimension should be plotted
  #### length(group) must equal length(dim(x)[3]) even if you plot just some individuals (the argument "ind")
  
  
  if(class(x)=="matrix"){x<-array(x,dim=c(dim(x)[1],dim(x)[2],1))}
  
  if(length(v)>2){
    
    if(!is.null(group)){
      
      
      for (i in ind){
        text3d(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),text=l,col=c(as.numeric(group)[i]),add=T,cex=cextext)
      }
      open3d()
      for (i in ind){       
        plot3D(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),col=c(as.numeric(group)[i]),add=T,bbox=F)
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),links,lwd=lwd,col=c(as.numeric(group)[i]),add=T)
        }else{NULL}}
      
      open3d()
      for (i in ind){       
        plot3D(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),col=c(as.numeric(group)[i]),add=T,bbox=F)
        text3d(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),text=l,col=c(as.numeric(group)[i]),add=T,cex=cextext)          
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),links,lwd=lwd,col=c(as.numeric(group)[i]),add=T)
        }else{NULL}
        
      }
      
      
      
    }else{
      
      if(is.null(col)){col=c(1:length(ind))}else{col=col}
      if(is.null(cex)){cex=rep(1,length(ind))}else{cex=cex}
      for (i in ind){
        
        
        text3d(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),text=l,col=col[which(ind==i)],add=T,cex=cextext)
      }
      open3d()
      for (i in ind){  
        
        
        plot3D(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),col=col[i],cex=cex[i],add=T,bbox=F)
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),links,lwd=lwd,col=col[which(ind==i)],add=T)
        }else{NULL}}
      
      
      open3d()
      for (i in ind){    
        
        
        plot3D(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),col=col[i],add=T,bbox=F,cex=cex[i])
        text3d(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),text=l,col=col[which(ind==i)],add=T,cex=cextext)          
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),links,lwd=lwd,col=col[which(ind==i)],add=T)
        }else{NULL}
        
      }
    }
  }
  
  
  
  else{
    
    if(!is.null(group)){
      if(is.null(pch)){pch==rep(1:length(ind))}else{pch=pch}
      if(is.null(cex)){cex=rep(1,length(ind))}else{cex=cex}
      
      for (i in ind){
        par(new=T)
        plot(x[l,v[1],i],x[l,v[2],i],col=c(as.numeric(group)[i]),xlim=xlim,ylim=ylim,pch=pch[i],asp=asp,cex=cex[i],xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt)
        if(txt==T){textxy(x[l,v[1],i],x[l,v[2],i],labs=l,col=c(as.numeric(group)[i]),cex=cextext)}
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i]),links,lwd=lwd,col=c(as.numeric(group)[i]))}else{NULL}
        
      }
      
      
    }else{
      
      if(is.null(col)){col=c(1:length(ind))}else{col=col}
      if(is.null(pch)){pch==rep(1:length(ind))}else{pch=pch}
      if(is.null(cex)){cex=rep(1,length(ind))}else{cex=cex}
      
      
      for (i in ind){
        plot(x[l,v[1],i],x[l,v[2],i],col=col[which(ind==i)],xlim=xlim,ylim=ylim,pch=pch[i],asp=asp,cex=cex[i],xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt)
        if(txt==T){textxy(x[l,v[1],i],x[l,v[2],i],labs=l,col=col[which(ind==i)],cex=cextext)}
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i]),links,lwd=lwd,col=col[which(ind==i)])}else{NULL}
        par(new=T)
      }
    }
  }
  par(new=F)
}

