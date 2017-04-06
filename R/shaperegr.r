#' shaperegr
#'
#' This function performs the multivariate regression between an array of shapes (2D or 3D) and a dependent variable or a matrix of independent variables. It also plot shapes predicted  at low and high independent variable values. If the independent variable is univariate Canonical Correlation Analysis (with optional group-structure) is also displayed. Optionally, the heatmap on deformation between shapes predicted at low and high independent variable values is computed. The function saves also a sequence of shapes predicted at equally spaced values (20 by default) within the range of independent variable.
#' @param shapearray array:	an array of shapes treated as dependent variable indep	a vector representing the independent variable or a matrix of independent variables. In this latter case the shapes are predicted at low and high values of each variable present in the matrix. 
#' @param mag numeric: magnification parameter for deformation visualization
#' @param frames numeric:	number of shapes predicted at equally spaced values within the range of independent variable
#' @param links numeric: links structure 
#' @param zlim numeric:	range of heatmap color map for 2D visualization. 
#' @param colcca numeric: colors for points in CCA plot
#' @param legend character:	legend for group structure
#' @param pchcca numeric:	pch symbols in CCA plot
#' @param lwd numeric: links width
#' @param heatmap logical: if TRUE the 2D heatmap color is displayed
#' @param triang list:	for 3D data an optional triangulation structure that is used for computing heatmap in 3D
#' @param group numeric:	group structure to be visualized in CCA plot
#' @param rampcolors character: color palette for heatmap
#' @param alpha numeric:	Transparency parameter for 3D visualization
#' @param from numeric:	Low range value for heatmap visualization in 3D
#' @param to numeric:	High range value for heatmap visualization in 3D
#' @return predmin matrix: shape predicted at low values of independent variable(s)
#' @return predmax matrix: predicted at high values of independent variable(s)
#' @return seqshapes array: sequence of shapes predicted at equally spaced values within the range of independent variable(s)
#' @return myseq numeric vector: equally spaced values of independent variable(s) at which shapes are predicted
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' ### only one plot 
#' data(macrogroup)
#' data(my2d)
#' linksdors<-list(c(1,2),c(37,7),c(12,4),c(27,28),c(25,21),c(38,40),c(9,10),c(2,3),c(3,4),c(1,7),c(1,6),c(3,5),c(6,40),c(5,9),c(40,8),c(8,9),c(1,7),c(7,6),c(3,4),c(4,5),c(39,38),c(38,35),c(35,37),c(37,39),c(35,34),c(34,33),c(33,32),c(32,31),c(31,30),c(30,29),c(29,37),c(37,36),c(36,29),c(28,31),c(28,30),c(13,10),c(10,11),c(11,12),c(12,13),c(13,14),c(14,16),c(16,17),c(17,20),c(20,19),c(19,18),c(18,12),c(18,15),c(15,12),c(21,19),c(21,20),c(24,25),c(25,26),c(26,27),c(27,24),c(26,24),c(24,23),c(23,22),c(22,8),c(8,2))                                                                                                                                                           
#' amy2d<-procSym(my2d)
#' shapearray<-procSym(my2d,CSinit=T,scale=F)$orpdata#### prova cambiando CSinit per il size and shape space
#' indep<-amy2d$size
#' mag=1
#' frames=20
#' links=linksdors
#' col=1
#' lwd=2
#' group=macrogroup
#' shaperegr(shapearray,indep,links=links)
#' shaperegr(shapearray,indep,links=links,group=group,colcca=as.numeric(group))
#' shaperegr(shapearray,indep,links=links,group=group,heatmap=T)
#' ##  in 3 dimensions
#' data(pri3d)
#' data(sur_ent)
#' data(linksbase)
#' data(linksface)
#' data(linksentire)
#' data=pri3d
#' my3d<-centershapes(data)
#' amy3d<-procSym(my3d)
#' shapearray<-procSym(my3d,CSinit=T,scale=F)$orpdata####provacambiandoCSinitperilsizeandshapespace
#' indep<-amy3d$size
#' triang<-t(sur_ent$it)
#' group<-factor(substr(dimnames(amy3d$orpdata)[[3]],1,7))
#' shaperegr(shapearray,indep,links=linksentire)
#' shaperegr(shapearray,indep,links=linksentire,group=group,colcca=as.numeric(group),pchcca=as.numeric(group))
#' prov<-shaperegr(shapearray,indep,links=linksentire,group=group,heatmap=T,triang=triang)
#' prov2<-shaperegr(shapearray,cbind(amy3d$size,rnorm(length(indep),0,1)),links=linksentire,group=group,heatmap=T,triang=triang)
#' }
#' @export
shaperegr<-function(shapearray,indep,mag=1,frames=20,links=NULL,zlim=NULL,colcca=NULL,legend=T,pchcca=NULL,lwd=2,heatmap=F,triang=NULL,group=NULL,rampcolors=c("blue4","cyan2","yellow","red4"),alpha=0.7,from=NULL,to=NULL){
  if(is.null(links)==T){links=c(1,1)}
  require(Morpho)
  require(shapes)
  require(vegan)
  
  nland<-dim(shapearray)[1]
  ndim<-dim(shapearray)[2]
  nind<-dim(shapearray)[3]
  
  if(ndim>3){
    if(is.null(triang)==F){
      if(ncol(triang)>3){stop("I want triangulation as Mx3 matrix")}
    }
  }
  
  print(RsquareAdj(vegan::rda(array2mat(shapearray,nind,nland*ndim)~indep)))
  print(anova(vegan::rda(array2mat(shapearray,nind,nland*ndim)~indep),step=2000))
  mylm<-lm(array2mat(shapearray,nind,nland*ndim)~indep)
  if(ncol(as.matrix(indep))>1){
    shapemaxindep<-read.inn(t(as.matrix(as.numeric(crossprod(coef(mylm),c(1,apply(indep,2,max)*mag))))),nland,ndim)
    shapeminindep<-read.inn(t(as.matrix(as.numeric(crossprod(coef(mylm),c(1,apply(indep,2,min)*mag))))),nland,ndim)
  }else{  
    shapemaxindep<-read.inn(t(as.matrix(as.numeric(crossprod(coef(mylm),c(1,max(indep)*mag))))),nland,ndim)
    shapeminindep<-read.inn(t(as.matrix(as.numeric(crossprod(coef(mylm),c(1,min(indep)*mag))))),nland,ndim)
  }
  
  allcss<-unlist(lapply(list(shapemaxindep,shapeminindep),cSize))
  themax<-list(shapemaxindep,shapeminindep)[which.max(allcss)][[1]][,,1]
  
  if(is.null(colcca)==T){colcca<-rep(1,nind)}else{colcca=colcca}
  if(is.null(pchcca)==T){pchcca<-rep(19,nind)}else{pchcca=pchcca}
  
  if(ndim<3){
    par(mfrow=c(2,2))
    if(heatmap==F){
      plotmyarrays(shapeminindep[,,1],xlim=range(themax[,1]),ylim=range(themax[,2]),links=links,txt=F,pch=19)
      title(main="At low x-values")
      plotmyarrays(shapemaxindep[,,1],xlim=range(themax[,1]),ylim=range(themax[,2]),links=links,txt=F,pch=19)
      title(main="At high x-values")
      tpsgridpaolo(shapeminindep[,,1],shapemaxindep[,,1],opt=1,ngrid=20,,linksYY=links,displ=F,mag=mag,linksTT=links,collinksTT=2,collinksYY=1,xlim=range(themax[,1]),ylim=range(themax[,2]),ext=0)
      title(main="Black: at high x-values")
      if(ncol(as.matrix(indep))<2){mycca(indep,array2mat(shapearray,nind,ndim*nland),legend=legend,group=group,col=colcca,pch=pchcca,xlab=NULL,ylab=NULL,posl="bottomright")
      title("CCA analysis")}
    }else{
      par(mfrow=c(2,2))
      plotmyarrays(shapeminindep[,,1],xlim=range(themax[,1]),ylim=range(themax[,2]),links=links,txt=F,pch=19)
      title(main="At low x-values")
      plotmyarrays(shapemaxindep[,,1],xlim=range(themax[,1]),ylim=range(themax[,2]),links=links,txt=F,pch=19)
      title(main="At high x-values")
      myhxpos<-heat2d(shapeminindep[,,1],shapemaxindep[,,1],tol=10,nadd=5000,graphics=F,constr=T,colors=rampcolors,linkss=links,zlim=zlim) 
      image.plot(xyz2img(cbind(myhxpos$interpcoords,myhxpos$pred),tolerance=myhxpos$tol),asp=1,xlim=range(themax[,1]),ylim=range(themax[,2]),col=myhxpos$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
      par(new=T)
      plot(myhxpos$mate,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="From low (in black) to high x-values (in red)")
      if(is.null(links)==F){lineplot(myhxpos$mate,links,col=2,lwd=lwd)}
      lines(myhxpos$tpsgrid$grid$ngrid,xlim=range(themax[,1]),ylim=range(themax[,2]))
      lines(myhxpos$tpsgrid$grid$ngrid2,xlim=range(themax[,1]),ylim=range(themax[,2]))
      par(new=T)
      plot(myhxpos$mate2,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
      if(is.null(links)==F){lineplot(myhxpos$mate2,links,col=1,lwd=lwd)}
      if(ncol(as.matrix(indep))<2){mycca(indep,array2mat(shapearray,nind,ndim*nland),legend=legend,group=group,col=colcca,pch=pchcca,xlab=NULL,ylab=NULL,,posl="bottomright")
      title("CCA analysis")}
    }
  }
  
  if(ndim>2){
    open3d(windowRect=c(100,100,1000,1000)) 
    mat <- matrix(1:8, ncol=2)
    layout3d(mat, model = "inherit")
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    plot3d(shapeminindep[,,1],bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
    if(!is.null(links)){lineplot(shapeminindep[,,1],links,lwd=1)}
    if(!is.null(triang)){shade3d(plotsurf(shapeminindep[,,1],t(triang),plot=F),alpha=0.7)}
    next3d()
    text3d(0,0,0,"At low x-values")
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    plot3d(shapemaxindep[,,1],bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
    if(!is.null(links)){lineplot(shapemaxindep[,,1],links,lwd=1)}
    if(!is.null(triang)){shade3d(plotsurf(shapemaxindep[,,1],t(triang),plot=F),alpha=0.7)}
    
    next3d()
    text3d(0,0,0,"At high x-values")
    if(heatmap==T){
      if(is.null(triang)==T){stop("Please input triangulation")}
      myhxpos<-diffonmesh(shapeminindep[,,1],shapemaxindep[,,1],t(triang),from=from,to=to,rampcolors=rampcolors,alphas=c(alpha,0.7),graph=F,plotsource=F)
      next3d()
      plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      shade3d(myhxpos$obm$colMesh,alpha=alpha)
      if(!is.null(links)){lineplot(shapemaxindep[,,1],links)}
      next3d()
      text3d(0,0,0,"Heatmap")
    }
    
    if(ncol(as.matrix(indep))<2){mycca(indep,array2mat(shapearray,nind,ndim*nland),legend=legend,group=group,col=colcca,pch=pchcca,xlab=NULL,ylab=NULL,posl="bottomright")
    title("CCA analysis")}
  }
  
  
  if(ncol(as.matrix(indep))<2){myseq<-seq(min(indep),max(indep),length.out=frames)}
    
  if(ncol(as.matrix(indep))>1){
      myseq<-NULL
      for(k in 1:ncol(indep)){
        myseqi<-seq(min(indep[,k]),max(indep[,k]),length.out=frames)
        myseq<-cbind(myseq,myseqi)
      }
    }
  
  
  if(ncol(as.matrix(indep))<2){
    seqshapes<-NULL
    for(i in myseq){
      shapepositionindep<-read.inn(t(as.matrix(as.numeric(crossprod(coef(mylm),c(1,i))))),nland,ndim)
      seqshapes<-abind::abind(seqshapes,shapepositionindep)
    }
  }else{
    seqshapes<-NULL
    for(i in 1:nrow(myseq)){
      shapepositionindep<-read.inn(t(as.matrix(as.numeric(crossprod(coef(mylm),c(1,myseq[i,]))))),nland,ndim)
      seqshapes<-abind::abind(seqshapes,shapepositionindep)
      
    }
  }
  out<-list(predmin=shapeminindep[,,1],predmax=shapemaxindep[,,1],seqshapes=seqshapes,myseq=myseq)
  return(out)
}
