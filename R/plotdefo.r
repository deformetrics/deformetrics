#' plotdefo
#'
#' This function displays in 2D the deformations associated typically to a procsym object from procSym() function in package "Morpho" (Stefan Schlager (2016). Morpho: Calculations and Visualisations  Related to Geometric Morphometrics. R package version 2.4.1.1.https://github.com/zarquon42b/Morphor). However, the input list could contain any type of reference matrix, ordination scores and eigenvectors named "mshape" "PCscores" and "PCs" respectively. 
#' @param procsymobject list:	a list with three objects with the following names: "mshape" a kx2 matrix which is deformed according to the values and eigenvectors specified in the other two arguments; "PCscores" ordination scores coming from an ordination analysis; "PCs" eigenvectors relative to "mshape" and "PCscores". By default the function takes the first three PCscores columns and the first three PCs column to build the plot. The list can be built manually in order to contain the desired PCscores and PCs combination. If links are provided the heatmap is computed in the constrained Delaunay triangulation defined by the links. Otherwise the Delaunay triangulation is not constrained. 
#' @param zlim numeric:	zlim range for heat (if TRUE) color.
#' @param linkss numeric:	links for source configuration; if heat=T these (as well as "linkst" argument) must identify a closed contour in order to allow the Delaunay triangulation constrained within the external border. 
#' @param linkst numeric:	links for target configuration
#' @param colss numeric:	landmarks color for source
#' @param colst numeric: landmarks color for target
#' @param lwds numeric: links width for source
#' @param lwdt numeric: links width for target
#' @param cexs numeric: cex parameter for landmarks of source
#' @param cext numeric: cex parameter for landmarks of target
#' @param displ logical: if TRUE displacement vectors are displayed
#' @param axtit character: if NULL it uses as defaults: c("PC1+", "PC1-","PC2+", "PC2-","PC3+", "PC3-"). 
#' @param mag numeric: magnification parameter for the illustrated deformation
#' @param heat logical: if TRUE the heatmap is used. Heatmap is computed by calculating the log(determinant) of Jacobian matrix calculated using the first derivative of Thin Plate Spline function.
#' @param mfrow numeric vector: "mfrow" argument in par() tipically this default serves to display meaningfully positive and negative extremes of ordination scores
#' @param mar numeric vector: "mar" argument in par()
#' @param mai numeric vector:	"mai" argument in par()
#' @param oma numeric vector:	"oma" argument in par()
#' @param colors numeric:	color palette for heatmap
#' @param alpha numeric:	transparency parameter for heatmpa color
#' @param magrange numeric: zoom parameter for subplots. When this parameter is smaller the subplot is bigger. 
#' @param exts numeric vector:	"ext" parameter in tpsgrid() function from "shapes" package
#' @param PB numeric vector: "PB" argument from pslg() function in RTriangle package. 
#' @param PA matrix: "PA" argument from pslg() function in RTriangle package.
#' @param S matrix: "S" argument from pslg() function in RTriangle package.
#' @param SB numeric vector: "SB" argument from pslg() function in RTriangle package.
#' @param H matrix: "H" argument from pslg() function in RTriangle package.
#' @param V numeric: "V" argument from triangulate() function in RTriangle package.
#' @param a numeric: "a" argument from triangulate() function in RTriangle package.
#' @param q numeric: "q" argument from triangulate() function in RTriangle package.
#' @param Y logical:	"Y" argument from triangulate() function in RTriangle package.
#' @param j logical:		"j" argument from triangulate() function in RTriangle package.
#' @param D logical:	"D" argument from triangulate() function in RTriangle package.
#' @param St numeric:		"S" argument from triangulate() function in RTriangle package.
#' @param Q logical: "Q" argument from triangulate() function in RTriangle package.
#' @param constr logical:	if TRUE a constrained triangulation according to links is performed
#' @param pholes list:	an optional list of vectors indicating points identifying holes in the input geoemtries
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' data(macrogroup)
#' data(my2d)
#' linksdors<-list(c(1,2),c(37,7),c(12,4),c(27,28),c(25,21),c(38,40),c(9,10),c(2,3),c(3,4),c(1,7),c(1,6),c(3,5),c(6,40),c(5,9),c(40,8),c(8,9),c(1,7),c(7,6),c(3,4),c(4,5),c(39,38),c(38,35),c(35,37),c(37,39),c(35,34),c(34,33),c(33,32),c(32,31),c(31,30),c(30,29),c(29,37),c(37,36),c(36,29),c(28,31),c(28,30),c(13,10),c(10,11),c(11,12),c(12,13),c(13,14),c(14,16),c(16,17),c(17,20),c(20,19),c(19,18),c(18,12),c(18,15),c(15,12),c(21,19),c(21,20),c(24,25),c(25,26),c(26,27),c(27,24),c(26,24),c(24,23),c(23,22),c(22,8),c(8,2))                                                                                                                                                           
#' amy2d<-procSym(my2d,CSinit=T,scale=F)
#' plot(amy2d$PCscores[,1:2],asp=1)
#' plotdefo(amy2d,linkss=linksdors,linkst=linksdors)
#' orbitR<-centroids(amy2d$mshape[c(18,19,20,17,16,14,13,12,15),])
#' orbitL<-centroids(amy2d$mshape[c(37,35,34,33,32,31,30,29,36),])
#' plotdefo(amy2d,linkss=linksdors,linkst=rbind(orbitR,orbitL),pholes=list(c(18,19,20,17,16,14,13,12,15),c(37,35,34,33,32,31,30,29,36)))
#' }
#' @export
plotdefo<-function(procsymobject,zlim=NULL,linkss=NULL,linkst=NULL,colss=1,colst=2,lwds=2,lwdt=2,cexs=0,cext=0,displ=F,axtit=NULL,mag=1,heat=T,mfrow=c(3,2),mar=c(0,0,0.3,0),mai=c(0,0,0.3,0),oma=c(0,0,3,0),colors=c("blue4","cyan2","yellow","red4"),alpha=1,magrange=1,exts=rep(0,6),PB=NA,PA=NA,S=NA,SB=NA,H=NA,V=3,a=NULL,q=NULL,Y=FALSE,j=FALSE,D=FALSE,St=Inf,Q=TRUE,constr=T,pholes=NA){
  
  mshape<-procsymobject$mshape
  pc1pos<-showPC(max(procsymobject$PCscores[,1]),procsymobject$PCs[,1],procsymobject$mshape)
  pc1neg<-showPC(min(procsymobject$PCscores[,1]),procsymobject$PCs[,1],procsymobject$mshape)
  pc2pos<-showPC(max(procsymobject$PCscores[,2]),procsymobject$PCs[,2],procsymobject$mshape)
  pc2neg<-showPC(min(procsymobject$PCscores[,2]),procsymobject$PCs[,2],procsymobject$mshape)
  pc3pos<-showPC(max(procsymobject$PCscores[,3]),procsymobject$PCs[,3],procsymobject$mshape)
  pc3neg<-showPC(min(procsymobject$PCscores[,3]),procsymobject$PCs[,3],procsymobject$mshape)
  
  allcss<-c(cSize(pc1pos),cSize(pc1neg),cSize(pc2pos),cSize(pc2neg),cSize(pc3pos),cSize(pc3neg))
  themax<-list(pc1pos,pc1neg,pc2pos,pc2neg,pc3pos,pc3neg)[which.max(allcss)][[1]]
  
  
  par(mfrow=mfrow,mar=mar,oma=oma,mai=mai)
  if(!is.null(linkss)&is.na(S)==T){S<-list2matrix(linkss)}
  
  if(is.null(linkss)==T){linkss<-conslinks(1)}
  if(is.null(linkst)==T){linkst<-conslinks(1)}
  
  if(is.null(axtit)==T){axtit<-c("PC1+","PC1-","PC2+","PC2-","PC3+","PC3-")}else{axtit=axtit}
  
  if(heat==F){
    ##### plot deformations associated to PC extremes; in red the Grand Mean of the gpa
    ### PC1+ 
    tpsgridpaolo(pc1pos,mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    title(axtit[1])
    ### PC1-
    tpsgridpaolo(pc1neg,mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[2],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    
    title(axtit[2])
    ### PC2+
    tpsgridpaolo(pc2pos,mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[3],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    title(axtit[3])
    ### PC2-
    tpsgridpaolo(pc2neg,mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[4],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    title(axtit[4])
    ### PC3+ 
    tpsgridpaolo(pc3pos,mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[5],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    title(axtit[5])
    ### PC3-
    tpsgridpaolo(pc3neg,mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[6],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    title(axtit[6])
  }else{
    
    myh<-heat2d(mshape,pc1pos,nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
    image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)  
    par(new=T)
    plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
    lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    par(new=T)
    plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
    title(axtit[1])
    if(!is.null(pholes)){
      for(i in 1:length(myh$pols)){
        polygon(myh$pols[[i]],col="white",border=colss)
        polygon(myh$pols2[[i]],col="white",border=colst)}
    }
    
    
    myh<-heat2d(mshape,pc1neg,nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[2],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
    image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)  
    par(new=T)
    plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
    lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    par(new=T)
    plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
    title(axtit[2])
    if(!is.null(pholes)){
      for(i in 1:length(myh$pols)){
        polygon(myh$pols[[i]],col="white",border=colss)
        polygon(myh$pols2[[i]],col="white",border=colst)}
    }
    
    
    myh<-heat2d(mshape,pc2pos,nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[3],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
    image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)  
    par(new=T)
    plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
    lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    par(new=T)
    plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
    title(axtit[3])
    if(!is.null(pholes)){
      for(i in 1:length(myh$pols)){
        polygon(myh$pols[[i]],col="white",border=colss)
        polygon(myh$pols2[[i]],col="white",border=colst)}
    }
    
    myh<-heat2d(mshape,pc2neg,nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[4],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
    image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)  
    par(new=T)
    plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
    lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    par(new=T)
    plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
    title(axtit[4])
    if(!is.null(pholes)){
      for(i in 1:length(myh$pols)){
        polygon(myh$pols[[i]],col="white",border=colss)
        polygon(myh$pols2[[i]],col="white",border=colst)}
    }
    
    myh<-heat2d(mshape,pc3pos,nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[5],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
    image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)  
    par(new=T)
    plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
    lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    par(new=T)
    plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
    title(axtit[5])
    if(!is.null(pholes)){
      for(i in 1:length(myh$pols)){
        polygon(myh$pols[[i]],col="white",border=colss)
        polygon(myh$pols2[[i]],col="white",border=colst)}
    }
    
    myh<-heat2d(mshape,pc3neg,nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[6],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
    image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)  
    par(new=T)
    plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
    lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
    par(new=T)
    plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
    title(axtit[6])
    if(!is.null(pholes)){
      for(i in 1:length(myh$pols)){
        polygon(myh$pols[[i]],col="white",border=colss)
        polygon(myh$pols2[[i]],col="white",border=colst)}
    }
  }
}