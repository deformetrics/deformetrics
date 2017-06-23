#' animdefo
#'
#' This function saves images of  2D the deformations of a reference 2d shape relatively to  sequence(s) of shapes contained in "seqlist" argument; this can store from 1  to 4 sequences of shapes. Plots are saved in the working directory and can be combined in animated GIF using external softwares.  By default the titles of each plot are "PC1", "PC2", "PC3", "PC4" but any title can be set (see examples). Heatmap is computed by default.
#' @param mshape matrix: a kx2 matrix which is the reference matrix whose deformation is evaluated relatively to the shapes presented in "seqlist" argument. If links are provided the heatmap is computed in the constrained Delaunay triangulation defined by the links. Otherwise the Delaunay triangulation is not constrained.
#' @param seqlist	list: a list with 1 to 4 elements each containing a sequence of shapes that are contrasted with the reference shape. If more than one sequence is provided the sequences must contain the same number of shapes. 
#' @param suffix characther: suffix name to attch to the saved .png files (default = "anim")
#' @param czlim logical: if more than one sequence is present in "seqlist" argument the zlim of all plot is common by default. 
#' @param czlimw logical: if TRUE and If more than one sequence is present in "seqlist" argument the zlim of all plot is different and common only within any specific sequence.  
#' @param linkss numeric:	links for source configuration; if heat=T these (as well as "linkst" argument) must identify a closed contour in order to allow the Delaunay triangulation constrained within the external border. 
#' @param linkst numeric:	links for target configuration
#' @param colss numeric: landmarks color for source
#' @param colst numeric: landmarks color for target
#' @param lwds numeric:	links width for source
#' @param lwdt numeric: links width for target
#' @param cexs numeric:	cex parameter for landmarks of source
#' @param cext numeric:	cex parameter for landmarks of target
#' @param displ logical: if TRUE displacement vectors are displayed
#' @param axtit character: if NULL it uses as defaults: c("PC1+", "PC1-","PC2+", "PC2-","PC3+", "PC3-"). 
#' @param mag numeric: magnification parameter for the illustrated deformation
#' @param heat logical: if TRUE the heatmap is used. Heatmap is computed by calculating the log(determinant) of Jacobian matrix calculated using the first derivative of Thin Plate Spline function.
#' @param mar numeric vector: "mar" argument in par()
#' @param mai numeric vector:	"mai" argument in par()
#' @param oma numeric vector:	"oma" argument in par()
#' @param colors character:	color palette for heatmap
#' @param alpha numeric: transparency parameter for heatmpa color
#' @param magrange numeric:	zoom parameter for subplots. When this parameter is smaller the subplot is bigger. 
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
#' @param pholes list:	an optional list of vectors indicating points identifying holes in the input geometries
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' ### only one plot 
#' data(macrogroup)
#' data(my2d)
#' linksdors<-list(c(1,2),c(37,7),c(12,4),c(27,28),c(25,21),c(38,40),c(9,10),c(2,3),c(3,4),c(1,7),c(1,6),c(3,5),c(6,40),c(5,9),c(40,8),c(8,9),c(1,7),c(7,6),c(3,4),c(4,5),c(39,38),c(38,35),c(35,37),c(37,39),c(35,34),c(34,33),c(33,32),c(32,31),c(31,30),c(30,29),c(29,37),c(37,36),c(36,29),c(28,31),c(28,30),c(13,10),c(10,11),c(11,12),c(12,13),c(13,14),c(14,16),c(16,17),c(17,20),c(20,19),c(19,18),c(18,12),c(18,15),c(15,12),c(21,19),c(21,20),c(24,25),c(25,26),c(26,27),c(27,24),c(26,24),c(24,23),c(23,22),c(22,8),c(8,2))                                                                                                                                                           
#' amy2d<-procSym(my2d,CSinit=T,scale=F)
#' plot(amy2d$PCscores[,1:2],asp=1)
#' pc1seq<-list(scoresequencearray(amy2d$mshape,min(amy2d$PCscores[,1]),max(amy2d$PCscores[,1]),1,amy2d$PCs))
#' animdefo(amy2d$mshape,pc1seq,linkss=linksdors,linkst=linksdors)
#' #### multiplot with common zlim value
#' pc1seq<-scoresequencearray(amy2d$mshape,min(amy2d$PCscores[,1]),max(amy2d$PCscores[,1]),1,amy2d$PCs)
#' pc2seq<-scoresequencearray( amy2d$mshape,min( amy2d$PCscores[,2]),max( amy2d$PCscores[,2]),2, amy2d$PCs)
#' pc3seq<-scoresequencearray( amy2d$mshape,min( amy2d$PCscores[,3]),max( amy2d$PCscores[,3]),3, amy2d$PCs)
#' pc4seq<-scoresequencearray( amy2d$mshape,min( amy2d$PCscores[,4]),max( amy2d$PCscores[,4]),4, amy2d$PCs)
#' animdefo(amy2d$mshape,list(pc1seq,pc2seq,pc3seq,pc4seq),czlim=T,linkss=linksdors,linkst=linksdors)
#' #### plot for allometry indicating size values associated to shape change
#' ##### Shape change associated to evolutionary allometry 
#' regr<-shaperegr(amy2d$orpdata,amy2d$size,links=linksdors,heatmap=T) ### 
#' myaxtits<-list(paste("Centroid Size=",round(regr$myseq,2)))###
#' #### salva i frames dell'animazione
#' animdefo(amy2d$mshape,list(regr$seqshapes),czlim=F,czlimw=T,a=0.01,linkss=linksdors,linkst=linksdors,axtit=myaxtits,magrange=1.5)
#' }
#' @export

animdefo<-function(mshape,seqlist,suffix="anim",czlim=T,czlimw=F,linkss=NULL,linkst=NULL,colss=1,colst=2,lwds=2,lwdt=2,cexs=0,cext=0,displ=F,axtit=NULL,mag=1,heat=T,mar=c(0,0,0.3,0),mai=c(0,0,0.3,0),oma=c(0,0,3,0),colors=c("blue4","cyan2","yellow","red4"),alpha=1,magrange=1,exts=NULL,PB=NA,PA=NA,S=NA,SB=NA,H=NA,V=3,a=NULL,q=NULL,Y=FALSE,j=FALSE,D=FALSE,St=Inf,Q=TRUE,constr=T,pholes=NA){
  allshapes<-NULL
  for(i in 1:length(seqlist)){
    allshi<-seqlist[[i]]
    allshapes<-abind::abind(allshapes,allshi)
  }
  allshapes<-abind::abind(mshape,allshapes)
  allcss<-apply(allshapes,3,cSize)
  
  themax<-allshapes[,,which.max(allcss)]
  
  if(length(seqlist)<2){mfrow=c(1,1)}
  if(length(seqlist)==2){mfrow=c(1,2)}
  if(length(seqlist)==3){mfrow=c(2,2)}
  if(length(seqlist)==4){mfrow=c(2,2)}
  
  if(is.null(exts)==T){exts=rep(0,length(seqlist))}else{exts=exts}
  
  if(!is.null(linkss)&is.na(S)==T){S<-list2matrix(linkss)}
  
  if(is.null(linkss)==T){linkss<-conslinks(1)}
  if(is.null(linkst)==T){linkst<-conslinks(1)}
  
  if(is.null(axtit)==T){
    
    axtit<-NULL
    for(s in 1:length(seqlist)){
      axtits<-paste("PC",rep(s,dim(seqlist[[1]])[3]),sep="")
      axtit<-c(axtit,list(axtits))
    }
    
  }else{
    if(is.character(axtit)==T&abs(length(axtit)-length(seqlist))>0){stop("axtit must be a vector long as the number of elements in seqlist or a list where each element is long as each element in seqlist")}
    if(is.character(axtit)&length(axtit)==length(seqlist)){
      axtitlist<-NULL
      for(s in 1:length(seqlist)){
        axtits<-paste(axtit[s],rep(s,dim(seqlist[[1]])[3]),sep="")
        axtitlist<-c(axtitlist,list(axtits))
      }
      axtit<-axtitlist}else{axtit<-axtit}
  }
  
  
  
  
  
  
  
  if(heat==F){
    if(length(seqlist)==4){
      for(i in 1:dim(seqlist[[j]])[3]){
        png(filename=paste(suffix,i,".png",sep=""))
        par(mfrow=mfrow,mar=mar,oma=oma,mai=mai)
        tpsgridpaolo(seqlist[[1]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[1]][i])
        
        tpsgridpaolo(seqlist[[2]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[2]][i])
        
        tpsgridpaolo(seqlist[[3]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[3]][i])
        
        tpsgridpaolo(seqlist[[4]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[4]][i])
        dev.off()
      }}
    
    if(length(seqlist)==3){
      for(i in 1:dim(seqlist[[1]])[3]){
        png(filename=paste(suffix,i,".png",sep=""))
        par(mfrow=mfrow,mar=mar,oma=oma,mai=mai)
        tpsgridpaolo(seqlist[[1]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[1]][i])
        
        tpsgridpaolo(seqlist[[2]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[2]][i])
        
        tpsgridpaolo(seqlist[[3]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[3]][i])
        
        dev.off()
      }
    }
    if(length(seqlist)==2){
      for(i in 1:dim(seqlist[[1]])[3]){
        png(filename=paste(suffix,i,".png",sep=""))
        par(mfrow=mfrow,mar=mar,oma=oma,mai=mai)
        tpsgridpaolo(seqlist[[1]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[1]][i])
        tpsgridpaolo(seqlist[[2]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[2]][i])
        
        dev.off()
      }
    }
    
    if(length(seqlist)==1){
      for(i in 1:dim(seqlist[[1]])[3]){
        png(filename=paste(suffix,i,".png",sep=""))
        par(mfrow=mfrow,mar=mar,oma=oma,mai=mai)
        tpsgridpaolo(seqlist[[1]][,,i],mshape,linksTT=linkss,linksYY=linkst,lwdtt=lwds,lwdyy=lwdt,axes2d=F,ext=exts[1],mag=mag,collinksTT=colss,collinksYY=colst,displ=displ,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        title(axtit[[1]][i])
        dev.off()
      }
    }
    
    
  }else{
    
    if(czlim==T&czlimw==T){stop("Or czlim=T or czlimw=T; not both")}
    
    if(czlim==T){
      if(length(seqlist)==4){
        myh1neg<-heat2d(mshape,seqlist[[1]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh1pos<-heat2d(mshape,seqlist[[1]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh2neg<-heat2d(mshape,seqlist[[2]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred 
        myh2pos<-heat2d(mshape,seqlist[[2]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh3neg<-heat2d(mshape,seqlist[[3]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh3pos<-heat2d(mshape,seqlist[[3]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh4neg<-heat2d(mshape,seqlist[[4]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh4pos<-heat2d(mshape,seqlist[[4]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        zlim<-range(myh1neg,myh1pos,myh2neg,myh2pos,myh3neg,myh3pos,myh4neg,myh4pos)
      }
      if(length(seqlist)==3){
        myh1neg<-heat2d(mshape,seqlist[[1]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh1pos<-heat2d(mshape,seqlist[[1]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh2neg<-heat2d(mshape,seqlist[[2]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred 
        myh2pos<-heat2d(mshape,seqlist[[2]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh3neg<-heat2d(mshape,seqlist[[3]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh3pos<-heat2d(mshape,seqlist[[3]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        zlim<-range(myh1neg,myh1pos,myh2neg,myh2pos,myh3neg,myh3pos)
      }
      if(length(seqlist)==2){
        myh1neg<-heat2d(mshape,seqlist[[1]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh1pos<-heat2d(mshape,seqlist[[1]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh2neg<-heat2d(mshape,seqlist[[2]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred 
        myh2pos<-heat2d(mshape,seqlist[[2]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        zlim<-range(myh1neg,myh1pos,myh2neg,myh2pos)
      }
      if(length(seqlist)==1){
        myh1neg<-heat2d(mshape,seqlist[[1]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh1pos<-heat2d(mshape,seqlist[[1]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        zlim<-range(myh1neg,myh1pos)
      }
    }else{zlim<-NULL}
    
    
    if(czlimw==T){
      if(length(seqlist)==4){
        myh1neg<-heat2d(mshape,seqlist[[1]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh1pos<-heat2d(mshape,seqlist[[1]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh2neg<-heat2d(mshape,seqlist[[2]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred 
        myh2pos<-heat2d(mshape,seqlist[[2]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh3neg<-heat2d(mshape,seqlist[[3]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh3pos<-heat2d(mshape,seqlist[[3]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh4neg<-heat2d(mshape,seqlist[[4]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh4pos<-heat2d(mshape,seqlist[[4]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        zlim1<-range(myh1neg,myh1pos)
        zlim2<-range(myh2neg,myh2pos)
        zlim3<-range(myh3neg,myh3pos)
        zlim4<-range(myh4neg,myh4pos)
      }
      if(length(seqlist)==3){
        myh1neg<-heat2d(mshape,seqlist[[1]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh1pos<-heat2d(mshape,seqlist[[1]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh2neg<-heat2d(mshape,seqlist[[2]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred 
        myh2pos<-heat2d(mshape,seqlist[[2]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh3neg<-heat2d(mshape,seqlist[[3]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh3pos<-heat2d(mshape,seqlist[[3]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        zlim1<-range(myh1neg,myh1pos)
        zlim2<-range(myh2neg,myh2pos)
        zlim3<-range(myh3neg,myh3pos)
      }
      if(length(seqlist)==2){
        myh1neg<-heat2d(mshape,seqlist[[1]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh1pos<-heat2d(mshape,seqlist[[1]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh2neg<-heat2d(mshape,seqlist[[2]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred 
        myh2pos<-heat2d(mshape,seqlist[[2]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        zlim1<-range(myh1neg,myh1pos)
        zlim2<-range(myh2neg,myh2pos)
      }
      if(length(seqlist)==1){
        myh1neg<-heat2d(mshape,seqlist[[1]][,,1],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        myh1pos<-heat2d(mshape,seqlist[[1]][,,dim(seqlist[[1]])[3]],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=NULL)$pred
        zlim1<-range(myh1neg,myh1pos)
      }
    }
    
    
    
    
    
    if(length(seqlist)==4){
      for(i in 1:dim(seqlist[[1]])[3]){
        png(filename=paste(suffix,i,".png",sep=""))
        par(mfrow=mfrow,mar=mar,oma=oma,mai=mai)
        myh<-heat2d(mshape,seqlist[[1]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim1)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[1]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        
        
        
        myh<-heat2d(mshape,seqlist[[2]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim2)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[2]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        
        myh<-heat2d(mshape,seqlist[[3]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim3)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[3]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        
        
        
        myh<-heat2d(mshape,seqlist[[4]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim4)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[4]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        dev.off()
      }
      
    }   
    
    
    
    if(length(seqlist)==3){
      for(i in 1:dim(seqlist[[1]])[3]){
        png(filename=paste(suffix,i,".png",sep=""))
        par(mfrow=mfrow,mar=mar,oma=oma,mai=mai)
        myh<-heat2d(mshape,seqlist[[1]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim1)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[1]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        
        
        
        myh<-heat2d(mshape,seqlist[[2]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim2)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[2]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        
        myh<-heat2d(mshape,seqlist[[3]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim3)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[3]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        dev.off()
      }
      
    }
    
    
    if(length(seqlist)==2){
      for(i in 1:dim(seqlist[[1]])[3]){
        png(filename=paste(suffix,i,".png",sep=""))
        par(mfrow=mfrow,mar=mar,oma=oma,mai=mai)
        myh<-heat2d(mshape,seqlist[[1]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim1)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[1]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        
        myh<-heat2d(mshape,seqlist[[2]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim2)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[2]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        dev.off()
      }
    }    
    
    if(length(seqlist)==1){
      for(i in 1:dim(seqlist[[1]])[3]){
        png(filename=paste(suffix,i,".png",sep=""))
        par(mfrow=mfrow,mar=mar,oma=oma,mai=mai)
        myh<-heat2d(mshape,seqlist[[1]][,,i],nadd=5000,alpha=alpha,ngrid=22,mag=mag,tol=10,graphics=F,colors=colors,ext=exts[1],PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,constr=constr,pholes=pholes,linkss=linkss,zlim=zlim) 
        if(czlimw==T){image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim1)}else{
          image.plot(xyz2img(cbind(myh$interpcoords,myh$pred),tolerance=myh$tol),asp=1,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,col=myh$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)} 
        par(new=T)
        plot(myh$mate,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkss)==F){lineplot(myh$mate,linkss,col=colss,lwd=lwds)}
        lines(myh$tpsgrid$grid$ngrid,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        lines(myh$tpsgrid$grid$ngrid2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange)
        par(new=T)
        plot(myh$mate2,xlim=range(themax[,1])*magrange,ylim=range(themax[,2])*magrange,asp=1,pch=19,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(linkst)==F){lineplot(myh$mate2,linkst,col=colst,lwd=lwdt)}
        title(axtit[[1]][i])
        if(!is.null(pholes)){
          for(b in 1:length(myh$pols)){
            polygon(myh$pols[[b]],col="white",border=colss)
            polygon(myh$pols2[[b]],col="white",border=colst)}
        }
        dev.off()
      }
      
    }   
  }
}
