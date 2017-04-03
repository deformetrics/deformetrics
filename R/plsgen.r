#' plsgen
#'
#' This function uses pls2B() function from "Morpho" package in order to illustrate results of two-blocks partial least squares applied to shape/shape, shape/matrix or matrix/matrix blocks. 3D/3D 2D/2D and 3D/2D visualizations are allowed. Heatmap is optionally computed. If a block is a matrix shape is not illustrated. 
#' @param block1 array: either a matrix of variables or an array of 3D or 2D shapes
#' @param block2 array: either a matrix of variables or an array of 3D or 2D shapes
#' @param plsn numeric: which pair of PLS axes should be plotted
#' @param links1 numeric:	links structure for shape 1
#' @param links2 numeric: links structure for shape 2
#' @param commonref logical: if TRUE shapes are plotted on a common reference system. This is useful to plot results of analyses performed on Size and Shape Space.
#' @param heatmap logical: if TRUE heatmap is computed
#' @param heatcolors character: color palette for heatmap
#' @param triang1 matrix: triangulation structure for shape 1 in 3D
#' @param triang2 matrix: triangulation structure for shape 2 in 3D
#' @param alpha numeric: transparency parameter for heatmap
#' @param S1=NA		"S" argument from pslg() function in RTriangle package for shape 1 (if it is in 2D)
#' @param S2=NA	"S" argument from pslg() function in RTriangle package for shape 1 (if it is in 2D)
#' @param from1 numeric: low range value for heatmap visualization in 3D for shape 1
#' @param to1 numeric: max range value for heatmap visualization in 3D for shape 1
#' @param from2 numeric: low range value for heatmap visualization in 3D for shape 2
#' @param to2 numeric: max range value for heatmap visualization in 3D for shape 2
#' @param rounds numeric: permutation for pls2B() function from "Morpho" package for testing correlation between PLS axes pairs 
#' @param sdx numeric: magnification parameter (in terms of standard deviation unit) for shape 1 visualization
#' @param sdy numeric: magnification parameter (in terms of standard deviation unit) for shape 1 visualization
#' @param labels character: labels to add to PLS plot
#' @param group character: group structure to be plotted in PLS plot and for computing means of PLS axes valeus
#' @param col numeric: points color in PLS plot
#' @param pch numeric: points pch in PLS plot
#' @param colgroup character: groups color in PLS plot
#' @param zlim2d numeeric: range of heatmap in 2D
#' @return thepls list:	A pls2B object from pls2B() function from "Morpho" package
#' @return allxscores matrix:	Block1 scores
#' @return allyscores	matrix: Block2 scores
#' @return xscoresmeans	matrix: if group is provided group means for x scores
#' @return yscoresmeans	matrix: if group is provided group means for y scores
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' data(pri3d)
#' data(surf)
#' data(linksbase)
#' data(linksface)
#' data(linksentire)
#' data(sur_ent)
#' data(sur_fac)
#' data(sur_bas)
#' data(my2d)
#' data=pri3d
#' my3d<-centershapes(data)
#' amy3d<-procSym(my3d)
#' linksdors<-list(c(1,2),c(37,7),c(12,4),c(27,28),c(25,21),c(38,40),c(9,10),c(2,3),c(3,4),c(1,7),c(1,6),c(3,5),c(6,40),c(5,9),c(40,8),c(8,9),c(1,7),c(7,6),c(3,4),c(4,5),c(39,38),c(38,35),c(35,37),c(37,39),c(35,34),c(34,33),c(33,32),c(32,31),c(31,30),c(30,29),c(29,37),c(37,36),c(36,29),c(28,31),c(28,30),c(13,10),c(10,11),c(11,12),c(12,13),c(13,14),c(14,16),c(16,17),c(17,20),c(20,19),c(19,18),c(18,12),c(18,15),c(15,12),c(21,19),c(21,20),c(24,25),c(25,26),c(26,27),c(27,24),c(26,24),c(24,23),c(23,22),c(22,8),c(8,2))                                                                                                                                                           
#' amy2d<-procSym(my2d)
#' block1<-procSym(my2d[,,1:30])$orpdata
#' block2<-procSym(my2d[,,31:60])$orpdata
#' plsgen(block1,block2,links1=linksdors,links2=linksdors,commonref=T)
#' plsgen(block1,block2,links1=linksdors,links2=linksdors,commonref=F)
#' plsgen(block1,block2,links1=linksdors,links2=linksdors,commonref=F,heatmap=T)
#' block1<-procSym(my2d[,,1:30],,CSinit=T,scale=F)$orpdata
#' block2<-procSym(my3d[18:32,,1:30],CSinit=T,scale=F)$orpdata
#' plsgen(block1,block2,links1=linksdors,links2=linksface,heatmap=T,triang2=t(sur_fac$it),sdx=3)
#' block1<-procSym(my3d[1:17,,],CSinit=T,scale=F)$orpdata
#' block2<-procSym(my3d[18:32,,],CSinit=T,scale=F)$orpdata
#' plsgen(block1,block2,links1=linksbase,links2=linksface,commonref=T)
#' plsgen(block1,block2,commonref=T,triang1=t(sur_bas$it),heatmap=T,triang2=t(sur_fac$it),from1=0,from2=0,to1=0.051,to2=0.051)
### Size and shape space 
#' block1<-procSym(my3d[1:17,,],CSinit=F,scale=F)$orpdata
#' block2<-procSym(my3d[18:32,,],CSinit=F,scale=F)$orpdata
#' plsgen(block1,block2,commonref=T,triang1=t(sur_bas$it),heatmap=T,triang2=t(sur_fac$it),from1=0,from2=0,to1=16,to2=16)
#' block1<-matrix(rnorm(1000,0,1),ncol=10)
#' block2<-matrix(rnorm(1000,0,1),ncol=10)
#' plsgen(block1,block2)
#' block1<-procSym(my3d[1:17,,1:30])$orpdata
#' block2<-as.matrix(cbind(procSym(my3d[1:17,,1:30])$size,rnorm(30,0,1)))
#' plsgen(block1,block2,links1=linksbase,heatmap=T,triang1=t(sur_bas$it))
#' }
#' @export

plsgen<-function(block1,block2,plsn=1,links1=NULL,links2=NULL,commonref=F,heatmap=F,heatcolors=c("blue4","cyan2","yellow","red4"),triang1=NULL,triang2=NULL,alpha=1,S1=NA,S2=NA,from1=NULL,to1=NULL,from2=NULL,to2=NULL,rounds=0,sdx=1,sdy=1,labels=NULL,group=NULL,col=1,colgroup=NULL,zlim2d=NULL){
  thepls<-pls2B(block1,block2,rounds=rounds)
  XScores<-thepls$Xscores
  YScores<-thepls$Yscores
  plot(XScores[, plsn],YScores[,plsn],asp=1)
  if(!is.null(labels)){textxy(XScores[, plsn],YScores[, plsn],labels)}
  
  if(!is.null(group)){plot2dhull(cbind(XScores[, plsn],YScores[, plsn]),group,1,col=col,colhull=colgroup,labels=labels)}else{NULL}
  
  if(!is.null(group)){
    xscoresmeans<-as.matrix(aggregate(XScores,by=list(group),mean)[,-1])
    rownames(xscoresmeans)<-aggregate(XScores,by=list(group),mean)[,1]
    yscoresmeans<-as.matrix(aggregate(YScores,by=list(group),mean)[,-1])
    rownames(yscoresmeans)<-aggregate(YScores,by=list(group),mean)[,1]
    x11()
    plot(xscoresmeans[,plsn],yscoresmeans[,plsn])
    textxy(xscoresmeans[,plsn],yscoresmeans[,plsn],rownames(xscoresmeans))
  }else{NULL}
  
  if(!is.null(links1)&is.na(S1)==T){S1<-list2matrix(links1)}
  if(!is.null(links2)&is.na(S2)==T){S2<-list2matrix(links2)}
  
  theshapes<-plsCoVar(thepls,plsn,sdx=sdx,sdy=sdy)
  if(is.matrix(block1)==F){
    plsnxpos<-theshapes$x[,,2]
    plsnxneg<-theshapes$x[,,1]
    css1<-c(cSize(plsnxpos),cSize(plsnxneg))
    mshape1<-arrMean3(block1)
    if(heatmap==T){
      if(dim(plsnxpos)[2]<3){
        if(!is.null(links1)){S1<-list2matrix(links1)}
        myhxpos<-heat2d(mshape1,plsnxpos,tol=10,nadd=5000,graphics=F,constr=T,colors=heatcolors,linkss=links1,S=S1,zlim=zlim) 
        myhxneg<-heat2d(mshape1,plsnxneg,tol=10,nadd=5000,graphics=F,constr=T,colors=heatcolors,linkss=links1,S=S1,zlim=zlim) 
      }else{
        if(is.null(triang1)){stop("You cannot want heatmap in 3d without triangulation for block1")}
        myhxpos<-diffonmesh(mshape1,plsnxpos,t(triang1),from=from1,to=to1,rampcolors=heatcolors,alphas=c(alpha,0.7),graph=F)
        myhxneg<-diffonmesh(mshape1,plsnxneg,t(triang1),from=from1,to=to1,rampcolors=heatcolors,alphas=c(alpha,0.7),graph=F)
      }
    }
    
    
  }else{NULL}
  
  
  if(is.matrix(block2)==F){
    plsnypos<-theshapes$y[,,2]
    plsnyneg<-theshapes$y[,,1]
    css2<-c(cSize(plsnypos),cSize(plsnyneg))
    mshape2<-arrMean3(block2)
    if(heatmap==T){
      if(dim(plsnypos)[2]<3){
        if(!is.null(links2)){S2<-list2matrix(links2)}
        myhypos<-heat2d(mshape2,plsnypos,tol=10,nadd=5000,graphics=F,constr=T,colors=heatcolors,linkss=links2,S=S2) 
        myhyneg<-heat2d(mshape2,plsnyneg,tol=10,nadd=5000,graphics=F,constr=T,colors=heatcolors,linkss=links2,S=S2) 
      }else{
        if(is.null(triang2)){stop("You cannot want heatmap in 3d without triangulation for block2")}
        myhypos<-diffonmesh(mshape2,plsnypos,t(triang2),from=from2,to=to2,rampcolors=heatcolors,alphas=c(alpha,0.7),graph=F)
        myhyneg<-diffonmesh(mshape2,plsnyneg,t(triang2),from=from2,to=to2,rampcolors=heatcolors,alphas=c(alpha,0.7),graph=F)
      }
    }
    
  }else{NULL}
  
  ###### forma matrice 
  if(length(dim(block1))>2&length(dim(block2))<3){
    if(css1[1]<css1[2]){themax=plsnxneg}
    if(css1[1]>css1[2]){themax=plsnxpos}
    if(css1[1]==css1[2]){themax=plsnxpos}
    if(dim(block1)[2]>2){
      open3d(windowRect=c(100,100,1000,1000)) 
      mat <- matrix(1:4, ncol=2)
      layout3d(mat, height = rep(c(2,1), 2), model = "inherit")
      plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      if(heatmap==T){
        shade3d(myhxpos$obm$colMesh,alpha=alpha)
        if(is.null(links1)==F){lineplot(plsnxpos,links1,col=1)}
      }else{
        plot3d(plsnxpos,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links1)==F){lineplot(plsnxpos,links1,col=1)}}
      next3d()
      text3d(0,0,0, "Block1 positive")
      next3d()
      plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      if(heatmap==T){
        shade3d(myhxneg$obm$colMesh,alpha=alpha)
        if(is.null(links1)==F){lineplot(plsnxneg,links1,col=1)}
      }else{
        plot3d(plsnxneg,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links1)==F){lineplot(plsnxneg,links1,col=1)}}
      next3d()
      text3d(0,0,0, "Block1 negative")
    }
    if(dim(block1)[2]<3){
      par(mfrow=c(1,2))
      if(heatmap==T){
        image.plot(xyz2img(cbind(myhxpos$interpcoords,myhxpos$pred),tolerance=myhxpos$tol),asp=1,xlim=range(themax[,1]),ylim=range(themax[,2]),col=myhxpos$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")  
        par(new=T)
        plot(myhxpos$mate,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links1)==F){lineplot(myhxpos$mate,links1)}
        lines(myhxpos$tpsgrid$grid$ngrid,xlim=range(themax[,1]),ylim=range(themax[,2]))
        lines(myhxpos$tpsgrid$grid$ngrid2,xlim=range(themax[,1]),ylim=range(themax[,2]))
        par(new=T)
        plot(myhxpos$mate2,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links1)==F){lineplot(myhxpos$mate2,links1,col=2,lwd=1)}
        title("Block1 positive")
        image.plot(xyz2img(cbind(myhxneg$interpcoords,myhxneg$pred),tolerance=myhxneg$tol),asp=1,xlim=range(themax[,1]),ylim=range(themax[,2]),col=myhxneg$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")  
        par(new=T)
        plot(myhxneg$mate,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=0.5,col=2,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links1)==F){lineplot(myhxneg$mate,links1)}
        lines(myhxneg$tpsgrid$grid$ngrid,xlim=range(themax[,1]),ylim=range(themax[,2]))
        lines(myhxneg$tpsgrid$grid$ngrid2,xlim=range(themax[,1]),ylim=range(themax[,2]))
        par(new=T)
        plot(myhxneg$mate2,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links1)==F){lineplot(myhxneg$mate2,links1,col=2,lwd=1)}
        title("Block1 negative")
        
      }else{
        plotmyarrays(plsnxpos,pch=19,links=links1,txt=F,xlim=range(themax[,1]),ylim=range(themax[,2]))
        title("Block1 positive")
        plotmyarrays(plsnxneg,pch=19,links=links1,txt=F,xlim=range(themax[,1]),ylim=range(themax[,2]))
        title("Block1 negative")
      }
    }
  }
  #########  matrice forma  
  if(length(dim(block1))<3&length(dim(block2))>2){
    if(css2[1]<css2[2]){themax=plsnyneg}
    if(css2[1]>css2[2]){themax=plsnypos}
    if(css2[1]==css2[2]){themax=plsnypos}
    if(dim(block2)[2]>2){
      open3d(windowRect=c(100,100,1000,1000)) 
      mat <- matrix(1:4, ncol=2)
      layout3d(mat, height = rep(c(2,1), 2), model = "inherit")
      plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      if(heatmap==T){
        shade3d(myhypos$obm$colMesh,alpha=alpha)
        if(is.null(links2)==F){lineplot(plsnypos,links2,col=1)}
      }else{
        plot3d(plsnypos,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links2)==F){lineplot(plsnypos,links2,col=1)}}
      next3d()
      text3d(0,0,0, "Block2 positive")
      next3d()
      plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      if(heatmap==T){
        shade3d(myhyneg$obm$colMesh,alpha=alpha)
        if(is.null(links2)==F){lineplot(plsnyneg,links2,col=1)}
      }else{
        plot3d(plsnyneg,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links2)==F){lineplot(plsnypos,links2,col=1)}}
      next3d()
      text3d(0,0,0, "Block2 negative")
    }
    if(dim(block2)[2]<3){
      par(mfrow=c(1,2))
      if(heatmap==T){
        image.plot(xyz2img(cbind(myhypos$interpcoords,myhypos$pred),tolerance=myhypos$tol),asp=1,xlim=range(themax[,1]),ylim=range(themax[,2]),col=myhypos$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")  
        par(new=T)
        plot(myhypos$mate,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links2)==F){lineplot(myhypos$mate,links2)}
        lines(myhypos$tpsgrid$grid$ngrid,xlim=range(themax[,1]),ylim=range(themax[,2]))
        lines(myhypos$tpsgrid$grid$ngrid2,xlim=range(themax[,1]),ylim=range(themax[,2]))
        par(new=T)
        plot(myhypos$mate2,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links2)==F){lineplot(myhypos$mate2,links2,col=2,lwd=2)}
        title("Block2 positive")
        image.plot(xyz2img(cbind(myhyneg$interpcoords,myhyneg$pred),tolerance=myhyneg$tol),asp=1,xlim=range(themax[,1]),ylim=range(themax[,2]),col=myhyneg$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")  
        par(new=T)
        plot(myhyneg$mate,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=0.5,col=2,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links2)==F){lineplot(myhyneg$mate,links2)}
        lines(myhyneg$tpsgrid$grid$ngrid,xlim=range(themax[,1]),ylim=range(themax[,2]))
        lines(myhyneg$tpsgrid$grid$ngrid2,xlim=range(themax[,1]),ylim=range(themax[,2]))
        par(new=T)
        plot(myhyneg$mate2,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links2)==F){lineplot(myhyneg$mate2,links2,col=2,lwd=2)}
        title("Block2 negative")
        
      }else{
        plotmyarrays(plsnypos,pch=19,links=links2,txt=F,xlim=range(themax[,1]),ylim=range(themax[,2]))
        title("Block2 positive")
        plotmyarrays(plsnyneg,pch=19,links=links2,txt=F,xlim=range(themax[,1]),ylim=range(themax[,2]))
        title("Block2 negative")}
    }
  }
  ########  forma forma
  if(length(dim(block1))>2&length(dim(block2))>2){
    if(dim(block1)[2]>2&dim(block2)[2]>2){
      allcss<-c(css1,css2)
      themax<-list(plsnxpos,plsnxneg,plsnypos,plsnyneg)[which.max(allcss)][[1]]
      themaxx<-list(plsnxpos,plsnxneg)[which.max(css1)][[1]]
      themaxy<-list(plsnypos,plsnyneg)[which.max(css2)][[1]]
      open3d(windowRect=c(100,100,1000,1000)) 
      mat <- matrix(1:8, ncol=2)
      layout3d(mat, height = rep(c(2,1), 2), model = "inherit")
      
      if(commonref==T){plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)}else{
        plot3d(themaxx*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      }
      
      if(heatmap==T){
        shade3d(myhxpos$obm$colMesh,alpha=alpha)
        if(is.null(links1)==F){lineplot(plsnxpos,links1,col=1)}
      }else{
        
        plot3d(plsnxpos,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links1)==F){lineplot(plsnxpos,links1,col=1)}}
      next3d()
      text3d(0,0,0, "Block1 positive")
      next3d()
      if(commonref==T){plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)}else{
        plot3d(themaxx*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      }
      if(heatmap==T){
        shade3d(myhxneg$obm$colMesh,alpha=alpha)
        if(is.null(links1)==F){lineplot(plsnxneg,links1,col=1)}
      }else{
        plot3d(plsnxneg,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links1)==F){lineplot(plsnxneg,links1,col=1)}}
      next3d()
      text3d(0,0,0, "Block1 negative")
      next3d()
      if(commonref==T){plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)}else{
        plot3d(themaxy*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      }
      if(heatmap==T){
        shade3d(myhypos$obm$colMesh,alpha=alpha)
        if(is.null(links2)==F){lineplot(plsnypos,links2,col=1)}
      }else{
        plot3d(plsnypos,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links2)==F){lineplot(plsnypos,links2,col=1)}}
      next3d()
      text3d(0,0,0, "Block2 positive")
      next3d()
      if(commonref==T){plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)}else{
        plot3d(themaxy*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      }
      if(heatmap==T){
        shade3d(myhyneg$obm$colMesh,alpha=alpha)
        if(is.null(links2)==F){lineplot(plsnyneg,links2,col=1)}
      }else{
        plot3d(plsnyneg,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links2)==F){lineplot(plsnyneg,links2,col=1)}}
      next3d()
      text3d(0,0,0, "Block2 negative")
    }
    
    if(dim(block1)[2]>2&dim(block2)[2]<3){
      open3d(windowRect=c(100,100,1000,1000))
      themaxx<-list(plsnxpos,plsnxneg)[which.max(css1)][[1]]
      themaxy<-list(plsnypos,plsnyneg)[which.max(css2)][[1]]
      mat <- matrix(1:8, ncol=2)
      layout3d(mat, height = rep(c(2,1), 2), model = "inherit")
      plot3d(themaxx*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      if(heatmap==T){
        shade3d(myhxpos$obm$colMesh,alpha=alpha)
        if(is.null(links1)==F){lineplot(plsnxpos,links1,col=1)}
      }else{
        plot3d(plsnxpos,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links1)==F){lineplot(plsnxpos,links1,col=1)}}
      next3d()
      text3d(0,0,0, "Block1 positive")
      next3d()
      plot3d(themaxx*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      if(heatmap==T){
        shade3d(myhxneg$obm$colMesh,alpha=alpha)
        if(is.null(links1)==F){lineplot(plsnxneg,links1,col=1)}
      }else{
        plot3d(plsnxneg,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links1)==F){lineplot(plsnxneg,links1,col=1)}}
      next3d()
      text3d(0,0,0, "Block1 negative")
      if(heatmap==T){
        png(file = "myplot.png", bg = "transparent",width = 780, height = 780)
        image.plot(xyz2img(cbind(myhypos$interpcoords,myhypos$pred),tolerance=myhypos$tol),asp=1,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),col=myhypos$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        par(new=T)
        plot(myhypos$mate,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="Block2 positive")
        if(is.null(links2)==F){lineplot(myhypos$mate,links2)}
        lines(myhypos$tpsgrid$grid$ngrid,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]))
        lines(myhypos$tpsgrid$grid$ngrid2,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]))
        par(new=T)
        plot(myhypos$mate2,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links2)==F){lineplot(myhypos$mate2,links2,col=2,lwd=2)}
        dev.off()
        show2d(filename="myplot.png",ignoreExtent = F)
        next3d()
        text3d(0,0,0, "Block2 positive")
        next3d()
        png(file = "myplot.png", bg = "transparent",width = 780, height = 780)
        image.plot(xyz2img(cbind(myhyneg$interpcoords,myhyneg$pred),tolerance=myhyneg$tol),asp=1,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),col=myhyneg$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        par(new=T)
        plot(myhyneg$mate,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="Block2 negative")
        if(is.null(links2)==F){lineplot(myhyneg$mate,links2)}
        lines(myhyneg$tpsgrid$grid$ngrid,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]))
        lines(myhyneg$tpsgrid$grid$ngrid2,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]))
        par(new=T)
        plot(myhyneg$mate2,xlim=range(themax[,1]),ylim=range(themaxy[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links2)==F){lineplot(myhyneg$mate2,links2,col=2,lwd=2)}
        dev.off()
        show2d(filename="myplot.png",ignoreExtent = F)
        next3d()
        text3d(0,0,0, "Block2 negative")
      }else{
        next3d()
        plot3d(cbind(themaxy,0)*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
        plot3d(cbind(plsnypos,0),bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links2)==F){lineplot(cbind(plsnypos,0),links2,col=1)}
        next3d()
        text3d(0,0,0, "Block2 positive")
        next3d()
        plot3d(cbind(themaxy,0)*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
        plot3d(cbind(plsnyneg,0),bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links2)==F){lineplot(cbind(plsnyneg,0),links2,col=1)}
        next3d()
        text3d(0,0,0, "Block2 negative")}
    } 
    
    if(dim(block1)[2]<3&dim(block2)[2]>2){
      open3d(windowRect=c(100,100,1000,1000)) 
      
      themaxx<-list(plsnxpos,plsnxneg)[which.max(css1)][[1]]
      themaxy<-list(plsnypos,plsnyneg)[which.max(css2)][[1]]
      
      mat <- matrix(1:8, ncol=2)
      layout3d(mat, height = rep(c(2,1), 2), model = "inherit")
      if(heatmap==T){
        png(file = "myplot.png", bg = "transparent",width = 780, height = 780)
        image.plot(xyz2img(cbind(myhxpos$interpcoords,myhxpos$pred),tolerance=myhxpos$tol),asp=1,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),col=myhypos$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        par(new=T)
        plot(myhxpos$mate,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="Block1 positive")
        if(is.null(links1)==F){lineplot(myhxpos$mate,links1)}
        lines(myhxpos$tpsgrid$grid$ngrid,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]))
        lines(myhxpos$tpsgrid$grid$ngrid2,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]))
        par(new=T)
        plot(myhxpos$mate2,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links1)==F){lineplot(myhxpos$mate2,links1,col=2,lwd=2)}
        dev.off()
        show2d(filename="myplot.png",ignoreExtent = F)
        next3d()
        text3d(0,0,0, "Block1 positive")
        next3d()
        png(file = "myplot.png", bg = "transparent",width = 780, height = 780)
        image.plot(xyz2img(cbind(myhxneg$interpcoords,myhxneg$pred),tolerance=myhxneg$tol),asp=1,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),col=myhxneg$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        par(new=T)
        plot(myhxneg$mate,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="Block1 negative")
        if(is.null(links1)==F){lineplot(myhxneg$mate,links1)}
        lines(myhxneg$tpsgrid$grid$ngrid,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]))
        lines(myhxneg$tpsgrid$grid$ngrid2,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]))
        par(new=T)
        plot(myhxneg$mate2,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        if(is.null(links1)==F){lineplot(myhxneg$mate2,links1,col=2,lwd=2)}
        dev.off()
        show2d(filename="myplot.png",ignoreExtent = F)
        next3d()
        text3d(0,0,0, "Block1 negative")
      }else{
        plot3d(cbind(themaxx,0)*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
        plot3d(cbind(plsnxpos,0),bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links1)==F){lineplot(cbind(plsnxpos,0),links1,col=1)}
        next3d()
        text3d(0,0,0, "Block1 positive")
        next3d()
        plot3d(cbind(themaxx,0)*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
        plot3d(cbind(plsnxneg,0),bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links2)==F){lineplot(cbind(plsnxneg,0),links1,col=1)}
        next3d()
        text3d(0,0,0, "Block1 negative")}
      next3d()
      plot3d(themaxy*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      if(heatmap==T){
        shade3d(myhypos$obm$colMesh,alpha=alpha)
        if(is.null(links2)==F){lineplot(plsnypos,links2,col=1)}
      }else{
        plot3d(plsnypos,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links2)==F){lineplot(plsnypos,links2,col=1)}}
      next3d()
      text3d(0,0,0, "Block2 positive")
      next3d()
      plot3d(themaxy*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      if(heatmap==T){
        shade3d(myhyneg$obm$colMesh,alpha=alpha)
        if(is.null(links2)==F){lineplot(plsnyneg,links2,col=1)}
      }else{
        plot3d(plsnyneg,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
        if(is.null(links2)==F){lineplot(plsnyneg,links2,col=1)}}
      next3d()
      text3d(0,0,0, "Block2 negative")
    }#### fine condizione 2d/3d
    
    if(dim(block1)[2]<3&dim(block2)[2]<3){
      
      themax<-list(plsnxpos,plsnxneg,plsnypos,plsnyneg)[which.max(c(css1,css2))][[1]]
      if(which.max(css1)<2){themaxx<-plsnxpos}else{themaxx<-plsnxneg}
      if(which.max(css2)<2){themaxy<-plsnypos}else{themaxy<-plsnyneg}
      
      par(mfrow=c(2,2))
      if(heatmap==F){
        if(commonref==T){plotmyarrays(plsnxpos,links=links1,txt=F,pch=19,xlim=range(themax[,1]),ylim=range(themax[,2]))}else{
          plotmyarrays(plsnxpos,links=links1,txt=F,pch=19,xlim=range(themaxx[,1]),ylim=range(themaxx[,2])) 
        }
        
        title("Block1 positive")
        
        if(commonref==T){plotmyarrays(plsnxneg,links=links1,txt=F,pch=19,xlim=range(themax[,1]),ylim=range(themax[,2]))}else{
          plotmyarrays(plsnxneg,links=links1,txt=F,pch=19,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]))
        }
        title("Block1 negative")
        
        if(commonref==T){plotmyarrays(plsnypos,links=links1,txt=F,pch=19,xlim=range(themax[,1]),ylim=range(themax[,2]))}else{
          plotmyarrays(plsnypos,links=links1,txt=F,pch=19,xlim=range(themaxy[,1]),ylim=range(themaxy[,2])) 
        }
        title("Block2 positive")
        
        if(commonref==T){plotmyarrays(plsnyneg,links=links1,txt=F,pch=19,xlim=range(themax[,1]),ylim=range(themax[,2]))}else{
          plotmyarrays(plsnyneg,links=links1,txt=F,pch=19,xlim=range(themaxy[,1]),ylim=range(themaxy[,2])) 
        }
        title("Block2 positive")}else{
          
          image.plot(xyz2img(cbind(myhxpos$interpcoords,myhxpos$pred),tolerance=myhxpos$tol),asp=1,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),col=myhxpos$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          par(new=T)
          plot(myhxpos$mate,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="Block1 positive")
          if(is.null(links1)==F){lineplot(myhxpos$mate,links1)}
          lines(myhxpos$tpsgrid$grid$ngrid,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]))
          lines(myhxpos$tpsgrid$grid$ngrid2,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]))
          par(new=T)
          plot(myhxpos$mate2,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(is.null(links1)==F){lineplot(myhxpos$mate2,links1,col=2,lwd=2)}
          image.plot(xyz2img(cbind(myhxneg$interpcoords,myhxneg$pred),tolerance=myhxneg$tol),asp=1,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),col=myhxneg$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          par(new=T)
          plot(myhxneg$mate,xlim=range(themaxx[,1]),ylim=range(themaxx[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="Block1 negative")
          if(is.null(links1)==F){lineplot(myhxneg$mate,links1)}
          lines(myhxneg$tpsgrid$grid$ngrid,xlim=range(themax[,1]),ylim=range(themax[,2]))
          lines(myhxneg$tpsgrid$grid$ngrid2,xlim=range(themax[,1]),ylim=range(themax[,2]))
          par(new=T)
          plot(myhxneg$mate2,xlim=range(themax[,1]),ylim=range(themax[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(is.null(links1)==F){lineplot(myhxneg$mate2,links1,col=2,lwd=2)}
          image.plot(xyz2img(cbind(myhypos$interpcoords,myhypos$pred),tolerance=myhypos$tol),asp=1,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),col=myhypos$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          par(new=T)
          plot(myhypos$mate,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="Block2 positive")
          if(is.null(links2)==F){lineplot(myhypos$mate,links2)}
          lines(myhypos$tpsgrid$grid$ngrid,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]))
          lines(myhypos$tpsgrid$grid$ngrid2,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]))
          par(new=T)
          plot(myhypos$mate2,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(is.null(links2)==F){lineplot(myhypos$mate2,links2,col=2,lwd=2)}
          image.plot(xyz2img(cbind(myhyneg$interpcoords,myhyneg$pred),tolerance=myhyneg$tol),asp=1,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),col=myhyneg$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          par(new=T)
          plot(myhyneg$mate,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",main="Block2 negative")
          if(is.null(links2)==F){lineplot(myhyneg$mate,links2)}
          lines(myhyneg$tpsgrid$grid$ngrid,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]))
          lines(myhyneg$tpsgrid$grid$ngrid2,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]))
          par(new=T)
          plot(myhyneg$mate2,xlim=range(themaxy[,1]),ylim=range(themaxy[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(is.null(links2)==F){lineplot(myhyneg$mate2,links2,col=2,lwd=2)}
          
        }
    } ##fine 2d/2d 
  }###   fine condizione forma-forma
  out<-list(plsob=thepls,allxscores=XScores,allyscores=YScores,if(!is.null(group)){xscoresmeans=xscoresmeans}else{ck=2},if(!is.null(group)){yscoresmeans=yscoresmeans}else{ck=2})
  if(length(out)>3){names(out)=c("pls","allxscores","allyscores","xscoresmeans","yscoresmeans")}else{NULL}
  return(out)
}
