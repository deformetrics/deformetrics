#'replotontraj
#' This function takes an object from finction ontrajplot() and uses it in order to re-plot rapidly the results of ontrajplot().
#' @param ontrajplotobject list: an object outputted by ontrajplot()
#' @param zlimp numeric: range for colors for heatmap in 2D for transported predictions.
#' @param zlimo numeric: range for colors for heatmap in 2D for original predictions.
#' @param fromp numric: lower limit for colors for heatmap in 2D for transported predictions.
#' @param top numeric: higher limit for colors for heatmap in 2D for transported predictions.
#' @param fromo numeric: lower limit for colors for heatmap in 2D for original predictions.
#' @param too numeric: higher limit for colors for heatmap in 2D for transported predictions.
#' @param linksscol numeric: color of links plotted on source. 
#' @param lwdt numeric: width of links for target.
#' @param lwds numeric: width of links for source. 
#' @param scaleramp logical: "scaleramp" argument in meshDist() function in "Morpho" package
#' @param heatcolors character: Color palette for heatmap (by default=c("blue4","cyan2","yellow","red4"))
#' @param polyn numeric: the degree of regression. This should equal the same argument in ptau6(). 
#' @param mar numeric: "mar" argument in par() (default=c(0.3,0.3,0.3,0.3))
#' @param mai numeric: "mai" argument in par() (default=c(0,0,0.3,0))
#' @param oma numeric: "oma" argument in par() (default=c(0,0,3,0))
#' @param mag numeric: magnification parameter
#' @param levelcex2d numeric: text size for group names for 2D data
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
####  2D example
#' data(group)
#' data(my2d)
#' mysel<-group%in%c("Macroscelides_proboscideus","Petrodromus_tetradactylus","Elephantulus_rozeti","Elephantulus_edwardii")
#' linksdors<-list(c(1,2),c(37,7),c(12,4),c(27,28),c(25,21),c(38,40),c(9,10),c(2,3),c(3,4),c(1,7),c(1,6),c(3,5),c(6,40),c(5,9),c(40,8),c(8,9),c(1,7),c(7,6),c(3,4),c(4,5),c(39,38),c(38,35),c(35,37),c(37,39),c(35,34),c(34,33),c(33,32),c(32,31),c(31,30),c(30,29),c(29,37),c(37,36),c(36,29),c(28,31),c(28,30),c(13,10),c(10,11),c(11,12),c(12,13),c(13,14),c(14,16),c(16,17),c(17,20),c(20,19),c(19,18),c(18,12),c(18,15),c(15,12),c(21,19),c(21,20),c(24,25),c(25,26),c(26,27),c(27,24),c(26,24),c(24,23),c(23,22),c(22,8),c(8,2))       
#' dors4<-my2d[,,mysel]
#' factordors4<-drop.levels(group[mysel],reorder=T)
#' factordors4<-factor(factordors4,levels=unique(factordors4))
#' adors4<-procSym(dors4,scale=F,pcAlign=F,reflect=F)
#' mypredictbook<-read.inn(predict(lm(array2mat(adors4$orpdata,105,80)~poly(adors4$size,1,raw=T)*factordors4)),40,2)
#' procbook<-procSym(mypredictbook,pcAlign=F)
#' # linear shift in the shape space
#' plot(procbook$PCscores[,1:2],pch=as.numeric(factordors4),asp=1)
#' objptau2d<-ptau6(dors4,factordors4,CSinit=T,perm=2)
#' ontra20<-ontrajplot(objptau2d,links=linksdors,grid2d=F)#### just shapes
#' ontra21<-ontrajplot(objptau2d,links=linksdors,heat=F,mag=2)#### tpsgrid
#' replotontraj(ontra21)##replot rapidly
#' ontra22<-ontrajplot(objptau,links=linksdors,heat=T,mag=2)#### tpsgrid and heatmap necessarily slow
#' replotontraj(ontra22,zlimo=c(-3,2)) ##replot rapidly by adjusting zlim for original per group predictions 
#' ########  3D example
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
#' group<-factor(substr(dimnames(amy3d$orpdata)[[3]],1,7))
#' objptau3<-ptau6(amy3d$orpdata,group,CSinit=T,perm=2)
#' ## not run
#' ontra31<-ontrajplot(objptau3, links=NULL)#### only shapes with links
#' ontra32<-ontrajplot(objptau3,triang=t(sur_ent$it)) ## only meshes no heatmap
#' replotontraj(ontra32,zlimo=c(-3,2))##replot rapidly
#' ontra33<-ontrajplot(objptau3,triang=t(sur_ent$it),heat=T,mag=1)#### mesh with heatmap in 3D; necessarily slow: about 20 minutes
#' replotontraj(ontra33)##replot rapidly the object
#' replotontraj(ontra33)##replot rapidly the object
#' }
#' @export 
replotontraj<-function(ontrajplotobject,zlimp=NULL,zlimo=NULL,fromp=NULL,top=NULL,fromo=NULL,too=NULL,linksscol=2,lwdt=1,lwds=1,scaleramp=F,heatcolors=c("blue4","cyan2","yellow","red4"),polyn=1,mar=c(0.3,0.3,0.3,0.3),mai=c(0,0,0.3,0),oma=c(0,0,3,0),mag=1,levelcex2d=1,levelcex3d=1){
  factorord<-ontrajplotobject$factorord
  links<-ontrajplotobject$links
  linkss<-ontrajplotobject$linkss
  heat<-ontrajplotobject$heat
  grid2d<-ontrajplotobject$grid2d
  triang<-ontrajplotobject$triang
  traspred<-ontrajplotobject$traspred
  origpred<-ontrajplotobject$origpred
  serfac<-ontrajplotobject$serfac
  themaxt<-ontrajplotobject$themaxt
  themaxo<-ontrajplotobject$themaxo
  traspred2<-ontrajplotobject$traspred2
  origpred2<-ontrajplotobject$origpred2
  serfac2<-ontrajplotobject$serfac2
  alphas<-ontrajplotobject$alphas
  alphas2<-ontrajplotobject$alphas2
  allobsp2<-ontrajplotobject$allobsp2
  allobso2<-ontrajplotobject$allobso2
  allithsp2<-ontrajplotobject$allithsp2
  allithso2<-ontrajplotobject$allithso2
  z<-ontrajplotobject$z
  m<-ontrajplotobject$m
  exts<-ontrajplotobject$exts
  ext2<-ontrajplotobject$exts2
  
if(m<3){
  if(is.null(zlimp)==T){zlimp=range(na.omit(c(unlist(subListExtract(allithsp2,"pred")),unlist(allobsp2))))}else{zlimp=zlimp}
  if(is.null(zlimo)==T){zlimo=range(na.omit(c(unlist(subListExtract(allithso2,"pred")),unlist(allobso2))))}else{zlimp=zlimo}
}

if(m>2){
  
  if(is.null(fromp)==T){fromp=min(na.omit(c(unlist(subListExtract(allithsp2,"obs")),unlist(allobsp2))))}else{fromp=fromp}
  if(is.null(fromo)==T){fromo=min(na.omit(c(unlist(subListExtract(allithso2,"obs")),unlist(allobso2))))}else{fromo=fromo}
  
  if(is.null(top)==T){top=max(na.omit(c(unlist(subListExtract(allithsp2,"obs")),unlist(allobsp2))))}else{top=top}
  if(is.null(too)==T){too=max(na.omit(c(unlist(subListExtract(allithso2,"obs")),unlist(allobso2))))}else{too=too}
  }
  
  
  
  
  if(m>2&heat==T){
    if(is.null(triang)==T){stop("Heatmap in 3D requires triangulation")}
    
    open3d(windowRect=c(100,100,1000,1000)) 
    mfrow3d(nlevels(factorord),z+1,sharedMouse = T,byrow=T)
    for(i in 1:length(serfac2)){
      plot3d(themaxt*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      meshDist(allithsp2[[i]]$obm$colMesh,distvec=allobsp2[[i]],from=fromp,to=top,scaleramp=scaleramp,add=T,alpha=alphas[i])
      text3d(0,0,0,rep(levels(factorord),each=(z+1))[i],alpha=alphas2[i],add=T,cex=levelcex3d)
      if(i<length(serfac2)){next3d()}
    }
    open3d(windowRect=c(100,100,1000,1000)) 
    mfrow3d(nlevels(factorord),z+1,sharedMouse = T,byrow=T)
    for(i in 1:length(serfac2)){
      plot3d(themaxo*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      meshDist(allithso2[[i]]$obm$colMesh,distvec=allobso2[[i]],from=fromo,to=too,scaleramp=scaleramp,add=T,alpha=alphas[i])
      text3d(0,0,0,rep(levels(factorord),each=(z+1))[i],alpha=alphas2[i],add=T,cex=levelcex3d)
      if(i<length(serfac2)){next3d()}
    }
  }
  
  if(m>2&heat==F){
    open3d(windowRect=c(100,100,1000,1000)) 
    mfrow3d(nlevels(factorord),z+1,sharedMouse = T,byrow=T)
    for(i in 1:length(serfac2)){
      plot3d(themaxt*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      plot3d(traspred2[,,i],bbox=F,type="s",asp=F,alpha=alphas[i],axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
      if(is.null(links)==F){if(i%in%c(1:length(serfac2))[firstsfac(serfac2)]){plot3d(themaxt*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F,add=T)}else{lineplot(traspred2[,,i],links)}}
      if(is.null(triang)==F){shade3d(plotsurf(traspred2[,,i],t(triang),plot=F),add=T,alpha=alphas[i]*0.5,col=2)}
      text3d(0,0,0,rep(levels(factorord),each=(z+1))[i],alpha=alphas2[i],add=T,cex=levelcex3d)
      if(i<length(serfac2)){next3d()}
    }
    open3d(windowRect=c(100,100,1000,1000)) 
    mfrow3d(nlevels(factorord),z+1,sharedMouse = T,byrow=T)
    for(i in 1:length(serfac2)){
      plot3d(themaxo*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      plot3d(origpred2[,,i],bbox=F,type="s",asp=F,alpha=alphas[i],axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
      if(is.null(links)==F){if(i%in%c(1:length(serfac2))[firstsfac(serfac2)]){plot3d(themaxo*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F,add=T)}else{lineplot(origpred2[,,i],links)}}
      if(is.null(triang)==F){shade3d(plotsurf(origpred2[,,i],t(triang),plot=F),add=T,alpha=alphas[i]*0.5,col=2)}
      text3d(0,0,0,rep(levels(factorord),each=(z+1))[i],alpha=alphas2[i],add=T,cex=levelcex3d)
      if(i<length(serfac2)){next3d()}
    }
  }
  
  
  if(m<3){
    if(heat==F&grid2d==F){
      par(mfrow=c(nlevels(factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:length(serfac2)){
        plot(traspred2[,,i],xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),col=alphas[i],pch=19,cex=0.5,bty="n",xaxt="n",yaxt="n",xlab="",ylab="",asp=1)
        if(!is.null(links)){lineplot(traspred2[,,i],links,col=alphas[i],lwd=lwdt)}
        if(linkss==T){lineplot(traspred2[,,1],links,lwd=lwds,col=makeTransparent(linksscol,alpha=alphas[i]))}
        
        text(0,0,rep(levels(factorord),each=(z+1))[i],col=makeTransparent(1,alpha=alphas2[i]),cex=levelcex2d)
      }
      title(main="Ontogenetic trajectories of transported predictions at equally spaced size values within per-group ranges",outer=T)
      
    }
    if(heat==F&grid2d==T){  
      par(mfrow=c(nlevels(factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:length(serfac2)){
        tpsgridpaolo(traspred2[,,1],traspred2[,,i],xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),collandsTT=alphas[i],collandsYY=alphas[i],pch=19,cex=0.5,displ=F,axes2d=F,mag=mag,colgrid=makeTransparent(1,alpha=alphas[i]),ext=exts2[i])
        if(!is.null(links)){lineplot(traspred2[,,i],links,col=alphas[i])}
        if(linkss==T){lineplot(traspred2[,,1],links,lwd=lwds,col=makeTransparent(linksscol,alpha=alphas[i]))}
        text(0,0,rep(levels(factorord),each=(z+1))[i],col=makeTransparent(1,alpha=alphas2[i]),cex=levelcex2d)
      }
      title(main="Ontogenetic trajectories of transported predictions at equally spaced size values within per-group ranges",outer=T)
      
    }
    
    if(heat==F&grid2d==F){
      par(mfrow=c(nlevels(factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:nlevels(serfac2)){
        for(j in 1:table(serfac2)[i]){
          plot(origpred2[,,as.numeric(serfac2)==i][,,j],xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),col=alphas[as.numeric(serfac2)==i][j],pch=19,cex=0.5,bty="n",xaxt="n",yaxt="n",xlab="",ylab="",asp=1)
          if(!is.null(links)){lineplot(origpred2[,,as.numeric(serfac2)==i][,,j],links,col=alphas[as.numeric(serfac2)==i][j],lwd=lwdt)}
          if(linkss==T){lineplot(origpred2[,,as.numeric(serfac2)==i][,,1],links,lwd=lwds,col=makeTransparent(linksscol,alpha=alphas[as.numeric(serfac2)==i][j]))}
          text(0,0,rep(levels(factorord),each=(z+1))[as.numeric(serfac2)==i][j],col=makeTransparent(1,alpha=alphas2[as.numeric(serfac2)==i][j]),cex=levelcex2d)
        }
      }
      title(main="Ontogenetic trajectories of original per-group predictions at equally spaced size values within per-group ranges",outer=T)
    }
    if(heat==F&grid2d==T){  
      par(mfrow=c(nlevels(factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:nlevels(serfac2)){
        for(j in 1:table(serfac2)[i]){
          tpsgridpaolo(origpred2[,,as.numeric(serfac2)==i][,,1],origpred2[,,as.numeric(serfac2)==i][,,j],xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),collandsTT=alphas[as.numeric(serfac2)==i][j],collandsYY=alphas[i],pch=19,cex=0.5,displ=F,axes2d=F,mag=mag,colgrid=makeTransparent(1,alpha=alphas[as.numeric(serfac2)==i][j]),ext=exts2[as.numeric(serfac2)==i][j])
          if(!is.null(links)){lineplot(origpred2[,,as.numeric(serfac2)==i][,,j],links,col=alphas[as.numeric(serfac2)==i][j])}
          if(linkss==T){lineplot(origpred2[,,as.numeric(serfac2)==i][,,1],links,lwd=lwds,col=makeTransparent(linksscol,alpha=alphas[as.numeric(serfac2)==i][j]))}
          text(0,0,rep(levels(factorord),each=(z+1))[as.numeric(serfac2)==i][j],col=makeTransparent(1,alpha=alphas2[as.numeric(serfac2)==i][j]),cex=levelcex2d)
        }
      }
      title(main="Ontogenetic trajectories of original per-group predictions at equally spaced size values within per-group ranges",outer=T)
      
    }
    
    if(heat==T){ 
      
      par(mfrow=c(nlevels(factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:length(serfac2)){
        if(i%in%c(1:length(serfac2))[firstsfac(serfac2)]){plot(allithsp2[[i]]$mate,bty="n",cex=0,asp=1,xlab="",ylab="",xaxt="n",yaxt="n",xlim=range(themaxt[,1]),ylim=range(themaxt[,2]))}else{
          image.plot(xyz2img(cbind(allithsp2[[i]]$interpcoords,allithsp2[[i]]$pred),tolerance=allithsp2[[i]]$tol),zlim=zlimp,asp=1,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),col=allithsp2[[i]]$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")  
          par(new=T)
          plot(allithsp2[[i]]$mate,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(is.null(links)==F){lineplot(allithsp2[[i]]$mate,links)}
          lines(allithsp2[[i]]$tpsgrid$grid$ngrid,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]))
          lines(allithsp2[[i]]$tpsgrid$grid$ngrid2,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]))
          par(new=T)
          plot(allithsp2[[i]]$mate2,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(linkss==T){lineplot(allithsp2[[i]]$mate2,links,col=2,lwd=lwds)}
        }
        text(0,0,rep(levels(factorord),each=(z+1))[i],col=makeTransparent(1,alpha=alphas2[i]),cex=levelcex2d)
      }
      title(main="Ontogenetic trajectories of transported predictions at equally spaced size values within per-group ranges",outer=T)
      
      
      par(mfrow=c(nlevels(factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:length(serfac2)){
        if(i%in%c(1:length(serfac2))[firstsfac(serfac2)]){plot(allithso2[[i]]$mate,bty="n",cex=0,asp=1,xlab="",ylab="",xaxt="n",yaxt="n",xlim=range(themaxt[,1]),ylim=range(themaxt[,2]))}else{
          image.plot(xyz2img(cbind(allithso2[[i]]$interpcoords,allithso2[[i]]$pred),tolerance=allithso2[[i]]$tol),zlim=zlimo,asp=1,xlim=range(themaxo[,1]),ylim=range(themaxo[,2]),col=allithso2[[i]]$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")  
          par(new=T)
          plot(allithso2[[i]]$mate,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(is.null(links)==F){lineplot(allithso2[[i]]$mate,links)}
          lines(allithso2[[i]]$tpsgrid$grid$ngrid,xlim=range(themaxo[,1]),ylim=range(themaxo[,2]))
          lines(allithso2[[i]]$tpsgrid$grid$ngrid2,xlim=range(themaxo[,1]),ylim=range(themaxo[,2]))
          par(new=T)
          plot(allithso2[[i]]$mate2,xlim=range(themaxo[,1]),ylim=range(themaxo[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(linkss==T){lineplot(allithso2[[i]]$mate2,links,col=2,lwd=lwds)}
        }
        text(0,0,rep(levels(factorord),each=(z+1))[i],col=makeTransparent(1,alpha=alphas2[i]),cex=levelcex2d)
      }
      title(main="Ontogenetic trajectories of original per-group predictions at equally spaced size values within per-group ranges",outer=T)
    }
  }
  
  
}
