#'ontrajplot
#' This function plots in 2D or 3D, with or without heatmap.
#' @param objptau	list: an object resulting from ptau6() function
#' @param z=5	numeric: the number of predictions evaluated at equally spaced size values on per-group ranges. 
#' @param triang=NULL	list: trinagulation structure for 3D data
#' @param links=NULL list: links structure
#' @param heat=F logical: if TRUE the heatmap is estimated. Necessarily slow. 
#' @param grid2d=T logical: if TRUE the TPS grid is calculated and displayed
#' @param linkss=T logical: if TRUE links on the source configuration, i.e. the most juvenile "starting shape, are plotted (only for 2D data). 
#' @param linksscol=2	numeric: links color for the source configuration.
#' @param lwdt=1 numeric: links width for target. 
#' @param lwds=1 numeric: links width for target.
#' @param scaleramp logical: "scaleramp" argument in meshDist() function in "Morpho" package
#' @param heatcolors character: Color palette for heatmap (by default=c("blue4","cyan2","yellow","red4"))
#' @param polyn numeric: the degree of regression. This should equal the same argument in ptau6(). 
#' @param mar numeric: "mar" argument in par() (default=c(0.3,0.3,0.3,0.3))
#' @param mai numeric: "mai" argument in par() (default=c(0,0,0.3,0))
#' @param oma numeric: "oma" argument in par() (default=c(0,0,3,0))
#' @param mag numeric: magnification parameter
#' @param levelcex2d numeric: text size for group names for 2D data
#' @return factorord the factor ordered by ptau() and outputted by that function.	
#' @return z the number of per-group predicitons to be saved. 
#' @return m the number of dimensions
#' @return links links if specified
#' @return triang triangulation if specified
#' @return linkss Logical. If TRUE links on source were specified. 				
#' @return heat Logical. If TRUE heatmap has been estimated. 	
#' @return grid2d Logical. If TRUE deformation grids in 2D have been estimated.
#' @return traspred	Predictions transported according to Linear shift Parallel Transport strategy and evaluated at equally spaced size values within per-group size ranges. 
#' @return origpred	Classical MANCOVA predictions evaluated at equally spaced size values within per-group size ranges.
#' @return serfac A factor with z individuals per each level in factorord
#' @return themaxt The largest configuration among all transported predictions. It is used for setting axes limits if analyses are done in Size and Shape Space. 
#' @return themaxo The largest configuration among all MANCOVA predictions. It is used for setting axes limits if analyses are done in Size and Shape Space.
#' @return traspred2 The same as traspred but with z+1 individuals for each level in factorord
#' @return origpred2 The same as origpred but with z+1 individuals for each level in factorord
#' @return serfac2 The same as serfac but with z+1 individuals for each level in factorord
#' @return alphas	A vector with 0 or 1 indicating the alpha parameter for each subscene in the main plots with z entries for each level in factorord
#' @return alphas2 The same as alpha with z+1 entries for each level in factorord
#' @return allobsp The log(determinant) of the jacobian matrix evaluated within the body of interest for transported predictions. This is a list with z entries for each level in factorord. Details for calculation can be found in Marquez et al (2012).
#' @return allithsp	A list with meshes colored with heatmap (in 3D) or images (in2D) for transported predictions.
#' @return allobso	The log(determinant) of the jacobian matrix evaluated within the body of interest for original per-group predictions (not transported). 
#' @return allithso	A list with meshes colored with heatmap (in 3D) or images (in2D) for original predictions.
#' @return allobsp2	The same as allobobsp but with z+1 entries for each level in factorord.
#' @return allithsp2 The same as allithsp but with z+1 entries for each level in factorord.
#' @return allobso2	The same as allobobso but with z+1 entries for each level in factorord.
#' @return allithso2 The same as allithso but with z+1 entries for each level in factorord.
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' ####  2D example
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
ontrajplot<-function(objptau,z=5,triang=NULL,links=NULL,heat=F,grid2d=T,linkss=T,linksscol=2,lwdt=1,lwds=1,scaleramp=F,heatcolors=c("blue4","cyan2","yellow","red4"),polyn=1,mar=c(0.3,0.3,0.3,0.3),mai=c(0,0,0.3,0),oma=c(0,0,3,0),mag=1,levelcex2d=1,levelcex3d=1,exts=NULL){
  if(is.null(triang)==F){if(ncol(triang)>3){stop("I need triangulation as nx3 matrix")}}
  
  
  traspred<-NULL
  origpred<-NULL
  for(i in 1:nlevels(objptau$factorord)){
    traspredi<-serpred(objptau$predictpure2[as.numeric(objptau$factorord)==i,],objptau$indepepure[as.numeric(objptau$factorord)==i],polyn=polyn,length.out=z)
    traspred<-rbind(traspred,traspredi)
    origpredi<-serpred(array2mat(objptau$predictpure[,,as.numeric(objptau$factorord)==i],table(objptau$factorord)[i],objptau$k*objptau$m),objptau$indepepure[as.numeric(objptau$factorord)==i],polyn=polyn,length.out=z)
    origpred<-rbind(origpred,origpredi)
  }
  traspred<-read.inn(traspred,objptau$k,objptau$m)
  origpred<-read.inn(origpred,objptau$k,objptau$m)
  serfac<-factor(rep(levels(objptau$factorord),each=z))
  serfac<-factor(serfac,levels=unique(serfac))
  if(is.null(exts)==T){exts<-rep(0.5,length(serfac))}else{exts=exts}
  
  
  
  traspred2<-NULL
  origpred2<-NULL
  for(i in 1:nlevels(serfac)){
    traspred2i<-abind::abind(traspred[,,as.numeric(serfac)==i][,,1,drop=F],traspred[,,as.numeric(serfac)==i])
    traspred2<-abind::abind(traspred2,traspred2i)
    origpred2i<-abind::abind(origpred[,,as.numeric(serfac)==i][,,1,drop=F],origpred[,,as.numeric(serfac)==i])
    origpred2<-abind::abind(origpred2,origpred2i)
  }
  
  serfac2<-factor(rep(levels(objptau$factorord),each=(z+1)))
  serfac2<-factor(serfac2,levels=unique(serfac2))
  alphas<-rep(1,length(serfac2))
  alphas[firstsfac(serfac2)]<-0
  alphas2<-rep(0,length(serfac2))
  alphas2[firstsfac(serfac2)]<-1
  if(is.null(exts)==T){exts2<-rep(0.5,length(serfac2))}else{
    exts2<-NULL
    for(i in 1:nlevels(serfac)){
      exts2i<-c(exts[as.numeric(serfac)==i][1],exts[as.numeric(serfac)==i])
      exts2<-c(exts2,exts2i)
    }
  }
  
  
  
  if(objptau$m>2&heat==T){
    
    if(is.null(triang)==T){stop("Heatmap in 3D requires triangulation")}
    allobsp<-NULL
    allithsp<-NULL
    for(i in 1:length(serfac)){
      ithp<-heat3d(objptau$CR,traspred[,,i],triang,graphics=F,colors=heatcolors,mag=mag)
      allithsp<-c(allithsp,list(ithp))
      allobspi<-ithp$obs3
      allobsp<-c(allobsp,list(allobspi))
      print(i)
    }
    
    allshapest<-abind::abind(list2array(subListExtract(allithsp,"mate2")),list2array(subListExtract(allithsp,"mate")))
    themaxt<-allshapest[,,which.max(apply(allshapest,3,cSize))]
    
    
    
    allobso<-NULL
    allithso<-NULL
    for(i in 1:nlevels(serfac)){
      for(j in 1:table(serfac)[i]){
        itho<-heat3d(origpred[,,as.numeric(serfac)==i][,,1],origpred[,,as.numeric(serfac)==i][,,j],triang,graphics=F,colors=heatcolors,mag=mag)
        allithso<-c(allithso,list(itho))
        allobsoi<-itho$obs3
        allobso<-c(allobso,list(allobsoi))
        print(paste(i,j,sep="_"))
      }
    }
    
    allshapeso<-abind::abind(list2array(subListExtract(allithso,"mate2")),list2array(subListExtract(allithso,"mate")))
    themaxo<-allshapeso[,,which.max(apply(allshapeso,3,cSize))]
    
    
    allobsp2<-NULL
    allithsp2<-NULL
    allobso2<-NULL
    allithso2<-NULL
    
    for(i in 1:nlevels(serfac)){
      allithsp2i<-c(allithsp[as.numeric(serfac)==i][1],allithsp[as.numeric(serfac)==i])
      allithsp2<-c(allithsp2,allithsp2i)
      allobsp2i<-c(allobsp[as.numeric(serfac)==i][1],allobsp[as.numeric(serfac)==i])
      allobsp2<-c(allobsp2,allobsp2i)
      
      allithso2i<-c(allithso[as.numeric(serfac)==i][1],allithso[as.numeric(serfac)==i])
      allithso2<-c(allithso2,allithso2i)
      allobso2i<-c(allobso[as.numeric(serfac)==i][1],allobso[as.numeric(serfac)==i])
      allobso2<-c(allobso2,allobso2i)
    }
    
    
    open3d(windowRect=c(100,100,1000,1000)) 
    mfrow3d(nlevels(objptau$factorord),z+1,sharedMouse = T,byrow=T)
    for(i in 1:length(serfac2)){
      plot3d(themaxt*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      meshDist(allithsp2[[i]]$obm$colMesh,distvec=allobsp2[[i]],from=min(unlist(allobsp2)),to=max(unlist(allobsp2)),scaleramp=scaleramp,add=T,alpha=alphas[i])
      text3d(0,0,0,rep(levels(objptau$factorord),each=(z+1))[i],alpha=alphas2[i],add=T,cex=levelcex3d)
      if(i<length(serfac2)){next3d()}
    }
    open3d(windowRect=c(100,100,1000,1000)) 
    mfrow3d(nlevels(objptau$factorord),z+1,sharedMouse = T,byrow=T)
    for(i in 1:length(serfac2)){
      plot3d(themaxo*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      meshDist(allithso2[[i]]$obm$colMesh,distvec=allobso2[[i]],from=min(unlist(allobso2)),to=max(unlist(allobso2)),scaleramp=scaleramp,add=T,alpha=alphas[i])
      text3d(0,0,0,rep(levels(objptau$factorord),each=(z+1))[i],alpha=alphas2[i],add=T,cex=levelcex3d)
      if(i<length(serfac2)){next3d()}
    }
  }
  
  if(objptau$m>2&heat==F){
    themaxt<-traspred[,,which.max(apply(traspred,3,cSize))]
    themaxo<-origpred[,,which.max(apply(origpred,3,cSize))]
    
    
    open3d(windowRect=c(100,100,1000,1000)) 
    mfrow3d(nlevels(objptau$factorord),z+1,sharedMouse = T,byrow=T)
    for(i in 1:length(serfac2)){
      plot3d(themaxt*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      plot3d(traspred2[,,i],bbox=F,type="s",asp=F,alpha=alphas[i],axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
      if(is.null(links)==F){if(i%in%c(1:length(serfac2))[firstsfac(serfac2)]){plot3d(themaxt*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F,add=T)}else{lineplot(traspred2[,,i],links)}}
      if(is.null(triang)==F){shade3d(plotsurf(traspred2[,,i],t(triang),plot=F),add=T,alpha=alphas[i]*0.5,col=2)}
      text3d(0,0,0,rep(levels(objptau$factorord),each=(z+1))[i],alpha=alphas2[i],add=T,cex=levelcex3d)
      if(i<length(serfac2)){next3d()}
    }
    open3d(windowRect=c(100,100,1000,1000)) 
    mfrow3d(nlevels(objptau$factorord),z+1,sharedMouse = T,byrow=T)
    for(i in 1:length(serfac2)){
      plot3d(themaxo*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
      plot3d(origpred2[,,i],bbox=F,type="s",asp=F,alpha=alphas[i],axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
      if(is.null(links)==F){if(i%in%c(1:length(serfac2))[firstsfac(serfac2)]){plot3d(themaxo*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F,add=T)}else{lineplot(origpred2[,,i],links)}}
      if(is.null(triang)==F){shade3d(plotsurf(origpred2[,,i],t(triang),plot=F),add=T,alpha=alphas[i]*0.5,col=2)}
      text3d(0,0,0,rep(levels(objptau$factorord),each=(z+1))[i],alpha=alphas2[i],add=T,cex=levelcex3d)
      if(i<length(serfac2)){next3d()}
    }
  }
  
  
  if(objptau$m<3){
    if(heat==F&grid2d==F){
      themaxt<-traspred[,,which.max(apply(traspred,3,cSize))]
      themaxo<-origpred[,,which.max(apply(origpred,3,cSize))]
      
      par(mfrow=c(nlevels(objptau$factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:length(serfac2)){
        plot(traspred2[,,i],xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),col=alphas[i],pch=19,cex=0.5,bty="n",xaxt="n",yaxt="n",xlab="",ylab="",asp=1)
        if(!is.null(links)){lineplot(traspred2[,,i],links,col=alphas[i],lwd=lwdt)}
        if(linkss==T){lineplot(traspred2[,,1],links,lwd=lwds,col=makeTransparent(linksscol,alpha=alphas[i]))}
        
        text(0,0,rep(levels(objptau$factorord),each=(z+1))[i],col=makeTransparent(1,alpha=alphas2[i]),cex=levelcex2d)
      }
      title(main="Ontogenetic trajectories of transported predictions at equally spaced size values within per-group ranges",outer=T)
      
    }
    if(heat==F&grid2d==T){  
      themaxt<-traspred[,,which.max(apply(traspred,3,cSize))]*mag
      themaxo<-origpred[,,which.max(apply(origpred,3,cSize))]*mag
      
      par(mfrow=c(nlevels(objptau$factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:length(serfac2)){
        tpsgridpaolo(traspred2[,,1],traspred2[,,i],xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),collandsTT=alphas[i],collandsYY=alphas[i],pch=19,cex=0.5,displ=F,axes2d=F,mag=mag,colgrid=makeTransparent(1,alpha=alphas[i]),ext=exts2[i])
        if(!is.null(links)){lineplot(traspred2[,,i],links,col=alphas[i])}
        if(linkss==T){lineplot(traspred2[,,1],links,lwd=lwds,col=makeTransparent(linksscol,alpha=alphas[i]))}
        text(0,0,rep(levels(objptau$factorord),each=(z+1))[i],col=makeTransparent(1,alpha=alphas2[i]),cex=levelcex2d)
      }
      title(main="Ontogenetic trajectories of transported predictions at equally spaced size values within per-group ranges",outer=T)
      
    }
    
    if(heat==F&grid2d==F){
      themaxt<-traspred[,,which.max(apply(traspred,3,cSize))]
      themaxo<-origpred[,,which.max(apply(origpred,3,cSize))]
      
      par(mfrow=c(nlevels(objptau$factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:nlevels(serfac2)){
        for(j in 1:table(serfac2)[i]){
          plot(origpred2[,,as.numeric(serfac2)==i][,,j],xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),col=alphas[as.numeric(serfac2)==i][j],pch=19,cex=0.5,bty="n",xaxt="n",yaxt="n",xlab="",ylab="",asp=1)
          if(!is.null(links)){lineplot(origpred2[,,as.numeric(serfac2)==i][,,j],links,col=alphas[as.numeric(serfac2)==i][j],lwd=lwdt)}
          if(linkss==T){lineplot(origpred2[,,as.numeric(serfac2)==i][,,1],links,lwd=lwds,col=makeTransparent(linksscol,alpha=alphas[as.numeric(serfac2)==i][j]))}
          text(0,0,rep(levels(objptau$factorord),each=(z+1))[as.numeric(serfac2)==i][j],col=makeTransparent(1,alpha=alphas2[as.numeric(serfac2)==i][j]),cex=levelcex2d)
        }
      }
      title(main="Ontogenetic trajectories of original per-group predictions at equally spaced size values within per-group ranges",outer=T)
    }
    if(heat==F&grid2d==T){  
      themaxt<-traspred[,,which.max(apply(traspred,3,cSize))]*mag
      themaxo<-origpred[,,which.max(apply(origpred,3,cSize))]*mag
      
      par(mfrow=c(nlevels(objptau$factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:nlevels(serfac2)){
        for(j in 1:table(serfac2)[i]){
          tpsgridpaolo(origpred2[,,as.numeric(serfac2)==i][,,1],origpred2[,,as.numeric(serfac2)==i][,,j],xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),collandsTT=alphas[as.numeric(serfac2)==i][j],collandsYY=alphas[i],pch=19,cex=0.5,displ=F,axes2d=F,mag=mag,colgrid=makeTransparent(1,alpha=alphas[as.numeric(serfac2)==i][j]),ext=exts2[as.numeric(serfac2)==i][j])
          if(!is.null(links)){lineplot(origpred2[,,as.numeric(serfac2)==i][,,j],links,col=alphas[as.numeric(serfac2)==i][j])}
          if(linkss==T){lineplot(origpred2[,,as.numeric(serfac2)==i][,,1],links,lwd=lwds,col=makeTransparent(linksscol,alpha=alphas[as.numeric(serfac2)==i][j]))}
          text(0,0,rep(levels(objptau$factorord),each=(z+1))[as.numeric(serfac2)==i][j],col=makeTransparent(1,alpha=alphas2[as.numeric(serfac2)==i][j]),cex=levelcex2d)
        }
      }
      title(main="Ontogenetic trajectories of original per-group predictions at equally spaced size values within per-group ranges",outer=T)
      
    }
    
    if(heat==T){ 
      allobsp<-NULL
      allithsp<-NULL
      for(i in 1:length(serfac)){
        ithp<-heat2d(objptau$CR,traspred[,,i],linkss=links,graphics=F,zlim=NULL,colors=heatcolors,mag=mag,ext=exts[i])
        allithsp<-c(allithsp,list(ithp))
        allobspi<-ithp$obs
        allobsp<-c(allobsp,list(allobspi))
        print(i)
      }
      
      allobso<-NULL
      allithso<-NULL
      for(i in 1:nlevels(serfac)){
        for(j in 1:table(serfac)[i]){
          itho<-heat2d(origpred[,,as.numeric(serfac)==i][,,1],origpred[,,as.numeric(serfac)==i][,,j],linkss=links,graphics=F,zlim=NULL,colors=heatcolors,mag=mag,ext=exts[i])
          allithso<-c(allithso,list(itho))
          allobsoi<-itho$obs
          allobso<-c(allobso,list(allobsoi))
          print(paste(i,j,sep="_"))
        }
      }
      
      allshapest<-abind::abind(list2array(subListExtract(allithsp,"mate2")),list2array(subListExtract(allithsp,"mate")))
      themaxt<-allshapest[,,which.max(apply(allshapest,3,cSize))]
      
      allshapeso<-abind::abind(list2array(subListExtract(allithso,"mate2")),list2array(subListExtract(allithso,"mate")))
      themaxo<-allshapeso[,,which.max(apply(allshapeso,3,cSize))]
      
      allobsp2<-NULL
      allithsp2<-NULL
      allobso2<-NULL
      allithso2<-NULL
      
      for(i in 1:nlevels(serfac)){
        allithsp2i<-c(allithsp[as.numeric(serfac)==i][1],allithsp[as.numeric(serfac)==i])
        allithsp2<-c(allithsp2,allithsp2i)
        allobsp2i<-c(allobsp[as.numeric(serfac)==i][1],allobsp[as.numeric(serfac)==i])
        allobsp2<-c(allobsp2,allobsp2i)
        
        allithso2i<-c(allithso[as.numeric(serfac)==i][1],allithso[as.numeric(serfac)==i])
        allithso2<-c(allithso2,allithso2i)
        allobso2i<-c(allobso[as.numeric(serfac)==i][1],allobso[as.numeric(serfac)==i])
        allobso2<-c(allobso2,allobso2i)
      }
      
      
      
      par(mfrow=c(nlevels(objptau$factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:length(serfac2)){
        if(i%in%c(1:length(serfac2))[firstsfac(serfac2)]){plot(allithsp2[[i]]$mate,bty="n",cex=0,asp=1,xlab="",ylab="",xaxt="n",yaxt="n",xlim=range(themaxt[,1]),ylim=range(themaxt[,2]))}else{
          image.plot(xyz2img(cbind(allithsp2[[i]]$interpcoords,allithsp2[[i]]$pred),tolerance=allithsp2[[i]]$tol),zlim=range(na.omit(c(allithsp2[[i]]$pred,unlist(allobsp2)))),asp=1,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),col=allithsp2[[i]]$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")  
          par(new=T)
          plot(allithsp2[[i]]$mate,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(is.null(links)==F){lineplot(allithsp2[[i]]$mate,links)}
          lines(allithsp2[[i]]$tpsgrid$grid$ngrid,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]))
          lines(allithsp2[[i]]$tpsgrid$grid$ngrid2,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]))
          par(new=T)
          plot(allithsp2[[i]]$mate2,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(linkss==T){lineplot(allithsp2[[i]]$mate2,links,col=2,lwd=lwds)}
        }
        text(0,0,rep(levels(objptau$factorord),each=(z+1))[i],col=makeTransparent(1,alpha=alphas2[i]),cex=levelcex2d)
      }
      title(main="Ontogenetic trajectories of transported predictions at equally spaced size values within per-group ranges",outer=T)
      
      
      par(mfrow=c(nlevels(objptau$factorord),z+1))
      par(mar=mar,mai=mai,oma=oma)
      for(i in 1:length(serfac2)){
        if(i%in%c(1:length(serfac2))[firstsfac(serfac2)]){plot(allithso2[[i]]$mate,bty="n",cex=0,asp=1,xlab="",ylab="",xaxt="n",yaxt="n",xlim=range(themaxt[,1]),ylim=range(themaxt[,2]))}else{
          image.plot(xyz2img(cbind(allithso2[[i]]$interpcoords,allithso2[[i]]$pred),tolerance=allithso2[[i]]$tol),zlim=range(na.omit(c(allithso2[[i]]$pred,unlist(allobso2)))),asp=1,xlim=range(themaxo[,1]),ylim=range(themaxo[,2]),col=allithso2[[i]]$cols,xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")  
          par(new=T)
          plot(allithso2[[i]]$mate,xlim=range(themaxt[,1]),ylim=range(themaxt[,2]),asp=1,pch=19,cex=1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(is.null(links)==F){lineplot(allithso2[[i]]$mate,links)}
          lines(allithso2[[i]]$tpsgrid$grid$ngrid,xlim=range(themaxo[,1]),ylim=range(themaxo[,2]))
          lines(allithso2[[i]]$tpsgrid$grid$ngrid2,xlim=range(themaxo[,1]),ylim=range(themaxo[,2]))
          par(new=T)
          plot(allithso2[[i]]$mate2,xlim=range(themaxo[,1]),ylim=range(themaxo[,2]),asp=1,pch=19,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
          if(linkss==T){lineplot(allithso2[[i]]$mate2,links,col=2,lwd=lwds)}
        }
        text(0,0,rep(levels(objptau$factorord),each=(z+1))[i],col=makeTransparent(1,alpha=alphas2[i]),cex=levelcex2d)
      }
      title(main="Ontogenetic trajectories of original per-group predictions at equally spaced size values within per-group ranges",outer=T)
    }
  }
  if(heat==F){
    allobsp<-NULL
    allobsp2<-NULL
    allithsp<-NULL
    allithsp2<-NULL
    allobso<-NULL
    allobso2<-NULL
    allithso<-NULL
    allithso2<-NULL
  }
  out<-list(factorord=objptau$factorord,z=z,m=objptau$m,links=links,triang=triang,linkss=linkss,heat=heat,grid2d=grid2d,traspred=traspred, origpred=origpred, serfac=serfac, themaxt=themaxt,themaxo=themaxo,traspred2=traspred2,origpred2=origpred2,serfac2=serfac2,alphas=alphas,alphas2=alphas2,allobsp=allobsp,allithsp=allithsp,allobso=allobso,allithso=allithso,allobsp2=allobsp2,allithsp2=allithsp2,allobso2=allobso2,allithso2=allithso2,exts=exts,exts2=exts2)
  out
}

