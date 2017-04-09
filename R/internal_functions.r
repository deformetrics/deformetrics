#' Internal_functions
#'
#' Here are reported a collection of internal functions
#' @author Paolo Piras
#' @export  
ptau6<-function(array,factor,CSinit=T,sepure=F,polyn=1,CR=NULL,locs=NULL,perm=999){
  library(Morpho)
  library(shapes)
  library(vegan)
  warning("WARNING: this function reorders data (if they are not) in increasing size order within each level")
  k<-dim(array)[1]
  m<-dim(array)[2]
  n<-dim(array)[3]
  
  newarray<-array[,,order(factor,centroid.size(array))]
  
  factor<-factor[order(factor)]
  
  
  cbind(dimnames(newarray)[[3]],as.character(factor))
  
  depepure<-array2mat(procSym(newarray,pcAlign=F,CSinit=CSinit,scale=F)$orpdata,n,k*m)
  if(sepure==T){depepure<-array2mat(sepgpa(newarray,factor,CSinit=CSinit,scale=F)$mountedorpdata,n,k*m)}
  indepepure<-centroid.size(newarray)
  print("Individual multivariate (linear) regression between shape and size" )
  print(manymultgr(depepure,indepepure,factor,steps=perm))
  print("Pairwise *linear* mancova p-values on original data") 
  print(pwpermancova(depepure,indepepure,factor,nperm=perm)$p_adonis_pred1_pred2)
  thedatapure<-data.frame(indepure=indepepure,depure=depepure)
  thelmlistpure<-NULL
  mypredictpure<-NULL
  myresidpure<-NULL
  for(i in 1:nlevels(factor)){
    thelmpurei<-lm(as.matrix(thedatapure[,-1][as.numeric(factor)==i,])~poly(indepure[as.numeric(factor)==i],degree=polyn,raw=TRUE),data=thedatapure)
    mypredictpurei<-predict(thelmpurei)
    myresidpurei<-resid(thelmpurei)
    thelmlistpure<-c(thelmlistpure,list(thelmpurei))
    mypredictpure<-rbind(mypredictpure,mypredictpurei)
    myresidpure<-rbind(myresidpure,myresidpurei)
  }
  
  mypredictpure<-read.inn(mypredictpure,k,m)
  myresidpure<-read.inn(myresidpure,k,m)
  

if(is.null(CR)==T){CR<-procSym(mypredictpure[,,firstsfac(factor)],scale=F,pcAlign=F,reflect=F,CSinit=F)$mshape}else{CR<-CR}
if(is.null(locs)==T){locs<-mypredictpure[,,firstsfac(factor)]}else{locs<-locs}

  prls<-lshift2(mypredictpure,factor,CSinit=CSinit,CR=CR,locs=locs)
  
  common<-array2mat(procSym(newarray,pcAlign=,scale=F,CSinit=CSinit)$orpdata,n,k*m)
  print("Pairwise multivariate (linear) regression between shape and size" )
  print(pwpermancova(common,indepepure,factor,nperm=perm)$p_adonis_pred1_pred2)
  space1<-prcomp(array2mat(prls$transported,n,k*m))
  space1mshape<-procSym(prls$transported,pcAlign=F)$mshape
  origtrasp<-prls$transported+myresidpure
  
  origproj<-predict(space1,array2mat(origtrasp,n,k*m))
  print("Pairwise multivariate (linear) regression between shape of transported data and size" )
  print(pwpermancova(array2mat(origtrasp,n,k*m),indepepure,factor,nperm=perm)$p_adonis_pred1_pred2)
  depepure2<-origtrasp
  
  print("Individual multivariate (linear) regression between shape of transported data and size" )
  print(manymultgr(array2mat(depepure2,n,k*m),indepepure,factor,steps=perm))
  
  thedatapure2<-data.frame(indepure2=indepepure,depure2=array2mat(depepure2,n,k*m))
  thelmlistpure2<-NULL
  mypredictpure2<-NULL
  myresidpure2<-NULL
  for(i in 1:nlevels(factor)){
    thelmpure2i<-lm(as.matrix(thedatapure2[,-1][as.numeric(factor)==i,])~poly(indepure2[as.numeric(factor)==i],degree=polyn,raw=TRUE),data=thedatapure2)
    mypredictpure2i<-predict(thelmpure2i)
    myresidpure2i<-resid(thelmpure2i)
    thelmlistpure2<-c(thelmlistpure2,list(thelmpure2i))
    mypredictpure2<-rbind(mypredictpure2,mypredictpure2i)
    myresidpure2<-rbind(myresidpure2,myresidpure2i)
  }
  
  
  out<-list(k=k,m=m,n=n,arrayord=newarray,factorord=factor,CR=CR,locs=locs,depepure=depepure,indepepure=indepepure,thelmlistpure=thelmlistpure,predictpure=mypredictpure,residpure=myresidpure,shifted=array2mat(prls$transported,n,k*m),thelmlistpure2=thelmlistpure2,predictpure2=array2mat(prls$transported,n,k*m),predictpure3=mypredictpure2,residpure2=myresidpure2,space1=space1,origtrasp=array2mat(origtrasp,n,k*m),origproj=origproj,space1mshape=space1mshape)                                   
}
#' export


plsgen<-function(block1,block2,plsn=1,links1=NULL,links2=NULL,commonref=F,heatmap=F,heatcolors=c("blue4","cyan2","yellow","red4"),triang1=NULL,triang2=NULL,alpha=1,S1=NA,S2=NA,from1=NULL,to1=NULL,from2=NULL,to2=NULL,rounds=0,sdx=1,sdy=1,labels=NULL,group=NULL,col=1,pch=19,colgroup=NULL,zlim2d=NULL){
  require(Morpho)
  require(shapes)
  
  thepls<-pls2B(block1,block2,rounds=rounds)
  XScores<-thepls$Xscores
  YScores<-thepls$Yscores
  plot(XScores[, plsn],YScores[,plsn],asp=1,col=col,pch=pch)
  if(!is.null(labels)){textxy(XScores[, plsn],YScores[, plsn],labels)}
  
  if(!is.null(group)){plot2dhull(cbind(XScores[, plsn],YScores[, plsn]),group,1,pch=pch,col=col,colhull=colgroup,labels=labels)}else{NULL}
  
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
#' export
plotptau5<-function(objptau,links=NULL,projorig=T,whichaxes=c(1,2),cexscale=round(max(centroid.size(objptau$arrayord)),digits=0),shapescale=10,mag=1,subplotdim=2,shiftnegy=1,shiftposy=1,shiftnegx=1,shiftposx=1,col=as.numeric(objptau$factorord),pch=19,triang=NULL,from=0,topcs=NULL,tore=NULL,plotsource=T){
  library(TeachingDemos)
  k<-dim(objptau$arrayord)[1]
  m<-dim(objptau$arrayord)[2]
  n<-dim(objptau$arrayord)[3]
  transp<-read.inn(objptau$shifted,k,m)
  
  origsize<-objptau$indepepure
  
if(projorig==T){
ordiwithshapes(objptau$space1mshape,objptau$space1$x,objptau$space1$rotation,procSym=F,whichaxes=whichaxes,addata=objptau$origproj[,whichaxes],asp=1,factraj=objptau$factorord,cex=origsize/cexscale,triang=triang,from=from,to=topcs,links=links,mag=mag,shiftnegy=1.5,shiftposy=2,col=col,pch=pch,subplotdim=subplotdim,plotsource=plotsource)
  title("LS on predictions in PCA space and original data re-projected")}else{

ordiwithshapes(objptau$space1mshape,objptau$space1$x,objptau$space1$rotation,procSym=F,whichaxes=whichaxes,asp=1,factraj=objptau$factorord,cex=origsize/cexscale,triang=triang,from=from,to=topcs,links=links,mag=mag,shiftnegy=1.5,shiftposy=2,col=col,pch=pch,subplotdim=subplotdim,plotsource=plotsource)
  title("LS on predictions in PCA space")}

  
  
  ranges<-apply(rbind(objptau$space1$x,objptau$origproj),2,range)
  
  ratesorig<-ratesbygroup(read.inn(objptau$predictpure2,k,m),objptau$factorord,objptau$indepepure)
  plot(objptau$indepepure,ratesorig,pch=pch,col=col)
  for(i in 1:nlevels(objptau$factorord)){
    lines(objptau$indepepure[as.numeric(objptau$factorord)==i][-firstsfac(objptau$factorord)],ratesorig[as.numeric(objptau$factorord)==i][-firstsfac(objptau$factorord)],col=col[firstsfac(objptau$factorord)][i])
  }
  title("Rates of shape change among consecutive predictions per unit size in original data")
  
  ratestrasp<-ratesbygroup(objptau$predictpure,objptau$factorord,objptau$indepepure)
  
  plot(objptau$indepepure,ratestrasp,pch=pch,col=col)
  for(i in 1:nlevels(objptau$factorord)){
    
    lines(objptau$indepepure[as.numeric(objptau$factorord)==i][-firstsfac(objptau$factorord)],ratestrasp[as.numeric(objptau$factorord)==i][-firstsfac(objptau$factorord)],col=col[firstsfac(objptau$factorord)][i])
  }
  
  
  
  title("Rates of shape change among consecutive predictions per unit size in transported data")
  
  
  
  plot(ratesorig,ratestrasp,col=col,pch=pch)
  title("Original rates of shape change per unit size vs those of transported data")
  
  adults<-NULL
  adultsdaplot<-NULL
  for(i in 1:nlevels(objptau$factorord)){
    adulti<-array(showPC2(objptau$space1$x[lastsfac(objptau$factorord),1:3][i,],objptau$space1$rotation[,1:3],objptau$space1mshape),dim=c(k,m,1))
    adultdaploti<-array(showPC2(objptau$space1$x[lastsfac(objptau$factorord),1:3][i,]*mag,objptau$space1$rotation[,1:3],objptau$space1mshape),dim=c(k,m,1))
    
    adults<-abind::abind(adults,adulti)
    adultsdaplot<-abind::abind(adultsdaplot,adultdaploti)
  }
  
  if(m<3){adultsdaplot<-abind::abind(adultsdaplot,array(rep(0,k*n),dim=c(k,1,nlevels(objptau$factorord))),along=2)}
  
  init<-array(showPC2(objptau$space1$x[firstsfac(objptau$factorord),1:3][1,],objptau$space1$rotation[,1:3],objptau$space1mshape),dim=c(k,m,1))
  initdaplot<-array(showPC2(objptau$space1$x[firstsfac(objptau$factorord),1:3][1,]*mag,objptau$space1$rotation[,1:3],objptau$space1mshape),dim=c(k,m,1))
  
  if(m<3){initdaplot<-abind::abind(initdaplot,array(rep(0,k*n),dim=c(k,1,1)),along=2)}
  
  if(is.null(triang)==F){
  adultsdaplotvis<-NULL
  for(i in 1:nlevels(objptau$factorord)){
  adultsdaplotvisi<-meshDist(plotsurf(adultsdaplot[,,i],triang,plot=F),plotsurf(initdaplot[,,1],triang,plot=F),to=tore,from=from,plot=F)
  daagg<-objptau$space1$x[1:n,1:3][lastsfac(objptau$factorord),][i,]
  adultsdaplotvisi<-translate3d(scalemesh(adultsdaplotvisi$colMesh,1/shapescale,center="none"),daagg[1],daagg[2],daagg[3])
 adultsdaplotvis<-c(adultsdaplotvis,list(adultsdaplotvisi))
  }
  }
  
  
  
  
  open3d()
  for(i in 1:nlevels(objptau$factorord)){
    lines3d(objptau$space1$x[1:n,1:3][as.numeric(objptau$factorord)==i,],col=i,add=T)
  }
  
  if(is.null(triang)==T){
  babedaagg<-rep.row(objptau$space1$x[firstsfac(objptau$factorord),1:3][1,],k)
  for(i in 1:nlevels(objptau$factorord)){
    daplottarei<-adultsdaplot[,,i]
    daagg<-rep.row(objptau$space1$x[1:n,1:3][lastsfac(objptau$factorord),][i,],k)
    plot3D((daplottarei/shapescale)+daagg,col=i,add=T,size=1,bbox=F)
    lineplot((daplottarei/shapescale)+daagg,links,col=i)
  }
  plot3D((initdaplot[,,1]/shapescale)+babedaagg,bbox=F,add=T)
  lineplot((initdaplot[,,1]/shapescale)+babedaagg,links)
  
  }else{
    
    for(i in 1:nlevels(objptau$factorord)){
      shade3d(adultsdaplotvis[[i]],add=T)
    }
    babedaagg<-rep.row(objptau$space1$x[firstsfac(objptau$factorord),1:3][1,],k)
  shade3d(plotsurf((initdaplot[,,1]/shapescale)+babedaagg,triang,plot=F),add=T,alpha=0.5)
    }
  title3d("LS data in the PCA; first three PC scores")
   decorate3d()
  out<-list(ratesorig=ratesorig,ratestrasp=ratestrasp)
  
}
#' export

benfrombygroup<-function(locs,array,factor,doopa=T){
  resu<-NULL
  for(j in 1:nlevels(factor)){
    for(k in 1:table(factor)[j]){
    resui<-jbe(locs[,,j],array[,,as.numeric(factor)==j][,,k],doopa=doopa)
    resu<-append(resu,resui)
    print(paste(j,k,sep="_"))
    }
  }
  resu
}
#' @export

plotancova<-function(y,x,group=NULL,pch=NULL,col=1,confint=T,cex=0.5,legend=T,labels=NULL,plot=T,xlab=NULL,ylab=NULL,xlim=range(x),ylim=range(y)){
  
  
  if(is.null(group)==T){group<-factor(c(rep(1,length(y))))}else{pch=pch}
  
  
  if(is.null(pch)==T){pch=19}else{pch=pch}
  dati<-data.frame(group,as.numeric(group),col,pch)
  
  if(is.null(xlab)==T){xlab<-c("x")}else{xlab<-xlab}
  if(is.null(ylab)==T){ylab<-c("y")}else{ylab<-ylab}
  if(plot==T){plot(x,y,xlim=xlim,ylim=ylim,col=col,pch=pch,cex=cex,xlab=xlab,ylab=ylab)}
  
  
  
  if(!is.null(labels)){textxy(x,y,labels)}
  
  
  lmlist<-NULL
  for(i in 1:nlevels(group)){ 
    lmi<-lm(y[as.numeric(group)==i]~x[as.numeric(group)==i])  
    lmlist<-c(lmlist,list(lmi))
  }
  names(lmlist)<-levels(group)
  
  
  datapredinflist<-NULL
  for(i in 1:nlevels(group)){
    datapredinfi<-cbind(x[as.numeric(group)==i], predict(lmlist[[i]], interval="confidence")[,2])
    datapredinfi<-datapredinfi[order(datapredinfi[,2]),]
    datapredinflist<-c(datapredinflist,list(datapredinfi))
  }  
  
  names(datapredinflist)<-levels(group)
  
  
  datapredsuplist<-NULL
  for(i in 1:nlevels(group)){
    datapredsupi<-cbind(x[as.numeric(group)==i], predict(lmlist[[i]], interval="confidence")[,3])
    datapredsupi<-datapredsupi[order(datapredsupi[,2]),]
    datapredsuplist<-c(datapredsuplist,list(datapredsupi))
  }    
  
  names(datapredsuplist)<-levels(group)
  
  if(plot==T){
    for(i in 1: length(lmlist)){
      abline(lmlist[[i]],col=dati[as.numeric(group)==i,][1,3],ylim=ylim,xlim=xlim)
    }   }                
  
  if(plot==T){
    if(confint==T){
      for(i in 1: length(lmlist)){
        if(anova(lmlist[[i]])$Pr[1]<0.05){
          lines(datapredinflist[[i]],cex=0.1,lty = 'dashed',col=dati[as.numeric(group)==i,][1,3],ylim=ylim,xlim=xlim)
        }}
      
      for(i in 1: length(lmlist)){
        if(anova(lmlist[[i]])$Pr[1]<0.05){
          lines(datapredsuplist[[i]],cex=0.1,lty = 'dashed',col=dati[as.numeric(group)==i,][1,3],ylim=ylim,xlim=xlim)
        }}
    }}
  
  if(legend==T){
    x11()
    plot(x,y,col="white")
    legend(min(x),max(y),unique(group), cex=1, col=unique(col), pch=unique(pch),box.col="white")}
  
  summarylist<-NULL
  for(i in 1:length(lmlist)){
    summarylisti<-summary(lmlist[[i]])  
    summarylist<-c(summarylist,list(summarylisti))
  }
  
  names(summarylist)<-levels(group)
  sint_results<-NULL
  
  for(i in 1:length(summarylist)){
    inti<-summarylist[[i]]$coefficients[1,1]
    p_inti<-summarylist[[i]]$coefficients[1,4]
    betai<-summarylist[[i]]$coefficients[2,1]
    p_betai<-summarylist[[i]]$coefficients[2,4]
    r_sqi<-summarylist[[i]]$r.squared
    
    sinti<-c(inti,p_inti,betai,p_betai,r_sqi)
    
    sint_results<-rbind(sint_results,sinti)
  }
  rownames(sint_results)<-levels(group)
  colnames(sint_results)<-c("Intercept","p-value Intercept","Beta","p-value Beta","R squared")
  
  
  print(sint_results)
  df<-data.frame(y=y,x=x,group=group)
  model<-lm(y~x*group)
  library(phia)

  slopint<-testInteractions(model, pairwise="group", slope="x") #####  TESTA SLOPE 
  slopintmat<- matrix(NA, nlevels(spe), nlevels(spe))
  slopintmat[lower.tri(slopintmat) ] <-round( slopint[1:(nrow(slopint)-1),5], 5)
  slopintmat<-t(slopintmat)
  colnames(slopintmat)<-levels(spe)
  rownames(slopintmat)<-levels(spe)
  slopintmat[lower.tri(slopintmat)]<-slopintmat[upper.tri(slopintmat)]

  
  
  elevint<-testInteractions(lm(y~x+group)) 
  elevintmat<- matrix(NA, nlevels(group), nlevels(group))
  elevintmat[lower.tri(elevintmat) ] <-round( elevint[1:(nrow(elevint)-1),5], 5)
  elevintmat<-t(elevintmat)
  colnames(elevintmat)<-levels(group)
  rownames(elevintmat)<-levels(group)
  elevintmat[lower.tri(elevintmat)]<-elevintmat[upper.tri(elevintmat)]

  out<-list(datapredsuplist=datapredsuplist,datapredinflist=datapredinflist,summarylist=summarylist,sint_results=sint_results,slopepval=slopintmat,elevpval=elevintmat)
return(out)
}
#' @export
heat3d<-function(source,target,triang,iter=3,linkss=NULL,linkst=NULL,plotlands=F,legend=T,cols=1,colt=2,plotsource=T,plottarget=T,collinkss=1,lwds=2,collinkst=2,lwdt=2,cexs=0.5,cext=0.5,colors=c("blue4","cyan2","yellow","red4"),alpha=1,ngrid=0,mag=1,graphics=T,to=NULL,from=NULL,scaleramp=F,lines=T){
  mate<-source
  mate2<-target
  mate2<-mate+(mate2-mate)*mag
  library(fields)
  tes<-tessell3d(triang,mate,iter)
  class(tes) <- "mesh3d"
  matr<-mate
  M<-t(tes$vb)[,1:3][-c(1:nrow(matr)),]
  
  
  tes2<-tessell3d(triang,mate2,iter)
  class(tes2) <- "mesh3d"
  
  M2<-t(tes2$vb)[,1:3][-c(1:nrow(matr)),]
  
  tpsgrid<-tpsgridpaolo(mate,mate2,linksTT=linkss,linksYY=linkss,graphics=F,ngrid=22)
  # 
  # veclist<-NULL
  # for(j in 1:nrow(M)){
  #   vecj<-NULL
  #   for(i in 1:nrow(matr)){
  #     vec1ji<--2*sqrt(abs(M[j,1]-matr[i,1]))
  #     vec2ji<--2*sqrt(abs(M[j,2]-matr[i,2]))
  #     vec3ji<--2*sqrt(abs(M[j,3]-matr[i,3]))
  #     vecji<-c(vec1ji,vec2ji,vec3ji)
  #     vecj<-rbind(vecj,vecji)
  #   }
  #   veclist<-c(veclist,list(vecj))
  # }
  # 
  
  veclist<-NULL
  for(j in 1:nrow(M)){
    vecj<-NULL
    for(i in 1:nrow(matr)){
      vec1ji<-2*(M[j,1]-matr[i,1])+2*(M[j,1]-matr[i,1])*log((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2)
      vec2ji<-2*(M[j,2]-matr[i,2])+2*(M[j,2]-matr[i,2])*log((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2)
      vec3ji<-2*(M[j,3]-matr[i,3])+2*(M[j,3]-matr[i,3])*log((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2)
      
      vecji<-c(vec1ji,vec2ji,vec3ji)
      vecj<-rbind(vecj,vecji)
    }
    veclist<-c(veclist,list(vecj))
  }
  
  
  jac<-NULL
  for(i in 1: length(veclist)){
    jaci<-t(tpsgrid$B)+t(tpsgrid$W)%*%veclist[[i]]
    jac<-c(jac,list(jaci))
  }
  
  jac2<-list2array(jac)
  mean11<-mean(jac2[1,1,])
  mean12<-mean(jac2[1,2,])
  mean13<-mean(jac2[1,3,])
  mean21<-mean(jac2[2,1,])
  mean22<-mean(jac2[2,2,])
  mean23<-mean(jac2[2,3,])
  mean31<-mean(jac2[3,1,])
  mean32<-mean(jac2[3,2,])
  mean33<-mean(jac2[3,3,])
  detmean<-det(matrix(c(mean11,mean21,mean31,mean12,mean22,mean32,mean31,mean32,mean33),ncol=3))
  myj<-unlist(lapply(jac,det))
  obs<-log2(myj/detmean)
  fit<- Tps(M, obs) 
  obs2<-predict(fit,M2)
  obs3<-c(obs2)
  
  obs3[is.na(obs3)]<-mean(obs3,na.rm=T)
  
  obm<-meshDist(plotsurf(rbind(mate2,M2),tes2$it,plot=F),distv=obs3,add=T,rampcolors =colors,to=to,from=from,scaleramp=scaleramp,plot=F,alpha=alpha)
  if(graphics==T){
    if(plotlands==T){deformGrid3d(mate,mate2,ngrid=ngrid,lines=lines,col1=cols,col2=colt)}
    shade3d(obm$colMesh,alpha=alpha)
    if(!is.null(linkst)){lineplot(mate2,linkst,col=collinkst,lwd=lwds)}
    if(plotsource==T){
      shade3d(plotsurf(source,t(triang),plot=F),alpha=0.5)
      if(!is.null(linkss)){lineplot(mate,linkss,col=collinkss,lwd=lwdt)}}
  }
  
  out<-list(mate=mate,mate2=mate2,centros=M,jacs=jac,detjac=myj,detmean=detmean,obs=obs,obs2=obs2,obs3=obs3,fit=fit,cols=makeTransparent(colorRampPalette(colors)(n = length(obs3)),alpha=alpha),tes=tes,obm=obm,sourcem=plotsurf(source,t(triang),plot=F))
  out
}
#' @export

plotsurf<-function(lmatrix,triang,col=1,alpha=0.5,plot=T){
thesurf_1=t(cbind(lmatrix,1))
 triangolazioni=triang
  thesurf=list(vb=thesurf_1,it=triangolazioni)
  class(thesurf) <- "mesh3d"
if(plot==T){shade3d(thesurf,col=col,alpha=alpha)}
return(thesurf)
}
#' @export

tessell3d=function(tri,ver,iter){
tri_t=tri
ver_t=ver
for(j in 1:iter){
nrows=nrow(tri_t)
numb=range(tri_t)[2]
for(i in 1:nrows){
tri_i=ver_t[tri_t[i,],]
cen_i=colMeans(tri_i) 
tri_1=c(tri_t[i,c(1,2)],numb+i)
tri_2=c(tri_t[i,1],numb+i,tri_t[i,3])
tri_3=c(tri_t[i,c(2,3)],numb+i)
tri_t=rbind(tri_t,tri_1,tri_2,tri_3)
ver_t=rbind(ver_t,cen_i)
}}
vertici_out=cbind(ver_t,1)
rownames(vertici_out)=NULL
triangoli_out=tri_t[(dim(tri)[1]+1):dim(tri_t)[1],]
rownames(triangoli_out)=NULL
mesh.out=list(vb=t(vertici_out),it=t(triangoli_out))
}
#' @export

plotdefo3d<-function(procsymobject,triang=NULL,links=NULL,collinks=1,lwd=1,axtit=NULL,mags=rep(1,3),heat=F,scaleramp=F,heateuc=F,sign=T,colors=c("blue4","cyan2","yellow","red4"),alpha=1,from=NULL,to=NULL,out=F,plot=T){
  if(is.null(axtit)==T){axtit=c("PC1+","PC2+","PC3+","PC1-","PC2-","PC3-")}else{axtit=axtit}
  texts<-axtit
  if(!is.null(triang)==T&&ncol(triang)>3){stop("I want triangles matrix in the form nx3")}
  
  mshape<-procsymobject$mshape
  pc1pos<-showPC(max(procsymobject$PCscores[,1])*mags[1],procsymobject$PCs[,1],procsymobject$mshape)
  pc1neg<-showPC(min(procsymobject$PCscores[,1])*mags[1],procsymobject$PCs[,1],procsymobject$mshape)
  pc2pos<-showPC(max(procsymobject$PCscores[,2])*mags[2],procsymobject$PCs[,2],procsymobject$mshape)
  pc2neg<-showPC(min(procsymobject$PCscores[,2])*mags[2],procsymobject$PCs[,2],procsymobject$mshape)
  pc3pos<-showPC(max(procsymobject$PCscores[,3])*mags[3],procsymobject$PCs[,3],procsymobject$mshape)
  pc3neg<-showPC(min(procsymobject$PCscores[,3])*mags[3],procsymobject$PCs[,3],procsymobject$mshape)
  
  if(heat==F&heateuc==F){
    
    
    allshapes<-rbind(pc1pos,pc1neg,pc2pos,pc2neg,pc3pos,pc3neg)
    allshapes<-matrix2arrayasis(allshapes,nrow(pc1pos))
    css<-apply(allshapes,3,cSize)
    themax<-allshapes[,,which.max(css)]
    
    if(plot==T){
    open3d(windowRect=c(100,100,1000,1000)) 
    mat <- matrix(1:12, ncol=2)
    layout3d(mat, height = rep(c(3,1), 3), model = "inherit")
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    plot3d(pc1pos,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
    if(is.null(links)==F){lineplot(pc1pos,links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[1])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    plot3d(pc2pos,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
    if(is.null(links)==F){lineplot(pc2pos,links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[2])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    plot3d(pc3pos,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
    if(is.null(links)==F){lineplot(pc3pos,links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[3])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    plot3d(pc1neg,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
    if(is.null(links)==F){lineplot(pc1neg,links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[4])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    plot3d(pc2neg,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
    if(is.null(links)==F){lineplot(pc2neg,links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[5])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    plot3d(pc3neg,bbox=F,type="s",asp=F,axes=F,box=F,size=0.6,xlab = "", ylab = "", zlab = "",add=T)
    if(is.null(links)==F){lineplot(pc3neg,links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[6])}
  }
  if (heat==T){
    
    myh1pos<-heat3d(pc1pos,mshape,triang,lines=F,linkss=links,graphics=F,from=from,to=to,colors=colors,alpha=alpha,scaleramp=scaleramp)
    myh1neg<-heat3d(pc1neg,mshape,triang,lines=F,linkss=links,graphics=F,from=from,to=to,colors=colors,alpha=alpha,scaleramp=scaleramp)
    myh2pos<-heat3d(pc2pos,mshape,triang,lines=F,linkss=links,graphics=F,from=from,to=to,colors=colors,alpha=alpha,scaleramp=scaleramp)
    myh2neg<-heat3d(pc2neg,mshape,triang,lines=F,linkss=links,graphics=F,from=from,to=to,colors=colors,alpha=alpha,scaleramp=scaleramp)
    myh3pos<-heat3d(pc3pos,mshape,triang,lines=F,linkss=links,graphics=F,from=from,to=to,colors=colors,alpha=alpha,scaleramp=scaleramp)
    myh3neg<-heat3d(pc3neg,mshape,triang,lines=F,linkss=links,graphics=F,from=from,to=to,colors=colors,alpha=alpha,scaleramp=scaleramp)
  }
  if (heateuc==T){
    myh1pos<-diffonmesh(mshape,pc1pos,t(triang),graph=F,from=from,to=to,rampcolors=colors,alphas=c(alpha,0.7),sign=sign)
    myh1neg<-diffonmesh(mshape,pc1neg,t(triang),graph=F,from=from,to=to,rampcolors=colors,alphas=c(alpha,0.7),sign=sign)
    myh2pos<-diffonmesh(mshape,pc2pos,t(triang),graph=F,from=from,to=to,rampcolors=colors,alphas=c(alpha,0.7),sign=sign)
    myh2neg<-diffonmesh(mshape,pc2neg,t(triang),graph=F,from=from,to=to,rampcolors=colors,alphas=c(alpha,0.7),sign=sign)
    myh3pos<-diffonmesh(mshape,pc3pos,t(triang),graph=F,from=from,to=to,rampcolors=colors,alphas=c(alpha,0.7),sign=sign)
    myh3neg<-diffonmesh(mshape,pc3neg,t(triang),graph=F,from=from,to=to,rampcolors=colors,alphas=c(alpha,0.7),sign=sign)
    
  }
  
  allshapes<-rbind(t(myh1pos$obm$colMesh$vb[-4,]),t(myh1neg$obm$colMesh$vb[-4,]),t(myh2pos$obm$colMesh$vb[-4,]),t(myh2neg$obm$colMesh$vb[-4,]),t(myh3pos$obm$colMesh$vb[-4,]),t(myh3neg$obm$colMesh$vb[-4,]))
  allshapes<-matrix2arrayasis(allshapes,nrow(t(myh1pos$obm$colMesh$vb[-4,])))
  css<-apply(allshapes,3,cSize)
  themax<-allshapes[,,which.max(css)]
  if(heat==T|heateuc==T){
    if(plot==T){
    open3d(windowRect=c(100,100,1000,1000)) 
    mat <- matrix(1:12, ncol=2)
    layout3d(mat, height = rep(c(3,1), 3), model = "inherit")
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    shade3d(myh1pos$obm$colMesh,alpha=alpha)
    if(is.null(links)==F){lineplot(t(myh1pos$obm$colMesh$vb[-4,]),links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[1])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    shade3d(myh2pos$obm$colMesh,alpha=alpha) 
    if(is.null(links)==F){lineplot(t(myh2pos$obm$colMesh$vb[-4,]),links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[2])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    shade3d(myh3pos$obm$colMesh,alpha=alpha) 
    if(is.null(links)==F){lineplot(t(myh3pos$obm$colMesh$vb[-4,]),links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[3])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    shade3d(myh1neg$obm$colMesh,alpha=alpha)
    if(is.null(links)==F){lineplot(t(myh1neg$obm$colMesh$vb[-4,]),links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[4])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    shade3d(myh2neg$obm$colMesh,alpha=alpha) 
    if(is.null(links)==F){lineplot(t(myh2neg$obm$colMesh$vb[-4,]),links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[5])
    next3d()
    plot3d(themax*1.2,box=F,axes=F,col="white",xlab="",ylab="",zlab="",type="s",size=0,aspect=F)
    shade3d(myh3neg$obm$colMesh,alpha=alpha) 
    if(is.null(links)==F){lineplot(t(myh3neg$obm$colMesh$vb[-4,]),links,col=collinks,lwd=lwd)}
    next3d()
    text3d(0,0,0, texts[6])}
  }
  dimnames(allshapes)[[3]]<-c("pc1pos","pc1neg","pc2pos","pc2neg","pc3pos","pc3neg")
  if(out==T){
    if(heat==T){outp<-list(pc1pos=myh1pos$obm$colMesh,pc1neg=myh1neg$obm$colMesh,pc2pos=myh2pos$obm$colMesh,pc2neg=myh2neg$obm$colMesh,pc3pos=myh3pos$obm$colMesh,pc3neg=myh3neg$obm$colMesh)}else{outp<-allshapes}
    }else{outp<-NULL}
  outp
}
#' @export

conslinks<-function(number,open=T){
k=seq(1:(number-1))
aw=NULL
for(i in k){
b=list(c(i,i+1))
aw<-c(aw,b)
}
if(open==T){return(aw)}else{
aw<-c(aw,list(c(1,number)))
}
aw
}
#' @export

list2array<-function(mylist){
require(abind)
final<-NULL
for(i in 1:length(mylist)){
temp<-array(mylist[[i]],dim=c(nrow(mylist[[i]]),ncol(mylist[[i]]),1))
final<-abind::abind(final,temp)
}
return(final)
}
#' @export

serpred<-function(dep,indep,polyn=1,length.out=10){
sizes<-seq(min(indep),max(indep),length.out=length.out)
sizesvetts<-NULL
thelm<-lm(dep~poly(indep,degree=polyn,raw=TRUE))
for(i in 1:length(sizes)){
if(polyn==1){sizesvettsi<-c(1,sizes[i])}else{sizesvettsi<-c(1,poly(c(sizes[i]),degree=polyn,raw=T))}
sizesvetts<-rbind(sizesvetts,sizesvettsi)}
thepredicts<-NULL
  for(j in 1:nrow(sizesvetts)){
  thepredictsj<-c(crossprod(coef(thelm),sizesvetts[j,]))
  thepredicts<-rbind(thepredicts,thepredictsj)
}
thepredicts
}
#' @export


plot2dhull<-function(matrix,group,scalevariable=1,asp=1,extl=F,legend=T,xl=NULL,yl=NULL,posl=c("topright"),labels=NULL,clabel=0,pch=19,lwd=0.5,colhull=NULL,col=as.numeric(group),xlab=NULL,ylab=NULL,grey=F,xlimi=range(matrix[,1])*1.3,ylimi=range(matrix[,2])*1.3,reorder=T,alpha=rep(0,nlevels(group))){
  library(MASS)
  library(compositions)
  library(spatstat)
  library(gdata)
  library(TeachingDemos)
  library(ade4)
  library(calibrate)
  if (!is.factor(group)) stop("'group' must be a factor")
  
  
  
  makeTransparent = function(..., alpha=0.5) {
    
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    
    alpha = floor(255*alpha)  
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    
    .makeTransparent = function(col, alpha) {
      rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    
    return(newColor)
    
  }
  
  x<-matrix[,1]
  y<-matrix[,2]
  
  if(reorder==T){group<-factor(group,levels=unique(group))}
  
  
  
  species<-as.numeric(group)
  
  col=col
  if(grey==T){col=col2grey(col)}
  
  
  
  coli<-data.frame(group,as.numeric(group),col,pch)
  if(is.null(colhull)){colhull=aggregate(coli[,-1],by=list(coli[,1]),mean)[,-3][,2]}else{colhull=colhull}
  
  colim<-cbind(aggregate(coli[,-1],by=list(coli[,1]),mean)[,-c(3,4)],colhull)
  
  
  
  
  
  plot(x,y,type="p",cex=scalevariable,col=as.character(coli[,3]),pch=pch,xlim=xlimi,ylim=ylimi,xlab=xlab,ylab=ylab,asp=asp)
  
  
  
  if(!is.null(labels)){textxy(x,y,labels)}else{NULL}
  for(i in 1:(max(species)))
  {
    abline(h=0,v=0)
    if(length(x[species==i])<3){NULL}else{
      hulli<-convexhull.xy(subset(cbind(x,y),species==i))
      par(new=T)
      plot(hulli,col=makeTransparent(colim[,3][i], 1,alpha=alpha[i]),lwd=lwd,lty=i,add=T,asp=asp)
      par(new=T)
      daplotx<-c(hulli$bdry[[1]][[1]][length(hulli$bdry[[1]][[1]])],hulli$bdry[[1]][[1]])
      daploty<-c(hulli$bdry[[1]][[2]][length(hulli$bdry[[1]][[2]])],hulli$bdry[[1]][[2]])
      par(new=T)
      plot(daplotx,daploty,lwd=lwd,type="l",lty=i,col=colim[,3][i],xlim=xlimi,ylim=ylimi,axes=F,xlab="",ylab="",asp=asp)
      
    }
  }
  
  
  
  opar <- par(mar = par("mar"))
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  on.exit(par(opar))
  
  
  meansx<-aggregate(matrix[,1],by=list(group),FUN=mean)
  meansy<-aggregate(matrix[,2],by=list(group),FUN=mean)
  
  
  if (clabel > 0) 
    for(i in 1:nlevels(group)){ scatterutil.eti(meansx[i,-1],meansy[i,-1], meansx[,1][i], clabel,coul=1)}
  
  
  if(legend==T){
    if(extl==T){
      x11()
      plot(x,y,col="white")
      
      legend(min(x),max(y),unique(coli[,1]), cex=1, col=unique(coli[,3]), pch=unique(coli[,4]),box.col="white")
      
    }else{
    if(is.null(xl)==T&is.null(yl)==T){
      
      legend(posl,legend=unique(coli[,1]), col=unique(coli[,3]), pch=unique(coli[,4]),bty='n')}else{
        
        legend(xl,yl,unique(coli[,1]),col=unique(coli[,3]),pch=unique(coli[,4]),bty='n')}
    
  }
  }
}

#' @export 

reparray<-function(array,n,type=c("sequ","block")){
  
  if(is.null(dimnames(array)[3])){dimnames(array)[3]<-list(c(1:dim(array)[3]))}else{dimnames(array)[3]<-dimnames(array)[3]}

require(abind)
myarray2<-NULL
steps<-n

if(type=="sequ"){
  
  
  for(i in 1:dim(array)[3]){
    temp1<-rep(list(array(array[,,i],dim=c(dim(array)[1],dim(array)[2],1))),steps)
    temp2<-list2array(temp1)
    myarray2<-abind::abind(myarray2,temp2)
  }
  return(myarray2)}else{NULL}
    
  
  if(type=="block"){
    temp1<-list2matrix(rep(array2list(array),n))
    temp2<-matrix2arrayasis(temp1,dim(array)[1])
    return(temp2)
  }else{NULL}
}
#' @export  

pt2dvaruno<-function(mua,mub,va,doopa=T,sss=T,tol=0.001){
  library(shapes)
  library(matrixcalc)
  if(doopa==T){
    theopamuamub<-procOPA(mub,mua,scale=F,reflect=F)
    mua<-theopamuamub$Bhat
    rotmua<-theopamuamub$R
    
    if(rotmua[1,1]+1<tol){mua<--mua}
    
    if(sss==F){
      va<-va-((matrix.trace(t(mua)%*%va)/matrix.trace(t(mua)%*%mua))*mua)
    }
    
    theopavamua<-procOPA(mua,va,scale=F)
    va<-theopavamua$Bhat
    rotva<-theopavamua$R
    
    if(rotva[1,1]+1<tol){va<--va}
    
  }else{
    if(sss==F){
      va<-va-((matrix.trace(t(mua)%*%va)/matrix.trace(t(mua)%*%mua))*mua)
    }
    rotmua<-NULL
    rotva<-NULL
  }
  muaz <- complex(real = mua[,1], imaginary = mua[,2])
  mubz <- complex(real = mub[,1], imaginary = mub[,2])
  vaz <- complex(real = va[,1], imaginary = va[,2])
  muazbarra<-muaz/(sqrt(muaz%*%(Conj(muaz))))
  mubzbarra<-mubz/(sqrt(mubz%*%(Conj(mubz))))
  if(sss==T){ 
    vbz<-vaz-1i*((Im(Conj(mubzbarra)%*%vaz))/(1+Conj(muazbarra)%*%mubzbarra))*(muazbarra+mubzbarra)}else{
      vbz<-vaz-1i*((Im(Conj(mubzbarra)%*%vaz))/(1+Conj(muazbarra)%*%mubzbarra))*(muazbarra+mubzbarra)-((Re(Conj(mubzbarra)%*%vaz))/(1+Conj(muazbarra)%*%mubzbarra))*(muazbarra+mubzbarra)}
  
  vazbarra<-vaz/(sqrt(vaz%*%(Conj(vaz))))
  eps<-Re(sqrt(2*Im(Conj(mubzbarra)%*%vazbarra)^2/(1+Conj(muazbarra)%*%mubzbarra)))
  
  eps2<-Re(sqrt(2*Mod(Conj(mubzbarra)%*%vazbarra)^2/(1+Conj(muazbarra)%*%mubzbarra)))
  
  vb<-cbind(Re(vbz),Im(vbz))
  out<-list(vb=vb,rotmua=rotmua,rotva=rotva,eps=eps,eps2=eps2)
  return(out)
}
#' @export

opaloop<-function(GM,array,scale=F,reflect=F){
library(abind)
k<-dim(array)[1]
  m<-dim(array)[2]
  n<-dim(array)[3]
  
  looped<-NULL
rots<-NULL
for(i in 1:n){

opai<-procOPA(GM,array[,,i],scale=scale,reflect=reflect)

loopedi<-array(opai$Bhat,dim=c(k,m,1))
looped<-abind::abind(looped,loopedi)  

rotsi<-opai$R
rots<-c(rots,list(rotsi))



}
out<-list(looped=looped,rots=rots)
out
}

#' @export

diffonmesh<-function(lmsource,lmtarget,triang,colsource=1,alphas=c(0.7,0.3),grid=F,aff=F,nonaff=F,distscale=1, scaleramp=F,from = NULL,to = NULL,tol=NULL,sign=F,title=T,displace = FALSE,plotsource=T, steps = 20,rampcolors = colorRamps::blue2green2red(steps - 1),graph=T){
  for ( i in c("Morpho","Rvcg","rgl")) {
    if (!require(i,character.only = TRUE))
      stop(paste0("please install package ",i))
  }
  lmsource <- rotonto(lmtarget,lmsource)$yrot
  thetarget <- t(cbind(lmtarget,1))
  thesource <- t(cbind(lmsource,1))
  getPlotNr <- length(which(c(aff,nonaff,grid) == TRUE))+1
  layout(matrix(1:getPlotNr,1,getPlotNr,byrow = T))
  triangolazioni <- triang
  thetarget_mesh <- list(vb=thetarget,it=triangolazioni)
  class(thetarget_mesh) <- "mesh3d"
  class(thetarget_mesh) <- "mesh3d"
  thesource_mesh<-list(vb=thesource,it=triangolazioni)
  class(thesource_mesh) <- "mesh3d"
  class(thesource_mesh) <- "mesh3d"
  if(grid==T){
    shade3d(thetarget_mesh)
    
    ###between pointclouds
    distvec <- sqrt(rowSums((lmsource-lmtarget)^2))/distscale
    if(graph==T){
    them<-meshDist(lmtarget,distvec = distvec,from=from,to=to,tol=tol,sign=sign,steps = 20,rampcolors =rampcolors,scaleramp=scaleramp)
    }else{
      them<-meshDist(lmtarget,distvec = distvec,from=from,to=to,tol=tol,sign=sign,steps = 20,rampcolors =rampcolors,shade=F,scaleramp=scaleramp)
      
    }
      
    ##if you want add a deformed cube
    if(graph==T){
    deformGrid3d(lmsource,lmtarget,type="p",size=10,ngrid = 5,add=T)
    title("Distance between point-clouds")}
  }
  
  ### between surfaces
  thetarget_warped <- tps3d(thesource_mesh,lmsource,lmtarget,lambda = 1e-8)
  distvec2 <- sqrt(rowSums((vert2points(thesource_mesh)-vert2points(thetarget_warped))^2))/distscale
  if(graph==T){
  if(plotsource==T){shade3d(thesource_mesh,col=colsource,alpha=alphas[2])}
  them<-meshDist(thetarget_mesh,thesource_mesh,distvec=distvec2,alpha=alphas[1],add=T,from=from,to=to,tol=tol,scaleramp=scaleramp,sign=sign,displace=displace,steps = 20,rampcolors =rampcolors)
  if(title==T){title3d("total")
    title("total")}
  }else{
    them<-meshDist(thetarget_mesh,distvec=distvec2,alpha=alphas[1],add=T,from=from,to=to,tol=tol,scaleramp=scaleramp,sign=sign,displace=displace,steps = 20,rampcolors =rampcolors,shade=F)
}
  
  
  
  if (nonaff || aff) {
    affinetrafo <- computeTransform(vert2points(thetarget_mesh),vert2points(thesource_mesh),type="a")
    affineshort <- applyTransform(thesource_mesh,affinetrafo)
    
    if(nonaff) {### visualize non-affine deform between surfaces 
      ##create affine transfrom to longnose.mesh to remove affine differences and calculate distance from affine transform to target
      
      distvec3 <- sqrt(rowSums((vert2points(affineshort)-vert2points(thetarget_warped))^2))/distscale
      if(graph==T){
      open3d()
      shade3d(thetarget_mesh,col=1,alpha=alphas[2])
      themna<-meshDist(thesource_mesh,distvec=distvec3,add=T,from=from,to=to,tol=tol,scaleramp=scaleramp,sign=sign,displace=displace,steps = 20,rampcolors =rampcolors)
      if(title==T){title3d("non affine")
        title("non-affine")}
      }else{
        themna<-meshDist(thesource_mesh,distvec=distvec3,add=T,from=from,to=to,tol=tol,scaleramp=scaleramp,sign=sign,displace=displace,steps = 20,rampcolors =rampcolors,shade=F)
        }
      
      
    }
    
    if(aff) {###  the affine transform looks like that:
      distaffnonaffon4<-sqrt(rowSums((vert2points(affineshort)-vert2points(thesource_mesh))^2))/distscale
      if(graph==T){
      open3d()
      shade3d(thetarget_mesh,col=1,alpha=alphas[2])
      thema<-meshDist(thesource_mesh,distvec=distaffnonaffon4,add=T,from=from,to=to,tol=tol,scaleramp=scaleramp,sign=sign,displace=displace,steps = 20,rampcolors =rampcolors)
      if(title==T){title3d("affine")
        title("affine")}
      }else{
        thema<-meshDist(thesource_mesh,distvec=distaffnonaffon4,add=T,from=from,to=to,tol=tol,scaleramp=scaleramp,sign=sign,displace=displace,steps = 20,rampcolors =rampcolors,shade=F)
      }
      
      
    }
  }else{themna<-NULL
  thema<-NULL
  }
  out<-list(obm=them,themna=themna,thema=thema)
}
#' @export

mycca<-function(x,y,pch=19,col=1,group=NULL,labels=NULL,extl=F,legend=T,xl=NULL,yl=NULL,posl=c("topright"),cex=1,xlab=NULL,ylab=NULL,asp=NULL){
  require(CCA)
  require(vegan)
  require(calibrate)
  
  myrda<-vegan::rda(y~x)
  rsqadj<-RsquareAdj(myrda)
  print(rsqadj)
  
  anovarda<-anova(myrda,step=2000)
  print(anovarda)
  
  if(is.null(xlab)==T){xlab="Independent"}else{xlab=xlab}
  if(is.null(ylab)==T){ylab="Dependent"}else{ylab=ylab}
  if (!is.null(group)){
    col=col
    species<-as.numeric(group)
    coli<-data.frame(group,as.numeric(group),col,pch)
    colim<-aggregate(coli[,-1],by=list(coli[,1]),mean)}else{
      coli<-data.frame(1,1,col,pch)
      colim<-aggregate(coli[,-1],by=list(coli[,1]),mean) 
    }
  thecca<-rcc(as.matrix(x),as.matrix(y),0.1,0.1)
  if(is.null(group)==T){plot(x,rcc(as.matrix(x),as.matrix(y),0.1,0.1)$scores$yscores[,1],col=col,pch=pch,cex=cex,xlab=xlab,ylab=ylab)}  
  if (!is.null(group)){
    
    plot2dhull(cbind(x,thecca$scores$yscores[,1]),group,scalevariable=cex,pch=pch,col=col,colhull=colim[,3],clabel=0,legend=F,reorder=F,xlimi=range(x),ylimi=range(thecca$scores$yscores[,1]),asp=asp,xlab=xlab,ylab=ylab)
  }
  if (!is.null(labels)){textxy(x,thecca$scores$yscores,labels)}
  
  if(legend==T){
    if(extl==T){
      x11()
      plot(x,y,col="white")
      
      legend(min(x),max(y),colim[,1], cex=1, col=colim[,3], pch=colim[,4],box.col="white")
      }else{
      if(is.null(xl)==T&is.null(yl)==T){
        legend(posl,legend=colim[,1], cex=1, col=colim[,3], pch=colim[,4],bty='n')}else{
        legend(xl,yl,colim[,1], cex=1, col=colim[,3], pch=colim[,4],bty='n')}
      }
  }
  
  if (!is.null(group)){
    
    if(dim(as.matrix(y))[2]>1){bygroup<-manymultgr(y,x,group)}else{
      bygroup<-plotancova(y,x,group,legend=F,plot=F)$sint_results}
  }else{bygroup=c("no group structure")}
  
  out<-list(xscores=thecca$scores$xscores,yscores=thecca$scores$yscores,AdjRsq=rsqadj,anovarda=anovarda,bygroup=bygroup)
  return(out)
}
#' @export

posfac<-function(factor){
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

pure}
#' @export

areasip<-function(matrix,links=NULL,PB=NA,PA=NA,S=NA,SB=NA,H=NA,V=3,a=NULL,q=NULL,Y=FALSE,j=FALSE,D=FALSE,St=Inf,Q=TRUE,graph=T,extpol=F){
  if(!is.null(links)){warning("Links **must** identify, among other structures, a closed contour")}
  library(geometry)
  library(tripack)
  library(geoR)
  library(gstat)
  library(sp)
  library(fields)
  library(RTriangle)
  mate<-matrix
  if(!is.null(links)&is.na(S)==T){S<-list2matrix(links)}
  
  if(extpol==T){
  step1<-pslg(mate,S=S,PB=PA,PA=PA,SB=SB,H=NA)
  posimpnoh<-RTriangle::triangulate(step1,V=V,a=NULL,S=St,q=q,Y=Y,j=j,D=D,Q=Q)
  step2<-matrix[unique(posimpnoh$E[ posimpnoh$EB%in%c(1)]),]
  maulinks<-posimpnoh$E[posimpnoh$EB%in%c(1),]
  newmaulinks<-maulinks
  newmaulinks[,1]<-c(1:nrow(maulinks))
  newmaulinks[,2]<-match(maulinks[,2],maulinks[,1])
  cycles2 <- function(links) {
    require(ggm)
    if (any(is.na(links))) stop("missing value/s in links")
    mat <- matrix(0L, max(links), max(links))
    mat[links] <- 1
    lapply(ggm::fundCycles(mat), function(xa) rev(xa[ , 1L]))
  }
  
  mypol<-step2[cycles2(newmaulinks)[[1]],]}else{mypol<-c("you do not want external polygon")}
  
  
  p<-pslg(mate,PB=PB,PA=PA,S=S,SB=SB,H=H)
  au<-RTriangle::triangulate(p,V=1,a=a,q=q,Y=Y,j=j,D=D,S=St,Q=Q)
  if(graph==T){
    plot(au,asp=1)
    lineplot(mate,links,col=3)
    points(au$P,pch=19)
    points(mate,pch=21,cex=2)
  }
  centros<-NULL
  areas<-NULL
  for(i in 1:nrow(au$T)){
    centrosi<-apply(au$P[au$T[i,],],2,mean)
    centros<-rbind(centros,centrosi)
    areasi<-convhulln(au$P[au$T[i,],],option="FA")$vol
    areas<-c(areas,areasi)
  }
  if(graph==T){points(centros,col=2)}
  M<-centros
  
  
  
  thecentr<-centroids(mate)
  dasomm<-NULL
  for(i in 1:nrow(M)){
    dasommi<-dist(rbind(thecentr,M[i,]))^2*areas[i]
    dasomm<-c(dasomm,dasommi)
  }
  ip<-sum(dasomm)
  are<-sum(areas)
  
  deltri<-au$T
  ptri<-au$P
  origcentros<-rbind(mate,M)
  
  out<-list(area=are,areas=areas,ip=ip,centros=centros,deltri=deltri,ptri=ptri,origcentros=origcentros,triangob=au,ext=mypol)
  out
}
#' @export


array2list<-function(array){
  

thelist<-NULL

for(i in 1:dim(array)[3]){

eli<-array[,,i]

thelist<-c(thelist,list(eli))
}
if(is.null(dimnames(array)[[3]])==F){names(thelist)<-dimnames(array)[[3]]}

return(thelist)
}
#' @export

array2mat<-function(array,inds,vars){
  if(class(array)=="matrix"){array<-array(array,dim=c(nrow(array),ncol(array),1))}
  X1 <-aperm(array,c(3,2,1))
  dim(X1)<- c(inds, vars)
  if(!is.null(dimnames(array)[3])){rownames(X1)<-unlist(dimnames(array)[3])}else{rownames(X1)<-c(1:nrow(X1))}
  return(X1)
}
#' @export

biharm.new<-function(vec1,vec2){
  dim<-length(vec1)
  n <- sqrt(sum((vec1-vec2)^2))
  if(dim<3){
  n^2*log(n^2)}else{
 -n
}
}
#' @export

centershapes<-function(array){
if(is.matrix(array)==T){array<-array(array,dim=c(nrow(array),ncol(array),1))}
  centros<-centroids(array)
k<-dim(array)[1]
m<-dim(array)[2]
n<-dim(array)[3]
prov<-array(rep(t(centros),each=k),dim=c(k,m,n))
array2<-array-prov
  return(array2)
}
#' @export

centroids<-function(array){
  meanmat<-function(mat){
    mean<-apply(mat,2,mean)
    mean
  }
  if(class(array)=="matrix"){array<-array(array,dim=c(dim(array)[1],dim(array)[2],1))}
centroidi<-t(apply(array,3,meanmat))
return(centroidi)
}
#' @export

firstsfac<-function(group){


firstsres<-NULL
for (i in levels(group)){
  
  firstresi<-min(which(group==i))
  firstsres<-rbind(firstsres,firstresi)
}

return(firstsres)

}
#' @export

heat2d<-function(init,fin,invc=F,constr=T,sen=F,logsen=T,nadd=5000,linkss=NULL,zlim=NULL,legend=T,plottarget=T,ext=0.1,collinkss=1,lwds=2,collinkst=2,lwdt=2,pchs=19,pcht=19,cexs=0.5,cext=0.5,colors=c("blue4","cyan2","yellow","red4"),alpha=1,ngrid=30,oma=c(5,5,5,5),mai=c(3,3,3,3),mar=c(3,3,3,3),mag=1,tol=0.1,graphics=T,PB=NA,PA=NA,S=NA,SB=NA,H=NA,V=3,a=NULL,q=NULL,Y=FALSE,j=FALSE,D=FALSE,St=Inf,Q=TRUE,pholes=NA){
  library(shapes)
  library(Morpho)
  library(matrixcalc)
  library(geometry)
  library(tripack)
  library(geoR)
  library(gstat)
  library(sp)
  library(fields)
  library(RTriangle)
  library(sp)
  library(alphahull)
  library(ggm)
  if(!is.na(H)&is.na(pholes)==T){stop("If you specify 'H' you need to specify 'pholes' also")}
  if(!is.na(pholes)&is.null(linkss)==T){stop("If you specify 'H' and 'pholes' you need to specify 'linkss' also")}
  
  
  mate<-fin
  mate2<-init
  
fin<-init+(fin-init)*mag
  
  if(!is.null(linkss)&is.na(S)==T){S<-list2matrix(linkss)}
  
  posimpfin<-areasip(fin,S=S,PB=PA,PA=PA,SB=SB,H=H,V=V,a=NULL,q=q,Y=Y,j=j,D=D,St=St,Q=Q,graph=F,extpol=F)
  posimp<-areasip(init,S=S,PB=PA,PA=PA,SB=SB,H=H,V=V,a=NULL,q=q,Y=Y,j=j,D=D,St=St,Q=Q,graph=F,extpol=F)
  pofin<-areasip(fin,S=S,PB=PA,PA=PA,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,graph=F,extpol=F)
  po<-areasip(init,S=S,PB=PA,PA=PA,SB=SB,H=H,V=V,a=a,q=q,Y=Y,j=j,D=D,St=St,Q=Q,graph=F,extpol=F)
  posimpnohfin<-areasip(fin,S=S,PB=PA,PA=PA,SB=SB,H=NA,V=V,a=NULL,q=q,Y=Y,j=j,D=D,St=St,Q=Q,graph=F,extpol=F)
  posimpnoh<-areasip(init,S=S,PB=PA,PA=PA,SB=SB,H=NA,V=V,a=NULL,q=q,Y=Y,j=j,D=D,St=St,Q=Q,graph=F,extpol=F)
 
  matr<-init
  M<-po$centros
  areas<-po$areas
  
tpsgrid<-tpsgridpaolo(init,fin,linksTT=linkss,linksYY=linkss,axes2d=T,ext=ext,graphics=F,ngrid=ngrid)
 
  veclist<-NULL
  for(j in 1:nrow(M)){
    vecj<-NULL
    for(i in 1:nrow(matr)){
      vec1ji<-2*(M[j,1]-matr[i,1])+2*(M[j,1]-matr[i,1])*log((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2)
      vec2ji<-2*(M[j,2]-matr[i,2])+2*(M[j,2]-matr[i,2])*log((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2)
      vecji<-c(vec1ji,vec2ji)
      vecj<-rbind(vecj,vecji)
    }
    veclist<-c(veclist,list(vecj))
  }
  jac<-NULL
  for(i in 1: length(veclist)){
    jaci<-t(tpsgrid$B)+t(tpsgrid$W)%*%veclist[[i]]
    jac<-c(jac,list(jaci))
  }
  
  deltus<-NULL
  dpl<-NULL
  sens<-NULL
  for(i in 1:length(jac)){
    deltusi<-jac[[i]]-diag(nrow(jac[[i]]))
    dpli<-0.5*(deltusi+t(deltusi)) 
    seni<-0.5*matrix.trace(t(dpli)%*%dpli)
    deltus<-c(deltus,list(deltusi))
    dpl<-c(dpl,list(dpli))
    sens<-c(sens,seni)
  }
  
  lsens<-log2(sens)
  sens2<-sens*areas
  
  jac2<-list2array(jac)
  mean11<-mean(jac2[1,1,])
  mean12<-mean(jac2[1,2,])
  mean21<-mean(jac2[2,1,])
  mean22<-mean(jac2[2,2,])
  
  detmean<-det(matrix(c(mean11,mean21,mean12,mean22),ncol=2))
  
  myj<-unlist(lapply(jac,det))
  
  
  if(sen==F){obs<-log2(myj/detmean)}else{
    
    if(logsen==T){obs=lsens}else{obs=sens}
    
  }
  
  fit<- Tps(M, obs) 
  summary(fit)
  if(constr==F){
    hull <- chull(fin)
    indx=hull
    if(invc==F){points <- fin[indx,]}else{points <- init[indx,]}
    points <- rbind(points,points[1,])
    sfe1<-spsample(Polygon(points),nadd,type="regular")
    sr1<-sfe1@coords}else{
      
      
      if(invc==F){mau<-fin[unique(posimpnoh$triangob$E[posimpnoh$triangob$EB%in%c(1)]),]}else{mau<-init[unique(posimpnoh$triangob$E[posimpnoh$triangob$EB%in%c(1)]),]}
      matepro<-fin
      matepro[1:nrow(fin)%in%unique(posimpnoh$triangob$E[posimpnoh$triangob$EB%in%c(1)])==F,]<-NA
      faclinks<-as.factor(is.na(matepro[,1]))
      newindlinks<-posfac(faclinks)
      maulinks<-posimpnoh$triangob$E[posimpnoh$triangob$EB%in%c(1),]
      #plotmyarrays(mate,links=array2list(matrix2arrayasis(maulinks,1)))
      newmaulinks<-maulinks
      newmaulinks[,1]<-c(1:nrow(maulinks))
      newmaulinks[,2]<-match(maulinks[,2],maulinks[,1])
      #plotmyarrays(mau,links=array2list(matrix2arrayasis(newmaulinks,1)),xlim=c(25,70))
      
      cycles2 <- function(links) {
        require(ggm)
        if (any(is.na(links))) stop("missing value/s in links")
        mat <- matrix(0L, max(links), max(links))
        mat[links] <- 1
        lapply(ggm::fundCycles(mat), function(xa) rev(xa[ , 1L]))
      }
      
      #plotmyarrays(mau[cycles2(newmaulinks)[[1]],],links=conslinks(nrow(mau),open=F))
      
      
      
      sfe1<-spsample(Polygon(rbind(mau[cycles2(newmaulinks)[[1]],],mau[cycles2(newmaulinks)[[1]],][1,])),nadd,type="regular")
      sr1<-sfe1@coords
    }
  #points(sr1)
  
  pred<-predict(fit,sr1)
  if(is.null(zlim)==T){zlim<-range(pred)}else{zlim<-zlim}
  if(graphics==T){
    par(oma=oma,mai=mai,mar=mar)
    if(legend==F){ima<-image(xyz2img(cbind(sr1,pred),tolerance=tol),asp=1,xlim=tpsgrid$grid$xlim,ylim=tpsgrid$grid$ylim,col=makeTransparent(colorRampPalette(colors)(n = length(obs)),alpha=alpha),xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)}else{
      ima<-image.plot(xyz2img(cbind(sr1,pred),tolerance=tol),asp=1,xlim=tpsgrid$grid$xlim,ylim=tpsgrid$grid$ylim,col=makeTransparent(colorRampPalette(colors)(n = length(obs)),alpha=alpha),xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)  
    }
    par(new=T)
    plot(fin,xlim=tpsgrid$grid$xlim,ylim=tpsgrid$grid$ylim,asp=1,pch=pchs,cex=cexs,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
    if(is.null(linkss)==F){lineplot(fin,linkss,col=collinkss,lwd=lwds)}
    lines(tpsgrid$grid$ngrid,xlim=tpsgrid$grid$xlim,ylim=tpsgrid$grid$ylim)
    lines(tpsgrid$grid$ngrid2,xlim=tpsgrid$grid$xlim,ylim=tpsgrid$grid$ylim)
    
    if(plottarget==T){
      par(new=T)
      plot(init,xlim=tpsgrid$grid$xlim,ylim=tpsgrid$grid$ylim,asp=1,pch=pcht,cex=cext,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
      lineplot(init,linkss,col=collinkst,lwd=lwdt)
    }
  }
  
  pols<-NULL
  pols2<-NULL
  if(!is.na(pholes)){
    for(i in 1:length(pholes)){
      myindex<-which(apply(matrix(list2matrix(linkss)%in%pholes[[i]],ncol=2),1,sum)>1)
      matili<-list2matrix(linkss)[myindex,]
      
      holi<-init[unique(sort(c(unique(matili[,1]),unique(matili[,2])))),]
      hol2i<-fin[unique(sort(c(unique(matili[,1]),unique(matili[,2])))),]
      si<-matrix(as.numeric(ordered(matili)),ncol=2)
      holtriangi<-RTriangle::triangulate(pslg(holi,S=si))
      holtriang2i<-RTriangle::triangulate(pslg(hol2i,S=si))
      
      holui<-holi[unique(holtriangi$E[ holtriangi$EB%in%c(1)]),]
      holu2i<-hol2i[unique( holtriang2i$E[ holtriang2i$EB%in%c(1)]),]
      
      holiproi<-holi
      holiproi[1:nrow(holi)%in%unique(holtriangi$E[holtriangi$EB%in%c(1)])==F,]<-NA
      faclinkshi<-as.factor(is.na(holiproi[,1]))
      newindlinkshi<-posfac(faclinkshi)
      hlinks<-holtriangi$E[holtriangi$EB%in%c(1),]
      newhlinks<-hlinks
      newhlinks[,1]<-c(1:nrow(hlinks))
      newhlinks[,2]<-match(hlinks[,2],hlinks[,1])
      pols<-c(pols,list(holui[cycles2(newhlinks)[[1]],]))
      pols2<-c(pols2,list(holu2i[cycles2(newhlinks)[[1]],]))
      
      if(graphics==T){
        polygon(holui[cycles2(newhlinks)[[1]],],col="white",border=collinkss)
        polygon(holu2i[cycles2(newhlinks)[[1]],],col="white",border=collinkst)}
    }}
  
  
  out<-list(mate=fin,mate2=init,centros=M,interpcoords=sr1,jacs=jac,detjac=myj,detmean=detmean,obs=obs,fit=fit,pred=pred,xlim=tpsgrid$grid$xlim,ylim=tpsgrid$grid$ylim,tol=tol,cols=makeTransparent(colorRampPalette(colors)(n = length(obs)),alpha=alpha),tpsgrid=tpsgrid,sumse=sum(sens2),pols=pols,pols2=pols2,areasipob=po)
  out
}




#' @export


helmert<-function(p)
{H<-matrix(0, p, p)
 diag(H)<--(0:(p-1)) * (-((0:(p-1))*((0:(p-1))+1))^(-0.5))
for (i in 2:p){H[i,1:(i-1)]<- -((i-1)*(i))^(-0.5)}
H[1,]<-1/sqrt(p)
 H}
#' @export

lastsfac<-function(group){
  lastsres<-NULL
  for (i in levels(group)){
    lastresi<-max(which(group==i))
    lastsres<-rbind(lastsres,lastresi)
  }
  return(lastsres)
}
#' @export

list2matrix<-function(mylist){
  final<-NULL
  for(i in 1:length(mylist)){
    temp<-mylist[[i]]
    final<-rbind(final,temp)
  }
  #if(is.null(names(mylist))==F){rownames(final)<-names(mylist)}
  return(final)
}               
#' @export               
               
makeTransparent = function(..., alpha=0.5) {

  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")

  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)

  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }

  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)

  return(newColor)

}
#' @export

manymultgr<-function(y,x,group,steps=5000){
res<-NULL
rsq<-NULL
adjrsq<-NULL
pval<-NULL

for(i in 1:max(as.numeric((group)))){
  
  rsqi<-RsquareAdj(vegan::rda(y[as.numeric(group)==i,],x[as.numeric(group)==i]))$r.squared
  adjrsqi<-RsquareAdj(vegan::rda(y[as.numeric(group)==i,],x[as.numeric(group)==i]))$adj.r.squared
  pvali<-anova(vegan::rda(y[as.numeric(group)==i,],x[as.numeric(group)==i]),step=steps)$Pr[1]
  
  rsq<-append(rsq,rsqi)
  adjrsq<-append(adjrsq,adjrsqi)
  pval<-append(pval,pvali)
}
Sig<-ifelse(pval<0.05,c("*"),c(""))
res<-data.frame(round(cbind(rsq,adjrsq,pval),digits=5),Sig)
rownames(res)<-levels(group)
return(res)

}
#' @export

mopa<-function(tt,yy,rot=c("mopa","opa"),CSinit = FALSE,volinit=FALSE,center=TRUE){
  warning("the SECOND input matrix will be rotated on the FIRST one")
  if(is.matrix(tt)==T){
    k<-nrow(tt)
    m<-ncol(tt)
    ttar<-array(tt,dim=c(k,m,1))
    yyar<-array(yy,dim=c(k,m,1))}else {
      ttar=tt
      yyar=yy
      k<-dim(tt)[1]
      m<-dim(tt)[2]
    }
  
  if(center==T){ttcs<-centershapes(ttar)[,,1]
  yycs<-centershapes(yyar)[,,1]} 
  
  if(CSinit==T){
yycs<-scaleshapes(yycs)
ttcs<-scaleshapes(ttcs)
    }else{yycs<-yycs;ttcs<-ttcs}
  
  if(volinit==T){
  at<-tpsdry2(ttcs,yycs,meth="mor",g11=F)$at
  detat<-det(at)
  if(ncol(yycs)>2){yycs<-yycs/(detat^(1/3))}else{yycs<-yycs/(detat^(1/2))}
  }
  
  if(rot=="opa"){rotprocstep1<-svd(t(ttcs)%*%yycs)}
  stef<-CreateL(ttcs)
  W<-(stef$Lsubk%*%yycs)
  appo<-stef$Linv[(k+1):(k+(m+1)),1:k]%*%yycs
  cT<-appo[1,]
  At<-appo[-1,]
  if(rot=="mopa"){
    rotprocstep1<-svd(At)}
  rotprocstep2<-rotprocstep1$v%*%t(rotprocstep1$u)
  opizzata1<-yycs%*%rotprocstep2
  opizzata2<-array(opizzata1,dim=c(k,m,1))[,,1,drop=F]####centershapes(array(opizzata1,dim=c(k,m,1)))[,,1,drop=F]
  Wafter<-stef$Lsubk%*%opizzata2[,,1]
  appoafter<-stef$Linv[(k+1):(k+(m+1)),1:k]%*%opizzata2[,,1]
  cTafter<-appoafter[1,]
  Atafter<-appoafter[-1,]
  out=list(opizzata=opizzata2,W=W,cT=cT,At=At,Wafter=Wafter,cTafter=cTafter,Atafter=Atafter)
  return(out)
  }
#' @export

newmb<-function(source,target){
k<-nrow(source)
m<-ncol(source)
h<-1*(helmert(k)[-1,])
cm<-t(h)%*%h
source<-cm%*%source
target<-cm%*%target
s1<-sigm.new(source)
source<-h%*%source
target<-h%*%target
s<-h%*%s1%*%t(h)
S<-s
Q<-source
Gam_11=solve(S)-solve(S)%*%Q%*%solve(t(Q)%*%solve(S)%*%Q)%*%t(Q)%*%solve(S)
Gam_21=solve(t(Q)%*%solve(S)%*%Q)%*%t(Q)%*%solve(S)
Gam_22=-solve(t(Q)%*%solve(S)%*%Q)

W=Gam_11%*%target
At=Gam_21%*%target
out<-list(gamma11=Gam_11,gamma21=Gam_21,gamma22=Gam_22,W=W,At=At,sold=s1,snew=s,h=h)
return(out)
}
#' @export


opaloop2<-function(GM,array,scale=F,reflect=F){
  library(abind)
library(Morpho)
  k<-dim(array)[1]
  m<-dim(array)[2]
  n<-dim(array)[3]
looped<-NULL
  rots<-NULL
  for(i in 1:n){
print(i)
opai<-rotonto(GM,array[,,i],scale=scale,reflection=reflect,signref=T)
loopedi<-array(opai$yrot,dim=c(k,m,1))
    looped<-abind::abind(looped,loopedi)  
 rotsi<-opai$gamm
    rots<-c(rots,list(rotsi))
}
  out<-list(looped=looped,rots=rots)
  out
}
#' @export

ordiwithshapes<-function(mshape,scores,rotation,whichaxes=c(1,2),addata=NULL,asp=1,xlab="PC1",ylab="PC2",triang=NULL,factraj=NULL,coltraj=NULL,distscale=1,scaleramp=F,from=0,to=NULL,mag=1,procSym=T,subplotdim=2,legendgroup=NULL,shiftposx=1.5,shiftnegx=1.5,shiftposy=1.5,shiftnegy=1.5,links=NULL,collinks=1,grouphull=NULL,colhull=NULL,labels=NULL,labelscex=1,pch=19,col=1,cex=1,mult=0.5,tit3d=T,plotsource=T,xlim=extendrange(extendrange(range(scores[,whichaxes[1]]),f=mult)),ylim=extendrange(extendrange(range(scores[,whichaxes[2]]),f=mult))){
  library(Morpho)
  library(spatstat)
  library(calibrate)
  
  if(procSym==T){
    
    pc1pos<-showPC(max(scores[,whichaxes[1]])*mag,rotation[,whichaxes[1]],mshape)
    pc1neg<-showPC(min(scores[,whichaxes[1]])*mag,rotation[,whichaxes[1]],mshape)
    pc2pos<-showPC(max(scores[,whichaxes[2]])*mag,rotation[,whichaxes[2]],mshape)
    pc2neg<-showPC(min(scores[,whichaxes[2]])*mag,rotation[,whichaxes[2]],mshape)}else{
      
      
      
      pc1pos<-showPC2(max(scores[,whichaxes[1]])*mag,rotation[,whichaxes[1]],mshape)
      pc1neg<-showPC2(min(scores[,whichaxes[1]])*mag,rotation[,whichaxes[1]],mshape)
      pc2pos<-showPC2(max(scores[,whichaxes[2]])*mag,rotation[,whichaxes[2]],mshape)
      pc2neg<-showPC2(min(scores[,whichaxes[2]])*mag,rotation[,whichaxes[2]],mshape)
    }
  
  
  if(dim(mshape)[2]<3){
    
    
    library(TeachingDemos)
    library(Morpho)
    library(calibrate)
    library(spatstat)
    library(shapes)
    plot(scores[,whichaxes[1]],scores[,whichaxes[2]],axes=F,xlim=xlim,ylim=ylim,pch=pch,col=col,cex=cex,asp=asp,xlab=xlab,ylab=ylab)
    axis(1,pos=0,at=round(seq(from=min(scores[,whichaxes[1]]),to=max(scores[,whichaxes[1]]),length.out=10),2))
    axis(2,pos=0,at=round(seq(from=min(scores[,whichaxes[2]]),to=max(scores[,whichaxes[2]]),length.out=10),2))
    if(!is.null(addata)){
      par(new=T)
      plot(addata,asp=asp,xlim=xlim,ylim=ylim,pch=pch,col=col,cex=cex,axes=F,xlab="",ylab="")}
    
    
    
    if(!is.null(factraj)){
      if(is.null(coltraj)){coltraj=rep(1,nlevels(factraj))}else{coltraj=coltraj}
      
      
      for (i in 1 :nlevels(factraj)){
        lines(scores[,whichaxes[1]][as.numeric(factraj)==i],scores[,whichaxes[2]][as.numeric(factraj)==i],xlim=xlim,ylim=ylim,col=coltraj[i])
        par(new=T)
      } 
    }
    
    
    
    
    
    if(!is.null(labels)==T){textxy(scores[,whichaxes[1]],scores[,whichaxes[2]],labs=labels,cex=labelscex)}
    
    if(!is.null(grouphull)==T){
      if(is.null(colhull)==T){colhull<-c(1:nlevels(grouphull))}else{colhull<-colhull}
      for(i in 1:(max(as.numeric(grouphull))))
      {if(length(scores[,1][as.numeric(grouphull)==i])<3){NULL} 
        else{
          hulli<-convexhull.xy(subset(cbind(scores[,whichaxes[1]],scores[,whichaxes[2]]),as.numeric(grouphull)==i)) 
          par(new=T)
          daplotx<-c(hulli$bdry[[1]][[1]][length(hulli$bdry[[1]][[1]])],hulli$bdry[[1]][[1]])
          daploty<-c(hulli$bdry[[1]][[2]][length(hulli$bdry[[1]][[2]])],hulli$bdry[[1]][[2]])
          
          plot(daplotx,daploty,lwd=2,type="l",lty=i,col=colhull[i],xlim=xlim,ylim=ylim,axes=F,asp=asp,xlab="",ylab="")
        }
      }}
    
    subplot(tpsgridpaolo(mshape,pc1pos,linksYY=links,collinksYY=collinks,axes2d=F,displ=F,collinksTT=makeTransparent("white",alpha=0)),max(scores[,whichaxes[1]])*shiftposx,0,size=c(subplotdim,subplotdim))
    subplot(tpsgridpaolo(mshape,pc1neg,linksYY=links,collinksYY=collinks,axes2d=F,displ=F,collinksTT=makeTransparent("white",alpha=0)),min(scores[,whichaxes[1]])*shiftnegx,0,size=c(subplotdim,subplotdim))
    subplot(tpsgridpaolo(mshape,pc2pos,linksYY=links,collinksYY=collinks,axes2d=F,displ=F,collinksTT=makeTransparent("white",alpha=0)),0,max(scores[,whichaxes[2]])*shiftposy,size=c(subplotdim,subplotdim))
    subplot(tpsgridpaolo(mshape,pc2neg,linksYY=links,collinksYY=collinks,axes2d=F,displ=F,collinksTT=makeTransparent("white",alpha=0)),0,min(scores[,whichaxes[2]])*shiftnegy,size=c(subplotdim,subplotdim))
    
  }else{
    
    
    
    if(is.null(triang)==F&is.null(links)==T){
      open3d(windowRect=c(100,100,1000,1000)) 
      mat <- matrix(1:8, ncol=2)
      layout3d(mat, height = rep(c(2,1), 2), model = "inherit")
      diffonmesh(mshape,pc1pos,triang,title=F,distscale=distscale,from=from,to=to,plotsource=plotsource,scaleramp=scaleramp)
      next3d()
      text3d(0,0,0,paste("Scores",whichaxes[1],"pos"))
      next3d()
      diffonmesh(mshape,pc2pos,triang,title=F,distscale=distscale,from=from,to=to,plotsource=plotsource,scaleramp=scaleramp)
      next3d()
      text3d(0,0,0,paste("Scores",whichaxes[2],"pos"))
      next3d() 
      diffonmesh(mshape,pc1neg,triang,title=F,distscale=distscale,from=from,to=to,plotsource=plotsource,scaleramp=scaleramp)
      next3d();text3d(0,0,0,paste("Scores",whichaxes[1],"neg"))
      next3d()
      diffonmesh(mshape,pc2neg,triang,title=F,distscale=distscale,from=from,to=to,plotsource=plotsource,scaleramp=scaleramp)
      next3d();text3d(0,0,0,paste("Scores",whichaxes[2],"neg"))
    }
    
    if(is.null(links)==F){
      open3d(windowRect=c(100,100,1000,1000)) 
      mat <- matrix(1:8, ncol=2)
      layout3d(mat, height = rep(c(2,1), 2), model = "inherit")
      plot3D(pc1pos,add=T,bbox=F)
      lineplot(pc1pos,links)
      if(tit3d==T){next3d();text3d(0,0,0,paste("Scores",whichaxes[1],"pos"))}
      next3d()
      plot3D(pc2pos,add=T,bbox=F)
      lineplot(pc2pos,links)
      if(tit3d==T){next3d();text3d(0,0,0,paste("Scores",whichaxes[2],"pos"))}
      next3d()
      plot3D(pc1neg,add=T,bbox=F)
      lineplot(pc1neg,links)
      if(tit3d==T){next3d();text3d(0,0,0,paste("Scores",whichaxes[1],"neg"))}
      
      next3d()
      plot3D(pc2neg,add=T,bbox=F)
      lineplot(pc2neg,links)
      if(tit3d==T){next3d();text3d(0,0,0,paste("Scores",whichaxes[2],"neg"))}}
    
    plot(scores[,whichaxes[1]],scores[,whichaxes[2]],axes=F,xlim=xlim,ylim=ylim,pch=pch,col=col,cex=cex,asp=asp,xlab=xlab,ylab=ylab)
    axis(1,pos=0,at=round(seq(from=min(scores[,whichaxes[1]]),to=max(scores[,whichaxes[1]]),length.out=10),2))
    axis(2,pos=0,at=round(seq(from=min(scores[,whichaxes[2]]),to=max(scores[,whichaxes[1]]),length.out=10),2))
    if(!is.null(addata)){
      par(new=T)
      plot(addata,xlim=xlim,ylim=ylim,pch=pch,col=col,cex=cex,axes=F,asp=asp,xlab="",ylab="")}
    
    
    if(!is.null(factraj)){
      if(is.null(coltraj)){coltraj=rep(1,nlevels(factraj))}else{coltraj=coltraj}
      
      
      for (i in 1 :nlevels(factraj)){
        lines(scores[,whichaxes[1]][as.numeric(factraj)==i],scores[,whichaxes[2]][as.numeric(factraj)==i],xlim=xlim,ylim=ylim,col=coltraj[i])
        par(new=T)
      } 
    }
    
    
    
    if(!is.null(labels)==T){textxy(scores[,whichaxes[1]],scores[,whichaxes[2]],labs=labels)}
    
    if(!is.null(grouphull)==T){
      if(is.null(colhull)==T){colhull<-c(1:nlevels(grouphull))}else{colhull<-colhull}
      for(i in 1:(max(as.numeric(grouphull))))
      {if(length(scores[,1][as.numeric(grouphull)==i])<3){NULL} 
        else{
          hulli<-convexhull.xy(subset(cbind(scores[,whichaxes[1]],scores[,whichaxes[2]]),as.numeric(grouphull)==i)) 
          par(new=T)
          daplotx<-c(hulli$bdry[[1]][[1]][length(hulli$bdry[[1]][[1]])],hulli$bdry[[1]][[1]])
          daploty<-c(hulli$bdry[[1]][[2]][length(hulli$bdry[[1]][[2]])],hulli$bdry[[1]][[2]])
          
          plot(daplotx,daploty,lwd=2,type="l",lty=i,col=colhull[i],xlim=xlim,ylim=ylim,axes=F,asp=asp)
        }
      }}}
  
  
  
  
  
  if(is.null(legendgroup)==F){
    x11()
    plot(scores[,1],scores[,2],col="white")
    legend(min(scores[,1]),max(scores[,2]),unique(legendgroup), cex=1, col=unique(col), pch=unique(pch),box.col="white")}
  
  
}
#' @export

plotptau5<-function(objptau,links=NULL,projorig=T,whichaxes=c(1,2),cexscale=round(max(centroid.size(objptau$arrayord)),digits=0),shapescale=10,mag=1,subplotdim=2,shiftnegy=1,shiftposy=1,shiftnegx=1,shiftposx=1,col=as.numeric(objptau$factorord),pch=19,triang=NULL,from=0,topcs=NULL,tore=NULL,plotsource=T){
  library(TeachingDemos)
  k<-dim(objptau$arrayord)[1]
  m<-dim(objptau$arrayord)[2]
  n<-dim(objptau$arrayord)[3]
  transp<-read.inn(objptau$shifted,k,m)
  
  origsize<-objptau$indepepure
  
if(projorig==T){
ordiwithshapes(objptau$space1mshape,objptau$space1$x,objptau$space1$rotation,procSym=F,whichaxes=whichaxes,addata=objptau$origproj[,whichaxes],asp=1,factraj=objptau$factorord,cex=origsize/cexscale,triang=triang,from=from,to=topcs,links=links,mag=mag,shiftnegy=1.5,shiftposy=2,col=col,pch=pch,subplotdim=subplotdim,plotsource=plotsource)
  title("LS on predictions on Shape Space and original data re-projected")}else{

ordiwithshapes(objptau$space1mshape,objptau$space1$x,objptau$space1$rotation,procSym=F,whichaxes=whichaxes,asp=1,factraj=objptau$factorord,cex=origsize/cexscale,triang=triang,from=from,to=topcs,links=links,mag=mag,shiftnegy=1.5,shiftposy=2,col=col,pch=pch,subplotdim=subplotdim,plotsource=plotsource)
  title("LS on predictions on Shape Space")}

  
  
  ranges<-apply(rbind(objptau$space1$x,objptau$origproj),2,range)
  
  ratesorig<-ratesbygroup(read.inn(objptau$predictpure2,k,m),objptau$factorord,objptau$indepepure)
  plot(objptau$indepepure,ratesorig,pch=pch,col=col)
  for(i in 1:nlevels(objptau$factorord)){
    lines(objptau$indepepure[as.numeric(objptau$factorord)==i][-firstsfac(objptau$factorord)],ratesorig[as.numeric(objptau$factorord)==i][-firstsfac(objptau$factorord)],col=col[firstsfac(objptau$factorord)][i])
  }
  title("Rates of shape change among consecutive predictions per unit size in original data")
  
  ratestrasp<-ratesbygroup(objptau$predictpure,objptau$factorord,objptau$indepepure)
  
  plot(objptau$indepepure,ratestrasp,pch=pch,col=col)
  for(i in 1:nlevels(objptau$factorord)){
    
    lines(objptau$indepepure[as.numeric(objptau$factorord)==i][-firstsfac(objptau$factorord)],ratestrasp[as.numeric(objptau$factorord)==i][-firstsfac(objptau$factorord)],col=col[firstsfac(objptau$factorord)][i])
  }
  
  
  
  title("Rates of shape change among consecutive predictions per unit size in transported data")
  
  
  
  plot(ratesorig,ratestrasp,col=col,pch=pch)
  title("Original rates of shape change per unit size vs those of transported data")
  
  adults<-NULL
  adultsdaplot<-NULL
  for(i in 1:nlevels(objptau$factorord)){
    adulti<-array(showPC2(objptau$space1$x[lastsfac(objptau$factorord),1:3][i,],objptau$space1$rotation[,1:3],objptau$space1mshape),dim=c(k,m,1))
    adultdaploti<-array(showPC2(objptau$space1$x[lastsfac(objptau$factorord),1:3][i,]*mag,objptau$space1$rotation[,1:3],objptau$space1mshape),dim=c(k,m,1))
    
    adults<-abind::abind(adults,adulti)
    adultsdaplot<-abind::abind(adultsdaplot,adultdaploti)
  }
  
  if(m<3){adultsdaplot<-abind::abind(adultsdaplot,array(rep(0,k*n),dim=c(k,1,nlevels(objptau$factorord))),along=2)}
  
  init<-array(showPC2(objptau$space1$x[firstsfac(objptau$factorord),1:3][1,],objptau$space1$rotation[,1:3],objptau$space1mshape),dim=c(k,m,1))
  initdaplot<-array(showPC2(objptau$space1$x[firstsfac(objptau$factorord),1:3][1,]*mag,objptau$space1$rotation[,1:3],objptau$space1mshape),dim=c(k,m,1))
  
  if(m<3){initdaplot<-abind::abind(initdaplot,array(rep(0,k*n),dim=c(k,1,1)),along=2)}
  
  if(is.null(triang)==F){
  adultsdaplotvis<-NULL
  for(i in 1:nlevels(objptau$factorord)){
  adultsdaplotvisi<-meshDist(plotsurf(adultsdaplot[,,i],triang,plot=F),plotsurf(initdaplot[,,1],triang,plot=F),to=tore,from=from,plot=F)
  daagg<-objptau$space1$x[1:n,1:3][lastsfac(objptau$factorord),][i,]
  adultsdaplotvisi<-translate3d(scalemesh(adultsdaplotvisi$colMesh,1/shapescale,center="none"),daagg[1],daagg[2],daagg[3])
 adultsdaplotvis<-c(adultsdaplotvis,list(adultsdaplotvisi))
  }
  }
  
  
  
  
  open3d()
  for(i in 1:nlevels(objptau$factorord)){
    lines3d(objptau$space1$x[1:n,1:3][as.numeric(objptau$factorord)==i,],col=i,add=T)
  }
  
  if(is.null(triang)==T){
  babedaagg<-rep.row(objptau$space1$x[firstsfac(objptau$factorord),1:3][1,],k)
  for(i in 1:nlevels(objptau$factorord)){
    daplottarei<-adultsdaplot[,,i]
    daagg<-rep.row(objptau$space1$x[1:n,1:3][lastsfac(objptau$factorord),][i,],k)
    plot3D((daplottarei/shapescale)+daagg,col=i,add=T,size=1,bbox=F)
    lineplot((daplottarei/shapescale)+daagg,links,col=i)
  }
  plot3D((initdaplot[,,1]/shapescale)+babedaagg,bbox=F,add=T)
  lineplot((initdaplot[,,1]/shapescale)+babedaagg,links)
  
  }else{
    
    for(i in 1:nlevels(objptau$factorord)){
      shade3d(adultsdaplotvis[[i]],add=T)
    }
    babedaagg<-rep.row(objptau$space1$x[firstsfac(objptau$factorord),1:3][1,],k)
  shade3d(plotsurf((initdaplot[,,1]/shapescale)+babedaagg,triang,plot=F),add=T,alpha=0.5)
    }
  title3d("LS data in the Shape Space; first three PC scores")
   decorate3d()
  out<-list(ratesorig=ratesorig,ratestrasp=ratestrasp)
  
}
#' @export

pwpermancova<-function(y, pred1,pred2,nperm=999){
library(vegan)

if(is.factor(pred1)==TRUE){pred1<-factor(pred1,levels=unique(pred1))}else{NULL}
if(is.factor(pred2)==TRUE){pred2<-factor(pred2,levels=unique(pred2))}else{NULL}

if(is.factor(pred1)==TRUE){species<-as.numeric(pred1)}else{NULL}
if(is.factor(pred2)==TRUE){species<-as.numeric(pred2)}else{NULL}

if(is.factor(pred1)==TRUE){group<-pred1}else{NULL}
if(is.factor(pred2)==TRUE){group<-pred2}else{NULL}



fat_species<-group
r_adonis_pred1<-matrix(0,nrow=max(species),ncol=max(species))
p_adonis_pred1<-matrix(0,nrow=max(species),ncol=max(species))

r_adonis_pred2<-matrix(0,nrow=max(species),ncol=max(species))
p_adonis_pred2<-matrix(0,nrow=max(species),ncol=max(species))

r_adonis_pred1_pred2<-matrix(0,nrow=max(species),ncol=max(species))
p_adonis_pred1_pred2<-matrix(0,nrow=max(species),ncol=max(species))





for(i in 1:(max(species)-1))
{
for(j in (i+1):max(species)){

if(ncol(as.matrix(pred1))>1){
ado<-adonis(y[species==i|species==j,]~pred1[species==i|species==j,]*pred2[species==i|species==j],permutations=nperm,method="euclidean")
}else{NULL}

if(ncol(as.matrix(pred2))>1){
ado<-adonis(y[species==i|species==j,]~pred1[species==i|species==j]*pred2[species==i|species==j,],permutations=nperm,method="euclidean")
}else{NULL}

if(ncol(as.matrix(pred1))<2&ncol(as.matrix(pred2))<2){
ado<-adonis(y[species==i|species==j,]~pred1[species==i|species==j]*pred2[species==i|species==j],permutations=nperm,method="euclidean")
}else{NULL}




r_adonis_pred1[i,j]<-ado$aov.tab$R2[1]
p_adonis_pred1[i,j]<-ado$aov.tab$Pr[1]

r_adonis_pred2[i,j]<-ado$aov.tab$R2[2]
p_adonis_pred2[i,j]<-ado$aov.tab$Pr[2]


r_adonis_pred1_pred2[i,j]<-ado$aov.tab$R2[3]
p_adonis_pred1_pred2[i,j]<-ado$aov.tab$Pr[3]



  }
}

rownames(r_adonis_pred1)=colnames(r_adonis_pred1)<-levels(group)
rownames(p_adonis_pred1)=colnames(p_adonis_pred1)<-levels(group)
rownames(r_adonis_pred2)=colnames(r_adonis_pred2)<-levels(group)
rownames(p_adonis_pred2)=colnames(p_adonis_pred2)<-levels(group)
rownames(r_adonis_pred1_pred2)=colnames(r_adonis_pred1_pred2)<-levels(group)
rownames(p_adonis_pred1_pred2)=colnames(p_adonis_pred1_pred2)<-levels(group)



out<-list(r_adonis_pred1=r_adonis_pred1,p_adonis_pred1=p_adonis_pred1,r_adonis_pred2=r_adonis_pred2,p_adonis_pred2=p_adonis_pred2,r_adonis_pred1_pred2=r_adonis_pred1_pred2,p_adonis_pred1_pred2=p_adonis_pred1_pred2)
return(out)

}
#' @export

ratesbygroup<-function(predictions,factor,indep,sss=F){
  prova3<-seqrate(predictions,factor,vector=indep,sss=sss)
  prova4<-NULL
  for(i in 1:nlevels(factor)){
    prova4i<-c(rep(0,nlevels(factor))[i],prova3$rate[as.numeric(factor[-firstsfac(factor)])==i])
    prova4<-append(prova4,prova4i)}
  
  prova5<-prova4
  plot(indep,prova5,col=as.numeric(factor),ylim=c(0,max(prova5)),xlim=range(indep))
  for(i in 1:nlevels(factor)){
    lines(indep[as.numeric(factor)==i][-firstsfac(factor)],prova4[as.numeric(factor)==i][-firstsfac(factor)],col=i)
  }
  return(prova5)
}
#' @export

rep.row<-function(x,n){
matrix(rep(x,each=n),nrow=n)
}
#' @export

scaleshapes<-function(array){
  k<-dim(array)[1]
  m<-dim(array)[2]
  n<-dim(array)[3]
  library(abind)
  library(shapes)
  library(Morpho)
  divmat<-function(mat){
    mat2<-mat/cSize(mat)
    mat2
  }
  if(is.matrix(array)==T){array<-array(array,dim=c(nrow(array),ncol(array),1))}
scaled<-vecx(t(apply(array,3,divmat)),revert=T,lmdim=m)
return(scaled)
}
#' @export


seqrate<-function(preds,factor,vector=centroid.size(preds),sss=F){
  mytable<-table(factor)
  chrate<-NULL
  sizes<-NULL
  seqdist<-NULL
  for(i in 1: nlevels(factor)){
    for(j in 2:mytable[i]){
      sizei<-(vector[as.numeric(factor)==i][j])-(vector[as.numeric(factor)==i][j-1])
      
      if(sss==T){
      seqdisti<-ssriemdist2(preds[,,as.numeric(factor)==i][,,j],preds[,,as.numeric(factor)==i][,,j-1])}else{
      
      seqdisti<-kendalldist(preds[,,as.numeric(factor)==i][,,j],preds[,,as.numeric(factor)==i][,,j-1])
      }
      
      
      
      chratei<-seqdisti/(sizei)
      sizes<-append(sizes,sizei)
      chrate<-append(chrate,chratei)
      seqdist<-append(seqdist,seqdisti)
    }
  }
  out<-list(rate=chrate,sizes=sizes,seqdist=seqdist)
  return(out)
}
#' @export



showPC2<-function(x,rotation,mshape) {
  k<-dim(mshape)[1]
m<-dim(mshape)[2]
  mshape2<-array2mat(array(mshape,dim=c(k,m,1)),1,k*m)
  dims <- dim(mshape)
  if (length(x) > 1){predPC<-c(rotation %*% x)}else{predPC<-rotation*x}
  modell <- read.inn(mshape2+predPC,k,m)[,,1]
  return(modell)
}
#' @export

sigm.new<-function(mat1,mat2=NULL){
  tt<-mat1
  if(is.null(mat2)){yy<-mat1}else{yy<-mat2}
  SGMr <- apply(tt,1,function(t)apply(yy,1,biharm.new,t))
  replace(SGMr,which(SGMr=="NaN",arr.ind=T),0)
}
#' @export

tpsgridpaolo<-function (TT, YY, xbegin = -999, ybegin = -999, xlim=NULL,ylim=NULL,xwidth = -999, opt = 1, ext = 0.0, asp=1,ngrid = 22, cex = 1, pch = 20,colshift=1,zslice = 0, graphics=T,mag = 1, axes3 = FALSE,linksTT=NULL,linksYY=NULL,collinksTT=1,collinksYY=2,lwdtt=2,lwdyy=2,colgrid=1,axes2d=T,collandsTT=collinksTT,collandsYY=collinksYY,displ=T,coldispl=4){
  ######### some SMALL changes from Ian Dryden's function from package "shapes"
  k <- dim(TT)[1]
  m <- dim(TT)[2]
  YY <- TT + (YY - TT) * mag
  bb <- array(TT, c(dim(TT), 1))
  aa <- defplotsize2(bb)
  if (xwidth == -999) {
    xwidth <- aa$width
  }
  if (xbegin == -999) {
    xbegin <- aa$xl
  }
  if (ybegin == -999) {
    ybegin <- aa$yl
  }
  if (m == 3) {
    zup <- max(TT[, 3])
    zlo <- min(TT[, 3])
    zpos <- zslice
    for (ii in 1:length(zslice)) {
      zpos[ii] <- (zup + zlo)/2 + (zup - zlo)/2 * zslice[ii]
    }
  }
  xstart <- xbegin
  ystart <- ybegin
  ngrid <- trunc(ngrid/2) * 2
  kx <- ngrid
  ky <- ngrid - 1
  l <- kx * ky
  step <- xwidth/(kx - 1)
  r <- 0
  X <- rep(0, times = kx)
  Y2 <- rep(0, times = ky)
  for (p in 1:kx) {
    ystart <- ybegin
    xstart <- xstart + step
    for (q in 1:ky) {
      ystart <- ystart + step
      r <- r + 1
      X[r] <- xstart
      Y2[r] <- ystart
    }
  }
  TPS <- bendingenergy(TT)
  gamma11 <- TPS$gamma11
  gamma21 <- TPS$gamma21
  gamma31 <- TPS$gamma31
  W <- gamma11 %*% YY
  ta <- t(gamma21 %*% YY)
  B <- gamma31 %*% YY
  WtY <- t(W) %*% YY
  trace <- c(0)
  for (i in 1:m) {
    trace <- trace + WtY[i, i]
  }
  if(m==2){benergy <- (16 * pi) * trace}else{benergy <- (8 * pi) * trace}
  l <- kx * ky
  phi <- matrix(0, l, m)
  s <- matrix(0, k, 1)
  
  for (islice in 1:length(zslice)) {
    if (m == 3) {
      refc <- matrix(c(X, Y2, rep(zpos[islice], times = kx * 
                                    ky)), kx * ky, m)
    }
    if (m == 2) {
      refc <- matrix(c(X, Y2), kx * ky, m)
    }
    for (i in 1:l) {
      s <- matrix(0, k, 1)
      for (im in 1:k) {
        s[im, ] <- shapes::sigma(refc[i, ] - TT[im, ])
      }
      phi[i, ] <- ta + t(B) %*% refc[i, ] + t(W) %*% s
    }
    if(graphics==T){
    if (m == 3) {
      if (opt == 2) {
        shapes3d(TT, color = 2, axes3 = axes3, rglopen = FALSE)
        shapes3d(YY, color = 4, rglopen = FALSE)
        for (i in 1:k) {
          lines3d(rbind(TT[i, ], YY[i, ]), col = 1)
        }
        for (j in 1:kx) {
          lines3d(refc[((j - 1) * ky + 1):(ky * j), ], 
                  color = 6)
        }
        for (j in 1:ky) {
          lines3d(refc[(0:(kx - 1) * ky) + j, ], color = 6)
        }
      }
      shapes3d(TT, color = collandsTT, axes3 = axes3, rglopen = FALSE)
      shapes3d(YY, color = collandsYY, rglopen = FALSE)
      for (i in 1:k) {
        lines3d(rbind(TT[i, ], YY[i, ]), col = colshift)
      }
      for (j in 1:kx) {
        lines3d(phi[((j - 1) * ky + 1):(ky * j), ], color = colgrid)
      }
      for (j in 1:ky) {
        lines3d(phi[(0:(kx - 1) * ky) + j, ], color = colgrid)
      }
    }
    
  }
    
    
  }
  
  if (m == 2) {
    par(pty = "s")
    if (opt == 2) {
      par(mfrow = c(1, 2))
      order <- linegrid(refc, kx, ky)
      
      if(is.null(xlim)==T){xlim = c(xbegin -xwidth * ext, xbegin + xwidth * (1 + ext))}else{xlim=xlim}
      if(is.null(ylim)==T){ylim = c(ybegin - (xwidth * ky)/kx * ext, ybegin + (xwidth * ky)/kx * (1 + ext))}else{ylim=ylim}
      if(graphics==T){
        plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = xlim, ylim = ylim, xlab = " ", ylab = " ",asp=asp,,axes=axes2d)
        lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l), 2], type = "l",col=colgrid,xlim = xlim, ylim = ylim,asp=asp)
        points(TT, cex = cex, pch = pch, col = collandsTT)
      }}
    if(is.null(xlim)==T){xlim = c(xbegin -xwidth * ext, xbegin + xwidth * (1 + ext))}else{xlim=xlim}
    if(is.null(ylim)==T){ylim = c(ybegin - (xwidth * ky)/kx * ext, ybegin + (xwidth * ky)/kx * (1 + ext))}else{ylim=ylim}
    order <- linegrid(phi, kx, ky)
    if(graphics==T){plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = xlim, ylim = ylim, xlab = " ", ylab = " ",col=colgrid,asp=asp,axes=axes2d)
      lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l), 2], type = "l",col=colgrid,xlim = xlim, ylim = ylim,asp=asp)}
    if(graphics==T){points(YY, cex = cex, pch = pch, col = collandsYY)}
    if(graphics==T){points(TT, cex = cex, pch = pch, col = collandsTT)}
    if(graphics==T){
      if(displ==T){
        for (i in 1:(k)) {
          arrows(TT[i, 1], TT[i, 2], YY[i, 1], YY[i, 2], col = coldispl, length = 0.1, angle = 20)
        }
      }}
    firstcol<-order[1:l,][1:kx,][-kx,]
    firstrow<-order[(l + 1):(2 * l),][((kx*ky-1)-ky):(kx*ky-1),][-c(1:2),]
    lastcol<-order[1:l,][((kx*ky)-ky):(kx*ky),][-1,]
    lastrow<-order[(l + 1):(2 * l),][2:ky,][order((ky:2)),]
    bound = rbind(firstcol,firstrow,lastcol,lastrow)
    them<-nrow(order[1:l,])/ngrid
  }
  if(m==2){grid<-list(ngrid=order[1:l,],ngrid2=order,m=them,n=ngrid,bound=bound,l=l,xlim=xlim,ylim=ylim,TT=TT,YY=YY)}else{
    grid<-list(n=ngrid,l=l,TT=TT,YY=YY)
  }
  out<-list(YY=YY,gamma11=gamma11,gamma21=gamma21,gamma31=gamma31,W=W,ta=ta,B=B,WtY=WtY,trace=trace,benergy=benergy,grid=grid,kx=kx,ky=ky)
  if(graphics==T){if(!is.null(linksTT)){lineplot(TT,linksTT,lwd=lwdtt,col=collinksTT)}}
  if(graphics==T){if(!is.null(linksYY)){lineplot(YY,linksYY,lwd=lwdyy,col=collinksYY)}}
  return(out)
}
#' @export
