#' Internal_functions
#'
#' Here are reported a collection of internal functions
#' @author Paolo Piras
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
  
  init<-fin+(init-fin)*mag
  
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

ordiwithshapes<-function(mshape,scores,rotation,whichaxes=c(1,2),addata=NULL,asp=1,xlab="PC1",ylab="PC2",triang=NULL,factraj=NULL,coltraj=NULL,distscale=1,from=0,to=NULL,mag=1,procSym=T,subplotdim=2,legendgroup=NULL,shiftposx=1.5,shiftnegx=1.5,shiftposy=1.5,shiftnegy=1.5,links=NULL,collinks=1,grouphull=NULL,colhull=NULL,labels=NULL,labelscex=1,pch=19,col=1,cex=1,mult=0.5,tit3d=T,plotsource=T,xlim=extendrange(extendrange(range(scores[,whichaxes[1]]),f=mult)),ylim=extendrange(extendrange(range(scores[,whichaxes[2]]),f=mult))){
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
      diffonmesh(mshape,pc1pos,triang,title=F,distscale=distscale,from=from,to=to,plotsource=plotsource)
      next3d()
      text3d(0,0,0,paste("Scores",whichaxes[1],"pos"))
      next3d()
      diffonmesh(mshape,pc2pos,triang,title=F,distscale=distscale,from=from,to=to,plotsource=plotsource)
      next3d()
      text3d(0,0,0,paste("Scores",whichaxes[2],"pos"))
      next3d() 
      diffonmesh(mshape,pc1neg,triang,title=F,distscale=distscale,from=from,to=to,plotsource=plotsource)
      next3d();text3d(0,0,0,paste("Scores",whichaxes[1],"neg"))
      next3d()
      diffonmesh(mshape,pc2neg,triang,title=F,distscale=distscale,from=from,to=to,plotsource=plotsource)
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

