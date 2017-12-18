#' newdt
#'
#' This function performs the TPS Parallel Transport according to Varano et al (accepted). Deformation occurring within groups are estimate relatively to per-group local references ("locs" argument") and applied to a Common Reference ("CR"). 
#' @param array numeric: an array kxmxn of landmark coordinates
#' @param group character vector: a factor with group affiliation
#' @param CR=NULL: if NULLL CR is estimated as the consensus of the array via Generalized Procrustes Analysis (GPA). Otherwise the specified CR is used.
#' @param pole=CR; by default the Riemannian pole coincides with CR. 
#' @param diffpole=F logical: if TRUE the pole is not coincident with CR. This is not supposed to be used by common users. 
#' @param locs=NULL: if NULL the local references ar estimated via separate per-group GPAs. Oterwise the specifiesed locs are used.
#' @param center=T logical: if TRUE shapes are centered. 
#' @param CSinit=F logical: if TRUE caling at unit size is performed. 
#' @param tolrot=1 numeric: tolerance for rotation during GPA. 
#' @param tol=1e-8 numeric: tolerance for eigenvalue decomposition
#' @param doopa=T logical: if TRUE locs are aligned with CR via Ordinary Procrustes Analyisis (OPA)
#' @param domopa=T logical: if TRUE each shape within each group is aligned with its proper loc via Modified Ordinary Procrustes Analysis introduced in Varano et al (accepted)
#' @return out matrix: a matrix containing the transported shapes. 
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' data(case1)
#' group<-factor(rep(1:5,each=21))
#' plottraj(procSym(case1,CSinit=T)$PCscores[,1:2],group,asp=1)######## a classic approach does not recover the cycle
#' dtcase1<-newdt(case1,group,CR=procSym(case1[,,firstsfac(group)],CSinit=F,reflect=F,scale=F,pcAlign=F)$mshape,locs=case1[,,firstsfac(group)],tolrot=10)
#' plottraj(procSym(dtcase1,CSinit=T)$PCscores[,1:2],group,asp=1)######## transported data recover the cycle
#' }
#' @export  
newdt<-function(array,group,CR=NULL,pole=CR,diffpole=F,locs=NULL,center=T,CSinit=F,tolrot=1,tol=1e-8,doopa=T,domopa=T,qid=F){
  if(center==T){array<-centershapes(array)}
  
  group<-factor(group,levels=unique(group))
  k<-dim(array)[1]
  m<-dim(array)[2]
  n<-length(group)
  ng<-table(group)
  if(is.null(CR)){CR<-procSym(array,scale=F,pcAlign=F,CSinit=CSinit)$mshape}else{CR=CR}
  
  if(is.null(locs)){
    locs<-NULL
    for(i in 1:nlevels(group)){
      locsi<-procSym(array[,,as.numeric(group)==i,drop=F],scale=F,CSinit=CSinit,pcAlign=F)$mshape
      locs<-abind::abind(locs,array(locsi,dim=c(k,m,1)))
    }}else(locs<-locs)
  
  if(doopa==T){
    loop<-opaloop2(CR,locs,reflect=F)
    pos1<-(acos(abs(list2array(loop$rots)[1,1,]))*180)/pi
    if(pos1=="NaN"){pos1<-tolrot+1e-4}
    if(max(pos1)>tolrot){locs<-loop$looped}
    print("Angles between CR and locs aligned with it")
    print(pos1)
    
    if(diffpole==T){
      loop<-opaloop2(pole,locs,reflect=F)
      pos1<-(acos(abs(list2array(loop$rots)[1,1,]))*180)/pi
      if(pos1=="NaN"){pos1<-tolrot+1e-4}
      if(max(pos1)>tolrot){locs<-loop$looped}
      print("Angles between CR and locs aligned with it")
      print(pos1)
    }}
  
  if(m>2){
    dets<-NULL
    for(i in 1:dim(locs)[[3]]){
      Atsi<-tpsdry2(locs[,,i],CR,doopa=F,meth="mor",g11=F)$at
      dets<-c(dets,list(det(Atsi)^((1/6))))
    }
  }else{dets<-replicate(dim(locs)[[3]], c(1), FALSE)}
  
  
  
  if(domopa==T){
    hiermopizedtps<-NULL
    hiermopized<-NULL
    for(i in 1:nlevels(group)){
      for(j in 1:(ng[i])){
        hiermopizedtpsij<-mopa(locs[,,i],array[,,as.numeric(group)==i,drop=F][,,j],rot=c("mopa"),CSinit=CSinit)
        hiermopizedtps<-c(hiermopizedtps,list(hiermopizedtpsij[5:7]))
        hiermopized<-abind::abind(hiermopized,hiermopizedtpsij$opizzata)
      }}
  }else{hiermopized<-array}
  
  
  dummygm<-newmb(CR,CR)
  eig11<-eigen(dummygm$gamma11)
  #datoglie<-which(Re(eig11$value)<tol)
  datoglie<-     (length(Re(eig11$value))-(m-1)):(length(Re(eig11$value)))
  Ub<-t(eig11$vectors)
  sqrts<-sqrt(eig11$values)
  sqrts[which(sqrts=="NaN")]<-0
  appo<-diag(sqrts)%*%Ub
  Ubmod<-Ub[-datoglie,]
  if(is.matrix(Ubmod)==F){Ubmod<-t(as.matrix(Ubmod))}
  if(diffpole==T){
    dummypole<-newmb(pole,pole)
    eigpole11<-eigen(dummypole$gamma11)
    datogliepole<-which(Re(eigpole11$value)<tol)
    datogliepole<-  (length(Re(eigpole11$value))-(m-1)):(length(Re(eigpole11$value)))
    Ubpole<-t(eigpole11$vectors)
    sqrteig11polval<-sqrt(eigpole11$values)
    sqrteig11polval[which( sqrteig11polval=="NaN")]<-0
    appopole<-diag(sqrteig11polval)%*%Ubpole
    Ubpolemod<-Ubpole[-datogliepole,]
    #Mb<-rbind(t(dummypole$h%*%CR),appo[-datoglie,]%*%dummypole$snew)#
    
    
    svdb<-svd(Ubpolemod%*%dummypole$snew%*%dummygm$snew%*%t(Ubmod))
    Rb<-svdb$v%*%t(svdb$u)
    Mb<-rbind(t(dummygm$h%*%CR),t(Rb)%*%appo[-datoglie,]%*%dummygm$snew)
  }else{
    
    Mb<-rbind(t(dummygm$h%*%CR),Re(appo)[-datoglie,]%*%dummygm$snew)
  }
  
  if(diffpole==F){
    transported<-NULL
    for(i in 1:nlevels(group)){
      for(j in 1:(ng[i])){
        vaij<-dummygm$h%*%(hiermopized[,,as.numeric(group)==i][,,j]-locs[,,i])
        newmbij<-newmb(locs[,,i],hiermopized[,,as.numeric(group)==i][,,j])
        eig11ij<-eigen(newmbij$gamma11)
        Uaij<-t(eig11ij$vectors)
        sqrteig11val<-sqrt(eig11ij$values)
        sqrteig11val[which(sqrteig11val=="NaN")]<-0
        appoij<-diag(sqrteig11val)%*%Uaij
        Uamodij<-Uaij[-datoglie,]
        
        
        if(is.matrix(Uamodij)==F){Uamodij<-t(as.matrix(Uamodij))}
        
        
        svdij<-svd(Ubmod%*%dummygm$snew%*%newmbij$snew%*%t(Uamodij))
        if(qid==F){Rij<-svdij$v%*%t(svdij$u)}else{Rij<-diag(k-m-1)}##############   qui è se vuoi usare l'identità
        Maij<-rbind(t(dummygm$h%*%locs[,,i]),t(Rij)%*%appoij[-datoglie,]%*%newmbij$snew)
        vbij<-(newmbij$h%*%CR%*%newmbij$gamma21+dets[[i]]*dummygm$snew%*%armaGinv(as.matrix(Mb))%*%Maij%*%newmbij$gamma11)%*%vaij
        transportedij<-CR+t(newmbij$h)%*%vbij
        transported<-abind::abind(transported,array(transportedij,dim=c(k,m,1)))
      }
    }
  }
  
  if(diffpole==T){
    transported<-NULL
    for(i in 1:nlevels(group)){
      for(j in 1:(ng[i])){
        vaij<-dummypole$h%*%(hiermopized[,,as.numeric(group)==i][,,j]-locs[,,i])
        newmbij<-newmb(locs[,,i],hiermopized[,,as.numeric(group)==i][,,j])
        eig11ij<-eigen(newmbij$gamma11)
        Uaij<-t(eig11ij$vectors)
        sqrteig11val<-sqrt(eig11ij$values)
        sqrteig11val[which(sqrteig11val=="NaN")]<-0
        appoij<-diag(sqrteig11val)%*%Uaij
        
        Uamodij<-Uaij[-datoglie,]
        if(is.matrix(Uamodij)==F){Uamodij<-t(as.matrix(Uamodij))}
        
        svdij<-svd(Ubpolemod%*%dummypole$snew%*%newmbij$snew%*%t(Uamodij))
        if(qid==F){Rij<-svdij$v%*%t(svdij$u)}else{Rij<-diag(k-m-1)}##############   qui è se vuoi usare l'identità
        Maij<-rbind(t(dummygm$h%*%locs[,,i]),t(Rij)%*%appoij[-datoglie,]%*%newmbij$snew)
        vbij<-(newmbij$h%*%CR%*%newmbij$gamma21+dets[[i]]*dummygm$snew%*%armaGinv(Mb)%*%Maij%*%newmbij$gamma11)%*%vaij
        transportedij<-CR+t(newmbij$h)%*%vbij
        transported<-abind::abind(transported,array(transportedij,dim=c(k,m,1)))
      }
    }
  }
  
  
  transported<-Re(transported)
  return(transported)
}
