#' pt2dvar
#'
#' This function performs the Levi Civita Riemannian Parallel Transport according to explicit formulation presented in Varano et al (accepted) [Varano V, Gabriele S,  Teresi L, Dryden I L,  Puddu P E, Torromeo C,  Piras P. Accepted. The TPS Direct Transport: a new method for transporting deformations in the Size-and-shape Space. International Journal of Computer Vision]. 
#' @param array	array: a kx2xn array of shapes 
#' @param group	character:	group structure
#' @param doopa logical: if TRUE Ordinary Procrustes Analysis is performed on the workhorse function. It is not supposed to be set at TRUE as OPAs are already performed outside that step. 
#' @param tol numeric: tolerance parameter for rotation
#' @param CR matrix:common Reference to which intra-group deformations are applied
#' @param locs array:	Local references (one for each level of group) from which intra-group deformations are computed
#' @param sss logical: if TRUE analyses are performed on the Size and Shape Space. It is not supposed to be set at FALSE
#' @return final list: a list with transported shapes
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' data(case1)
#' groupsimul<-factor(rep(1:5,each=21))
#' plottraj(procSym(case1,CSinit=T)$PCscores[,1:2],group,asp=1)######## a classic approach does not recover the cycle
#' ptcase1<-pt2dvar(case1,groupsimul,CR=procSym(case1[,,firstsfac(groupsimul)],CSinit=F,reflect=F,scale=F,pcAlign=F)$mshape,locs=case1[,,firstsfac(groupsimul)])
#' plottraj(procSym(ptcase1$final,CSinit=T)$PCscores[,1:2],groupsimul,asp=1) ####  with this kind Parallel Transport the cycle is not perfectly recovered as with the new approach of Direct Transport (Varano et al [accepted]). 
#' }
#' @export

pt2dvar<-function(array,group,doopa=F,tol=0.000001,CR=NULL,locs=NULL,sss=T){
    k<-dim(array)[1]
    m<-dim(array)[2]
    n<-length(group)
    ng<-table(group)
    if(is.null(locs)==T){locs<-procSym(meanarray(array,group,reorder=F),CSinit=F,scale=F,pcAlign=F,reflect=F)$rotated}else{locs<-locs}
    if(is.null(CR)==T){CR<-procSym(locs,CSinit=F,scale=F,pcAlign=F,reflect=F)$mshape}else{CR<-CR}
    
    locsopas<-opaloop(CR,locs,reflect=F)
    locsop<-locsopas$looped
    par(new=F)
  

    
    locsrot<-locsopas$rots
    specop<-NULL
    specrots<-NULL
    for(i in 1:nlevels(group)){
      specopas<-opaloop(locsop[,,i],array[,,as.numeric(group)==i])
      specopi<-specopas$looped
      specop<-abind::abind(specop,specopi)
      specrotsi<-specopas$rots
      specrots<-c(specrots,list(specrotsi))
    }
    specrots2<-list2array(unlist(specrots,recursive=F))
    
    
    specdiff<-NULL
    for(i in 1:nlevels(group)){
      for(j in 1:(ng[i])){
        specdiffi<-array(specop[,,as.numeric(group)==i,drop=F][,,j]-locsop[,,i],dim=c(k,m,1))
        specdiff<-abind::abind(specdiff,specdiffi)
      }}
    
    vbs<-NULL
    epsilons<-NULL
    epsilons2<-NULL
    trasprotva<-NULL
    trasprotmua<-NULL
    for(i in 1:nlevels(group)){
      for(j in 1:(ng[i])){ 
        
        pt2d<-pt2dvaruno(locsop[,,i],CR,specdiff[,,as.numeric(group)==i][,,j],doopa=doopa,tol=tol,sss=sss)
        vbsij<-array(pt2d$vb,dim=c(k,m,1))
        vbs<-abind::abind(vbs,vbsij)
        epsilonsij<-pt2d$eps
        epsilons2ij<-pt2d$eps2
        epsilons<-append(epsilons,epsilonsij)
        epsilons2<-append(epsilons2,epsilons2ij)
        trasprotmuaij<-pt2d$rotmua
        trasprotmua<-c(trasprotmua,list(trasprotmuaij))
        trasprotvaij<-pt2d$rotva
        trasprotva<-c(trasprotva,list(trasprotvaij))
      }}
    if(doopa==T){trasprotmua2<-list2array(trasprotmua)}else{trasprotmua2<-NULL}
    if(doopa==T){trasprotva2<-list2array(trasprotva)}else{trasprotva2<-NULL}
    final<-reparray(array(CR,dim=c(k,m,1)),n)+vbs
if(sss==F){final<-scaleshapes(final)}
    out<-list(final=final,trasprotsmua=trasprotmua2,trasprotsva=trasprotva2,epsilons=epsilons,epsilons2=epsilons2,locsop=locsop,locsrot=locsrot,specop=specop,specrots=specrots2)
    return(out)
  }
  