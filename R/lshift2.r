#' lshift2
#'
#' This function performs the Euclidean Parallel Transport, i.e. the Linear Shift strategy, described in Piras et al (2016). It applies a series of deformations occurring in different groups to a Common Reference by properly managing rotations. 
#' @param array numeric: an array kxmxn of landmark coordinates
#' @param factor character vector: a factor with group affiliation
#' @param CR logical: if NULLL CR is estimated as the consensus of the array via Generalized Procrustes Analysis (GPA). Otherwise the specified CR is used.
#' @param locs=NULL: if NULL the local references are estimated via separate per-group GPAs. Oterwise the specified locs are used.
#' @param CSinit=F logical: if TRUE scaling at unit size is performed.
#' @param scale=F logical “scale” argument in procSym() from package “Morpho”
#' @param reflect logical:  “reflect” argument in procSym() from package “Morpho”
#' @param reorder logical:  If TRUE the factor is reordered
#' @param sc logical: This is not supposed to be called by common user  
#' @return transported list: a list of length 1 containing the array with the transported shapes. 
#' @author Paolo Piras
#' @references Piras P., Teresi L., Traversetti L., Varano V, Gabriele S., Kotsakis T., Raia P., Puddu P.E., Scalici M. (2016). The conceptual framework of ontogenetic trajectories: Parallel Transport allows the recognition and visualization of pure deformation patterns. Evolution and Development 18: 182-200. doi: 10.1111/ede.12186
#' @examples
#' \dontrun{ 
#' data(case1)
#' group<-factor(rep(1:5,each=21))
#' linksdors<-conslinks(21)
#' library(gdata)
#' lscase1<-lshift2(case1,group,CR=procSym(case1[,,firstsfac(group)],CSinit=F,reflect=F,scale=F,pcAlign=F)$mshape,locs=case1[,,firstsfac(group)])
#' plottraj(procSym(lscase1$transported,CSinit=T)$PCscores[,1:2],group,asp=1) #### very similar to Levi Civita
#' }  
#' @export  


lshift2<-function(array,factor,CR=NULL,locs=NULL,CSinit=F,scale=F,reflect=F,reorder=F,sc=F){
  warning("individuals belonging to each factor must be consecutive")
  if(reorder==T){factor<-factor(factor,levels=unique(factor))}
  
  k<-dim(array)[1]
  m<-dim(array)[2]
  n<-dim(array)[3]
  ng<-table(factor)
  
  
  if(is.null(locs)==T){
    
    sepas<-sepgpa(array,factor,CSinit=CSinit,scale=scale,reflect=reflect)
    
    locs<-list2array(subListExtract(sepas$listofgpas,"mshape"))
    
  }else{locs<-locs}
  
  if(is.null(CR)==T){CR<-procSym(locs,CSinit=CSinit,scale=scale,reflect=reflect,pcAlign=F)$mshape}else{CR<-CR}
  
  if(CSinit==T){
    CR<-scaleshapes(CR)[,,1]
    locs<-scaleshapes(locs)
    array<-scaleshapes(array)
  }
  
  locsopas<-opaloop2(CR,locs,reflect=F)
  locsop<-locsopas$looped
  
  
  locsrot<-locsopas$rots
  specop<-NULL
  specrots<-NULL
  for(i in 1:nlevels(factor)){
    print(i)
    specopas<-opaloop2(locsop[,,i],array[,,as.numeric(factor)==i])
    specopi<-specopas$looped
    specop<-abind::abind(specop,specopi)
    specrotsi<-specopas$rots
    specrots<-c(specrots,list(specrotsi))
  }
  
  specrots2<-array(unlist(specrots),dim=c(m,m,n))
  
  specdiff2provv<-rep(array2list(locsop),ng)
  specdiff2<-specop-array(unlist( specdiff2provv), dim = c(k, m, n))
  
 
if(sc==T){
specdiff2<-multiplarray(specdiff2,c(1/rep(apply(locs,3,cSize),ng)))
valeu1<-(array(unlist(rep(list(CR),sum(ng))),dim=c(k,m,n)))+multiplarray(specdiff2,rep(cSize(CR),n))
}else{
  valeu1<-(array(unlist(rep(list(CR),sum(ng))),dim=c(k,m,n)))+specdiff2}
 out<-list("transported"=valeu1)
}
