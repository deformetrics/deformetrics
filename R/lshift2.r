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
#' data(my2d)
#' data(macrogroup)
#' linksdors<-list(c(1,2),c(37,7),c(12,4),c(27,28),c(25,21),c(38,40),c(9,10),c(2,3),c(3,4),c(1,7),c(1,6),c(3,5),c(6,40),c(5,9),c(40,8),c(8,9),c(1,7),c(7,6),c(3,4),c(4,5),c(39,38),c(38,35),c(35,37),c(37,39),c(35,34),c(34,33),c(33,32),c(32,31),c(31,30),c(30,29),c(29,37),c(37,36),c(36,29),c(28,31),c(28,30),c(13,10),c(10,11),c(11,12),c(12,13),c(13,14),c(14,16),c(16,17),c(17,20),c(20,19),c(19,18),c(18,12),c(18,15),c(15,12),c(21,19),c(21,20),c(24,25),c(25,26),c(26,27),c(27,24),c(26,24),c(24,23),c(23,22),c(22,8),c(8,2))                                                                                                                                                           
#' library(gdata)
#' mysel<-macrogroup%in%c("Macroscelides_proboscideus","Petrodromus_tetradactylus","Elephantulus_rozeti","Elephantulus_edwardii")
#' dors4<-my2d[,,mysel]
#' factordors4<-drop.levels(macrogroup[mysel],reorder=T)
#' factordors4<-factor(factordors4,levels=unique(factordors4))
#' adors4<-procSym(dors4,scale=F,pcAlign=F,reflect=F)
#' mypredictbook<-read.inn(predict(lm(array2mat(adors4$orpdata,105,80)~poly(adors4$size,1,raw=T)*factordors4)),40,2)### PCA on predictions from separate per-group regressions between shape and size
#' procbook<-procSym(mypredictbook,pcAlign=F)
#' plot(procbook$PCscores[,1:2],pch=as.numeric(factordors4),asp=1)
#' ordiwithshapes(procbook$mshape,procbook$PCscores,procbook$PCs,links=linksdors,subplotdim=1,pch=as.numeric(macrogroup))
#' objptau<-ptau6(dors4,factordors4,CSinit=T)### linear shift in the shape space
#' ##plot results from objptau; there you can find all necessary objects for plotting
#' ########  morphological appreciation is dramatically different!
#' plotptau5(objptau,linksdors,pch=as.numeric(objptau$factorord),col=rep(1,length(objptau$factorord)),round(max(centroid.size(objptau$arrayord)),digits=0)/2,shiftnegy=2,mag=2,subplotdim=1)
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
