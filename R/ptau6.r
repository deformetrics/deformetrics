#' ptau6
#'
#' This function perform the Euclidean Parallel transport for allometric investigations using the Linear Shift startegy described in Piras et al (2016). The Linear Shift is performed on MANCOVA predictions and original data are projected on PCA space computed on them.  
#' @param array numeric: an array kxmxn of landmark coordinates
#' @param factor character: variable factor that affiliates shapes to group levels 
#' @param CSinit logical: if TRUE shapes are scaled at unit size (default=TRUE)
#' @param sepur logical: if TRUE separate per-group multivariate regression between shape and size are performed on shapes aligned after separate GPAs  (default=FALSE)
#' @param polyn numeric: default=1 the degree of regression
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' library(Morpho)
#' library(gdata)
#' library(rgl)
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
#' objptau<-ptau6(dors4,factordors4,CSinit=T)
#' #  morphological apprciation is dramatically different!
#' plotptau5(objptau,linksdors,pch=as.numeric(objptau$factorord),col=rep(1,length(objptau$factorord)),round(max(centroid.size(objptau$arrayord)),digits=0)/2,shiftnegy=2,mag=2,subplotdim=1)
#' ## End(Not run)
#' }
#'
#' @export

ptau6<-function(array,factor,CSinit=T,sepure=F,polyn=1){
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
  print(manymultgr(depepure,indepepure,factor))
  print("Pairwise *linear* mancova p-values on original data") 
  print(pwpermancova(depepure,indepepure,factor)$p_adonis_pred1_pred2)
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
  
  prls<-lshift2(mypredictpure,factor,CSinit=CSinit,CR=procSym(mypredictpure[,,firstsfac(factor)],scale=F,pcAlign=F,reflect=F,CSinit=F)$mshape,locs=mypredictpure[,,firstsfac(factor)])
  
  common<-array2mat(procSym(newarray,pcAlign=,scale=F,CSinit=CSinit)$orpdata,n,k*m)
  print("Pairwise multivariate (linear) regression between shape and size" )
  print(pwpermancova(common,indepepure,factor)$p_adonis_pred1_pred2)
  space1<-prcomp(array2mat(prls$transported,n,k*m))
  space1mshape<-procSym(prls$transported,pcAlign=F)$mshape
  origtrasp<-prls$transported+myresidpure
  
  origproj<-predict(space1,array2mat(origtrasp,n,k*m))
  print("Pairwise multivariate (linear) regression between shape of transported data and size" )
  print(pwpermancova(array2mat(origtrasp,n,k*m),indepepure,factor)$p_adonis_pred1_pred2)
  depepure2<-origtrasp
  
  print("Individual multivariate (linear) regression between shape of transported data and size" )
  print(manymultgr(array2mat(depepure2,n,k*m),indepepure,factor))
  
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
  
  
  out<-list(arrayord=newarray,factorord=factor,depepure=depepure,indepepure=indepepure,thelmlistpure=thelmlistpure,predictpure=mypredictpure,residpure=myresidpure,shifted=array2mat(prls$transported,n,k*m),thelmlistpure2=thelmlistpure2,predictpure2=array2mat(prls$transported,n,k*m),predictpure3=mypredictpure2,residpure2=myresidpure2,space1=space1,origtrasp=array2mat(origtrasp,n,k*m),origproj=origproj,space1mshape=space1mshape)                                   
}
