#' scoresequencearray
#'
#' This function computes a sequence of shapes computed by deformed a reference shape ("mshape" argument) according to the eigenvectors specified in "rotation" argument. The sequence is computed in the interval between "from" and "to" argument. Usually the "rotation" argument comes from a PCA but any ordination method could be used.
#' @param mshape numeric: a kx2 matrix which will be deformed according to "from" and "to" values and to "rotation" eigenvectors. 
#' @param from numeric: a starting value for the computation od the sequence of shapes
#' @param to numeric: the final value for the computation od the sequence of shapes
#' @param pcn numeric: column number of "rotation" argument to be used to deform "mshape"
#' @param rotation numeric: rotation
#' @param mag numeric: magnification (default=1)
#' @param frames numeric: number of frames to be saved (default=15)
#' @return myarray array: an kxmxn array of predited shapes
#' @author Paolo Piras
#' @examples
#' \dontrun{ 
#' data(my2d)
#' linksdors<-list(c(1,2),c(37,7),c(12,4),c(27,28),c(25,21),c(38,40),c(9,10),c(2,3),c(3,4),c(1,7),c(1,6),c(3,5),c(6,40),c(5,9),c(40,8),c(8,9),c(1,7),c(7,6),c(3,4),c(4,5),c(39,38),c(38,35),c(35,37),c(37,39),c(35,34),c(34,33),c(33,32),c(32,31),c(31,30),c(30,29),c(29,37),c(37,36),c(36,29),c(28,31),c(28,30),c(13,10),c(10,11),c(11,12),c(12,13),c(13,14),c(14,16),c(16,17),c(17,20),c(20,19),c(19,18),c(18,12),c(18,15),c(15,12),c(21,19),c(21,20),c(24,25),c(25,26),c(26,27),c(27,24),c(26,24),c(24,23),c(23,22),c(22,8),c(8,2))                                                                                                                                                           
#' my2d<-read.inn(as.matrix(read.table("C:/docu paolo ac/LAVORI/ontogeny package/data/macroscelidea.txt",header=F,row.names=1)),40,2)
#' amy2d<-procSym(my2d,CSinit=T,scale=F)
#' pc1seq<-scoresequencearray(amy2d$mshape,min(amy2d$PCscores[,1]),max(amy2d$PCscores[,1]),1,amy2d$PCs)
#' par(mfrow=c(4,4))
#' for(i in 1:dim(pc1seq)[3]){
#' plotmyarrays(pc1seq[,,i],links=linksdors,cex=0,txt=F)}
#' }
#' @export
scoresequencearray<-function(mshape,from,to,pcn,rotation,mag=1,frames=15,type=1){
  require(abind)
  require(Morpho)
  if(length(pcn)<2){
    seqvector<-seq(from=from,to=to,length.out=frames)
    myarray<-NULL
    for(i in 1:length(seqvector)){
      
      
      if(type==1){
        shapeposition<-showPC(seqvector[i]*mag,rotation[,pcn],mshape)}else{shapeposition<-showPC2(seqvector[i]*mag,rotation[,pcn],mshape)}
      
      
      
      shapearr<-array(shapeposition,dim=c(nrow(shapeposition),ncol(shapeposition),1))
      myarray<-abind::abind(myarray,shapearr)
    }}else{
      
      binda<-rbind(from,to)
      seqmatrix<-NULL
      for(i in 1:ncol(binda)){
        seqmatrixi<-seq(binda[1,i],binda[2,i],length.out=frames)
        seqmatrix<-cbind(seqmatrix,seqmatrixi)
      }
      myarray<-NULL
      for(i in 1:nrow(seqmatrix)){
        
        if(type==1){
          shapeposition<-showPC(seqmatrix[i,]*mag,rotation[,pcn],mshape)}else{shapeposition<-showPC(seqmatrix[i,]*mag,rotation[,pcn],mshape)}
        
        
        
        shapearr<-array(shapeposition,dim=c(nrow(shapeposition),ncol(shapeposition),1))
        myarray<-abind::abind(myarray,shapearr)
      }}
  return(myarray)
}
  
  
