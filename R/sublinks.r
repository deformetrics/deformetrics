#' sublinks
#' This function creates a new links file starting from a links file and a vector of lansmarks selection
#' @param links list: A list of vectors each with two elements indicating the landmarks connected by a link
#' @param subvec list:  A vector ordered in ascending order indicating the subselection of landmarks for which the new link file should be created
#' @return res list: A list with new links
#' @author Paolo Piras
#' @examples 
#' \ dont run{
#' data(linksentire)
#' data(pri3d)
#' newlinks<-sublinks(linksentire,c(1:17))
#' plotmyarrays(pri3d[1:17,,1],links=newlinks)
#' }
#' @export
sublinks<-function(links,subvec){
  warning("landmarks subselection in <subvec> argument must be ordered in increasing order; do the same for landmarks subselection in shape data")
  vec<-sort(unique(unlist(links)))
  appaia<-function(vec,subvec){
    subvecord<-sort(subvec)
    subvecext1<-rep(NA,length(vec))
    vec%in%subvecord
    subvecext1[vec%in%subvecord]<-subvecord
    res<-cbind(vec,subvecext1,as.numeric(factor(subvecext1)))
    res
  }
  appamatrix<-appaia(vec,subvec)
  matinit<-list2matrix(links) 
  col1<-matinit[,1]%in%appamatrix[,2]
  col2<-matinit[,2]%in%appamatrix[,2]
  colmatr<-cbind(col1,col2,matinit)
  matinitpul<-matinit[apply(colmatr[,1:2],1,sum)==2,]
  appa2<-na.omit(appamatrix) 
  
  for(i in 1:length(matinitpul[,1])){
    for(j in 1:nrow(appa2)){
      if(matinitpul[i,1]==appa2[j,2]){matinitpul[i,1]<-appa2[j,3]} 
      if(matinitpul[i,2]==appa2[j,2]){matinitpul[i,2]<-appa2[j,3]} 
    }
  }
  res<-array2list(matrix2arrayasis(matinitpul,1))
  res
}
