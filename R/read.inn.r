#' read.inn
#'
#' This function convert a matrix in an array
#' @param matrix numeric: an array kxmxn of landmark coordinates
#' @param k numeric: number of landmarks
#' @param mtype numeric: if 0, if 1
#' @author Paolo Piras
#' @export  
read.inn<-function (matrix, k, m,type=0){
  ### "k"  equals the number of landmarks, while "m" the number of dimensions
  if(!is.matrix(matrix)){stop("read.inn requires a matrix")}
  if(!is.null(rownames(matrix))){rown<-rownames(matrix)}else{rown<-c(1:nrow(matrix))}
  if(type==0){
    vectorr<-as.vector(t(matrix))
    tem <- vectorr
    n <- length(tem)/(k * m)
    tem <- array(tem, c(m, k, n))
    tem <- aperm(tem, c(2, 1, 3))
    xup <- tem
    dimnames(xup)[3]<-list(rown)
  }else{NULL}
  
  if(type==1){
    
    finalue<-NULL
    for(i in 1:nrow(matrix)){      
      finali<-matrix(matrix[i,],k,m)
      finalue<-rbind(finalue,finali)   
    }
    finalue<-matrix2arrayasis(finalue,k)
    finalue
    dimnames(finalue)[3]<-list(rown)
    
  }else{NULL}
  
  
  if(exists("xup")==T){xup}else{NULL}
  if(exists("finalue")==T){finalue}else{xup}
  
}

