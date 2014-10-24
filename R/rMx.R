
rMx <-
function(Input){
  if(is.vector(Input)){Output<-t(as.matrix(Input))}  
  if(!is.vector(Input)){Output<-as.matrix(Input)}
  Output
}
