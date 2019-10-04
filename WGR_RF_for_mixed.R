
rf = function(y,X){
  rfs = ranger::ranger(y~.,data.frame(y=c(y),X=X))
  hat = rfs$predictions
  out = list( b=rfs, hat=hat )
  return(out)
}
