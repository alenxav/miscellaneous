# Fixed Categorical & Continuous Term
#
FCT = function(y,x){
  if(is.factor(x)){ b = tapply(y,x,mean); e = y-b[x]
  }else{ b = crossprod(y,x)/crossprod(x); e = y-x*b }
  return(list(b=b,e=e))}
#
FCT_update = function(y,x,b){
  if(is.factor(x)){ m = tapply(y,x,mean); e = y-m[x]
  }else{ m = crossprod(y,x)/crossprod(x); e = y-x*m }
  b=m+b; return(list(b=m,e=e))}
