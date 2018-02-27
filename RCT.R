# Random Categorical & Continuous Term
#
RCT = function(y,x,lmb){
  Mean = function(x,lmb) sum(x)/(length(x)+lmb)
  if(is.factor(x)){ b = tapply(y,x,Mean,lmb=lmb); e = y-b[x]
  }else{ b = crossprod(y,x)/(crossprod(x)+lmb); e = y-x*b }
  return(list(b=b,e=e))}
#
RCT_update = function(y,x,b,lmb){
  Mean = function(x,lmb) sum(x)/(length(x)+lmb)
  if(is.factor(x)){ m = tapply(y,x,Mean,lmb=lmb); e = y-m[x]
  }else{ m = crossprod(y,x)/(crossprod(x)+lmb); e = y-x*m }
  b=m+b; return(list(b=m,e=e))}
