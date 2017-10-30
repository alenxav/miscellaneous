# Structured Random Term
#
SRT = function(y,x,gen,msx,Ve,model=c('NOR','DER'),w=FALSE,...){
  # y = vector of response variable
  # x = factor indicating ID
  # gen = genotypic matrix
  # w = weighted (LOGICAL)
  # msx = output of MSX
  # Ve = Residual variance
  # model = either Normal ('NOR') ou Laplace ('DER')
  #
  # Calculate weights
  if(w) W = tapply(x,x,length) else W = rep(1,length(levels(x)))
  # Compress and fit
  if(w){
    # weigheted
    PHE = tapply(y,x,sum)[rownames(gen)]
  }else{
    # unweigheted
    PHE = tapply(y,x,mean)[rownames(gen)]
  }
  if(model=='NOR') FIT = NOR2(PHE,gen,msx$cxx,Ve,msx$xx,...)
  if(model=='DER') FIT = DER2(PHE,gen,msx$cxx,Ve,msx$xx,...)
  # Decompress
  hat = FIT$h/W; names(hat) = rownames(gen)
  HAT = hat[x]
  E = y-HAT
  # Return
  RET = list(B=FIT$b,Va=FIT$v,fit=hat,Hat=HAT,E=E)
  return(RET)
}
