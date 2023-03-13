# Get some data
data(tpod,package = 'bWGR')
Z = gen-1
set.seed(123)
h0 = bWGR::SimY(Z,50,h2 = 0.2)
Y = h0$Y
Y[sample(length(Y),length(Y)*0.2)] = NA

# MegaLMM function
MegaLmm = function(Y,Z,NumLoadings=5,TOI=NULL){
  require(mas)
  require(RSpectra)
  k = ncol(Y)
  cat('Step 1\n')
  UvBeta = sapply(1:k,function(i){
    y = Y[,i]
    w = which(!is.na(y))
    yy = y[w]
    xx = Z[w,]
    beta = MLM(matrix(yy),X = matrix(1,length(yy)),Z=xx)
    return( c(beta$u) )})
  cat('Step 2\n')
  if(is.null(TOI)){ toi = 1:k  }else{ toi = TOI } 
  MvBeta = sapply(toi,function(i){
    y = Y[,i]
    w = which(!is.na(y))
    yy = y[w]
    xx = Z[w,]
    zz = xx %*% UvBeta
    E = svds(zz,min(NumLoadings,ncol(Y)))
    X = cbind(1,zz %*% E$v)
    beta = MLM(matrix(yy),X,xx)
    betaf = c( UvBeta %*% E$v %*% beta$b[-1,] ) + c(beta$u)
    return(betaf)})
  out = list(UvLMM=UvBeta,MegaLMM=MvBeta)
}

# Benchmark
system.time(fit <- MegaLmm(Y,Z,2))[3]
Hat = lapply(fit, function(x) Z%*%x )
system.time(Hat$MvViaGS <- MRR3(Y,Z)$hat)[3]
acc = sapply(Hat, function(x) diag(cor(x,h0$tbv)) )
colMeans(acc)

