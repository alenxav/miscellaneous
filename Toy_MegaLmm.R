
# Get some data
data(tpod,package = 'bWGR')
Z = gen-1
set.seed(123)
h0 = bWGR::SimY(Z,50,h2 = 0.2)
Y = h0$Y
Y[sample(length(Y),length(Y)*0.2)] = NA

# MegaLMM function
MegaLmm = function(Y,Z,...){
  require(bWGR)
  k = ncol(Y)
  cat('Step 1\n')
  pb = txtProgressBar(style = 3)
  UvBeta = sapply(1:k,function(i){
    y = Y[,i]
    w = which(!is.na(y))
    yy = y[w]
    xx = Z[w,]
    #beta = emML(yy,xx)
    beta = MRR3(matrix(yy),xx,...)
    setTxtProgressBar(pb, i/k)
    return( c(beta$b) )})
  close(pb)
  cat('Step 2\n')
  pb = txtProgressBar(style = 3)
  MvBeta = sapply(1:k,function(i){
    y = Y[,i]
    w = which(!is.na(y))
    yy = y[w]
    xx = Z[w,]
    zz = xx %*% UvBeta[,-i]
    v = apply(zz,2,sd)
    zz = apply(zz,2,scale)
    ww = cbind(zz,xx)
    #beta = emML(yy,ww)
    beta = MRR3(matrix(yy),ww,...)
    w = 1:(k-1)
    betaf = c(UvBeta[,w] %*% (c(beta$b)[w]*v)) + c(beta$b)[-w]
    setTxtProgressBar(pb, i/k)
    return(betaf)})
  close(pb)
  out = list(UvLMM=UvBeta,MegaLMM=MvBeta)
}

# Benchmark
system.time(fit <- MegaLmm(Y,Z,NLfactor=1))[3]
Hat = lapply(fit, function(x) Z%*%x )
Hat$MvViaGS = MRR3(Y,Z)$hat
acc = sapply(Hat, function(x) diag(cor(x,h0$tbv)) )
colMeans(acc)


