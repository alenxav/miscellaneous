
# Get some data
data(tpod,package = 'bWGR')
Z = gen-1
set.seed(123)
h0 = SimY(Z,10,h2 = 0.2)
Y = h0$Y
Y[sample(length(Y),length(Y)*0.2)] = NA

# MegaLMM function
MegaLmm = function(Y,Z){
  require(bWGR)
  k = ncol(Y)
  cat('Step 1\n')
  pb = txtProgressBar(style = 3)
  UvBeta = sapply(1:k,function(i){
    y = Y[,i]
    w = which(!is.na(y))
    yy = y[w]
    xx = apply(Z[w,],2,function(x)x-mean(x))
    beta = emML(yy,xx)
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
    zz = xx %*% UvBeta[,-k]
    beta = emML2(yy,xx,zz)
    betaf = c(UvBeta[,-k] %*% beta$b2) + beta$b1
    setTxtProgressBar(pb, i/k)
    return(betaf)})
  close(pb)
  out = list(UvLMM=UvBeta,MegaLMM=MvBeta)
}

# Benchmark
fit = MegaLmm(Y,Z)
Hat = lapply(fit, function(x) Z%*%x )
Hat$MvViaGS = mrr(Y,Z)$hat
acc = sapply(Hat, function(x) diag(cor(x,h0$tbv)) )
colMeans(acc)


