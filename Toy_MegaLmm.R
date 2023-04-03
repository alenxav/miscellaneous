# Get some data and simulate phenotypes
data(tpod,package = 'bWGR')
Z = gen-1
set.seed(123)
h0 = bWGR::SimY(Z, k = 50, h2 = 0.2, GC = 0.5)
Y = h0$Y
Y[sample(length(Y),length(Y)*0.2)] = NA

# MegaLMM function
MegaLmm = function(Y,Z,TOI=NULL,...){
  require(bWGR)
  k = ncol(Y)
  cat('Step 1\n')
  pb = txtProgressBar(style = 3)
  Mu = colMeans(Y,na.rm=T)
  UvBeta = sapply(1:k,function(i){
    y = Y[,i]
    w = which(!is.na(y))
    yy = y[w]
    xx = Z[w,]
    beta = MRR3(matrix(yy),xx,...)
    setTxtProgressBar(pb, i/k)
    return( c(beta$b) )})
  close(pb)
  cat('Step 2\n')
  pb = txtProgressBar(style = 3)
  if(is.null(TOI)){ toi = 1:k  }else{ toi = TOI } 
  kk = length(toi)
  MvBeta = sapply(toi,function(i){
    y = Y[,i]
    w = which(!is.na(y))
    yy = y[w]
    xx = Z[w,]
    zz = xx %*% UvBeta
    v = apply(zz,2,sd)
    zz = apply(zz,2,scale)
    ww = cbind(zz,xx)
    beta = MRR3(matrix(yy),ww,...)
    q = 1:k
    betaf = c(UvBeta[,q] %*% (c(beta$b)[q]*v)) + c(beta$b)[-q]
    setTxtProgressBar(pb, i/kk)
    return(betaf)})
  close(pb)
  # Prepare output
  rownames(UvBeta) = rownames(MvBeta) = colnames(Z)
  colnames(UvBeta) = colnames(Y)
  if(is.null(TOI)){ colnames(MvBeta) = colnames(Y)  }else{ colnames(MvBeta) = colnames(Y)[TOI] } 
  out = list(UvLMM=UvBeta,MegaLMM=MvBeta,Mu=Mu)
}

# Benchmark
system.time(fit <- MegaLmm(Y,Z,NLfactor=1))[3]
Hat = lapply(fit[1:2], function(x) Z%*%x )
Hat$MvViaGS = MRR3(Y,Z)$hat
acc = sapply(Hat, function(x) diag(cor(x,h0$tbv)) )
round(colMeans(acc),2)
