require(bWGR)
data(tpod)
X1 = gen[1:100,]
X2 = gen[101:196,]

Acc2 = function(X1,X2,h2=0.4){
  alpha = 1/sqrt(mean(apply(X1,1,crossprod)))
  V = EigenBDCSVD(crossprod(X1)/alpha,T)$V
  Z1 = X1 %*% V *alpha
  Z2 = X2 %*% V *alpha
  D = apply(Z1,2,crossprod)
  ve = ((1-h2)/h2)
  VarBhat = 1-1/(D/ve+1)
  sqrtVarBhat = sqrt(VarBhat)
  Acc = apply(Z2,1,function(z) sqrt( crossprod(z*sqrtVarBhat)/(crossprod(z)) )  )
  return(Acc)}

h2=0.5
tmp = EigenAcc(X1,X2,h2)
a = Acc2(X1,X2,h2)
plot(a,tmp)
