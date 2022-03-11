
require(bWGR)
data(tpod)

X = CNT(gen)
Y = SimY(X)$Y

OUT = TIME = c()
dn = list(paste0('ENV',1:ncol(Y)),c('Unstructured','CmpSymmetry','FactorAnalytics'))

for(i in 1:5){
  cat('Cross-validation round',i,'\n')
  w = sample(nrow(Y),nrow(Y)*0.2)
  Y0 = Y; Y0[w,]=NA
  time = c(
    f0=(system.time(f0 <- MRR3(Y0,X))[3]),
    f1=(system.time(f1 <- MRR3(Y0,X,HCS=T))[3]),
    f2=(system.time(f2 <- MRR3(Y0,X,XFA2=T))[3])
  )
  out = c(f0=diag(cor(f0$hat[w,],Y[w,])),
          f1=diag(cor(f1$hat[w,],Y[w,])),
          f2=diag(cor(f2$hat[w,],Y[w,])))
  OUT = rbind(OUT,out)
  TIME = rbind(TIME,time)
  print(matrix(colMeans(OUT),ncol=3,dimnames = dn))
}

colMeans(matrix(colMeans(OUT),ncol=3,dimnames = dn))

