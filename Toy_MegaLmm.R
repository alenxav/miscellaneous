# Preconditioning on FA imputed with GEBVs
MEGA_PCFA = function(Y,Z,npc=2*sqrt(ncol(Y))){
  mu0 = colMeans(Y,na.rm=T)
  sd0 = apply(Y,2,sd,na.rm=T)
  for(i in 1:ncol(Y)) Y[,i]= (Y[,i]-mu0[i])/sd0[i]
  Yh = bWGR:::ZSEMF(Y,Z)
  Y2 = Y; Y2[is.na(Y2)] = Yh$hat[is.na(Y2)]
  E = bWGR::EigenBDCSVD(Y2)
  FA = E$U[,1:npc] %*% diag(E$D[1:npc]) %*% t(E$V[,1:npc])
  FAh = bWGR:::ZSEMF(E$U[,1:npc],Z)
  Jh = bWGR:::ZSEMF(Y-FA,Z)
  b = t(sd0 * t( FAh$b %*% diag(E$D[1:npc]) %*% t(E$V[,1:npc])  + Jh$b))
  gebv = Z %*% b
  hat = t(mu0 + sd0 * t(FA + Jh$hat))
  h2 = sapply(1:ncol(Y),function(i) cor(Y[,i],gebv[,i],use='p')^2 )
  colnames(gebv) = names(h2) = names(mu0) = colnames(Y)
  colnames(b) = colnames(Y); rownames(b) = colnames(Z)
  return(list(mu=mu0, b=b, 
              gebv=gebv, hat=hat, 
              GC=cor(gebv), h2=h2))}

# Iterative FA via SVD
IFA_MEGA = function(Y,Z,...){
  FA_dec = function(Y,iter=10,npc=3,verb=FALSE,std=TRUE){
    if(std){
      mu0 = colMeans(Y,na.rm=T)
      sd0 = apply(Y,2,sd,na.rm=T)
      Y2 = apply(Y,2,scale)
      Y2[is.na(Y)] = 0  }else{
      Y2 = apply(Y,2,function(x){x[is.na(x)]=mean(x,na.rm=T);x})}
    for(i in 1:iter){
      miss0 = Y2[is.na(Y)]
      E = bWGR::EigenBDCSVD(Y2)
      Yhat = E$U[,1:npc] %*% diag(E$D[1:npc]) %*% t(E$V[,1:npc])
      Y2[is.na(Y)] = Yhat[is.na(Y)]
      miss1 = Y2[is.na(Y)]
      res = round(c(-log10(crossprod(miss0-miss1))),2)
      if(verb) cat('it',i,' = ',res,'\n') }
    if(std) for(i in 1:ncol(Y)) Y2[,i] = Y2[,i]*sd0[i]+mu0[i]
    return(bWGR::EigenBDCSVD(Y2))}
  npc = round(sqrt(ncol(Y)))
  test = FA_dec(Y,...)
  F0 = test$U[,1:npc] %*% diag(test$D[1:npc])
  F1 = FUVBETA(F0, Z )
  FA = F0 %*% t(test$V[,1:npc])
  J0 = Y - FA
  J1 = FUVBETA(J0, Z )
  B = F1 %*% t(test$V[,1:npc]) + J1
  hat = Z %*% B
  mu = colMeans(Y-hat,na.rm = T)
  for(i in 1:ncol(hat)) hat[,i]=hat[,i]+mu[i]
  return(list(hat=hat,b=B,mu=mu))}


# Simulate some data
Z = bWGR::SimZ()
GC = bWGR::SimGC()
S = bWGR::SimY(Z,h2 = 0.3, PercMiss = 0.5,GC=GC)
Y = S$Y

# Fit model, test accuracy
fit1 = MEGA_PCFA(Y,Z)
mean(diag(cor(fit1$gebv,S$tbv)))
              
fit2 = IFA_MEGA(Y,Z)
mean(diag(cor(fit2$hat,S$tbv)))
              
