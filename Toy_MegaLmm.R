# Preconditioning on FA
MEGA_PCFA = function(Y,Z,npc=3){
  mu0 = colMeans(Y,na.rm=T)
  sd0 = apply(Y,2,sd,na.rm=T)
  for(i in 1:ncol(Y)) Y[,i]= (Y[,i]-mu0[i])/sd0[i]
  Yh = bWGR:::ZSEMF(Y,Z)
  Y2 = Y; Y2[is.na(Y2)] = Yh$hat[is.na(Y2)]
  E = bWGR::EigenBDCSVD(Y2)
  FA = E$U[,1:npc] %*% diag(E$D[1:npc]) %*% t(E$V[,1:npc])
  FAh = bWGR:::ZSEMF(FA,Z)
  Jh = bWGR:::ZSEMF(Y-FA,Z)
  mu = FAh$mu + Jh$mu + mu0
  b = t(sd0 * t(FAh$b + Jh$b))
  gebv = Z %*% b
  hat = t(mu0 + sd0 * t(FA + Jh$hat))
  h2 = FAh$h2 + Jh$h2
  colnames(gebv) = names(h2) = names(mu) = colnames(Y)
  colnames(b) = colnames(Y); rownames(b) = colnames(Z)
  GC = cor(gebv)
  return(list(mu=mu, b=b, gebv=gebv, hat=hat, GC=GC, h2=h2))
}

# Simulate some data
Z = bWGR::SimZ()
GC = bWGR::SimGC()
S = bWGR::SimY(Z,h2 = 0.3, PercMiss = 0.5,GC=GC)
Y = S$Y

# Fit model, test accuracy
fit = MEGA_PCFA(Y,Z)
mean(diag(cor(fit$gebv,S$tbv)))
