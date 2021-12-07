
# GWAS function
QuantGwas = function(y,Z,Tau=0.75){
  
  require(quantreg)
  require(NAM)
  
  # Centralize markers
  Z = apply(Z,2,function(x)x-mean(x))
  
  # Find alpha that maximized REML
  fit0 = suppressWarnings(reml(y,Z=Z))
  alpha = fit0$VC[['Vg']]/fit0$VC[['Ve']]
  V = tcrossprod(Z)*alpha + diag(length(y))
  EVD = eigen(V,symmetric = T)
  W = EVD$vectors %*% diag(sqrt(EVD$values))
  rm(V,fit0,EVD)
  gc()
  
  # Absorb random effects
  X = matrix(1,length(y))
  Wy = c( t(W) %*% y )
  WX = t(W) %*% X
  WZ = t(W) %*% Z
  
  # Fit null statistical model (no marker)
  fit_NullModel = suppressWarnings(rq(Wy~WX-1,tau=Tau))
  
  # Run GWAS
  GWAS = matrix(NA,nrow=ncol(Z),ncol=2,
                dimnames=list(colnames(Z),c('pvalue','effect')))
  pb = txtProgressBar(style = 3)
  for(j in 1:ncol(Z)){
    x = WZ[,j]
    fit_AlternModel = rq(Wy~WX+x-1,tau=Tau)
    effect = fit_AlternModel$coefficients[2]
    stat_test = suppressWarnings(anova(fit_NullModel,fit_AlternModel,se="iid"))
    pvalue = stat_test$table$pvalue
    GWAS[j,] = c(pvalue,effect)
    setTxtProgressBar(pb, j/ncol(Z))
  }
  close(pb)
  GWAS = data.frame(GWAS)
  GWAS[['logPval']] = -log10(GWAS$pvalue)
  
  # Output
  return(GWAS)
  
}


# Example
if(F){
  data(tpod,package = 'NAM')
  fit = QuantGwas(y,gen,Tau=0.6)
  plot(fit$logPval,xlab='SNP',type='h',pch=20,ylab='-log(p-value)',main='Quantile GWAS, Tau=0.6')
}


