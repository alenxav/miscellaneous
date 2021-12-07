
QuantGwas = function(y,Z,Tau=0.75){
  
  require(quantreg)
  require(NAM)
  
  # Find alpha that maximized REML
  fit0 = suppressWarnings(reml(y,Z=Z))
  alpha = fit0$VC[['Vg']]/fit0$VC[['Ve']]
  V = tcrossprod(Z)*alpha + diag(length(y))
  EVD = eigen(V,symmetric = T)
  W = EVD$vectors %*% diag(sqrt(EVD$values))
  rm(V,fit0,EVD)
  gc()
  
  # Absorb random effects
  Wy = c( t(W) %*% y )
  WZ = t(W) %*% Z

  fit_NullModel = suppressWarnings(rq(Wy~1,tau=Tau))
  
  # SMA function
  gwa_test = function(x){
    fit_AlternModel = rq(Wy~x,tau=Tau)
    effect = fit_AlternModel$coefficients[2]
    stat_test = suppressWarnings(anova(fit_NullModel,fit_AlternModel))
    pvalue = stat_test$table$pvalue
    return(c(pvalue=pvalue,effect=effect))
  }
  
  # Run GWAS
  gwas = data.frame(t(apply(WZ,2,gwa_test)))
  gwas[['logPval']] = -log10(gwas[,1])
  
  # Output
  return(gwas)
  
}


# Example
if(F){
  data(tpod,package = 'NAM')
  fit = QuantGwas(y,gen,Tau=0.6)
  plot(fit$logPval,xlab='SNP',type='h',pch=20,ylab='-log(p-value)',main='Quantile GWAS, Tau=0.6')
}


