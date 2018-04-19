tlm = function(formula,data,Tukey=TRUE){
  
  # fit the model
  fnew = update(formula,~.-1)
  fit = lm(formula=fnew,data=data)
  
  # recreate formula
  data_new = fit$model
  nl = length(levels(data_new[[2]]))
  fit$terms = formula(terms.formula(formula)[-1])
  fit$coefficients = fit$coefficients[-c(1:nl)]
  fit$effects = fit$effects[-c(1:nl)]
  fit$rank = fit$rank-nl
  fit$assign = fit$assign[-c(1:nl)]-1
  fit$xlevels = fit$xlevels[-1]
  if(length(fit$contrasts)>1){
    fit$contrasts = fit$contrasts[-1]
    }else{fit$contrasts=NULL}
  
  # Rearrange data with conditional phenotype
  suppressWarnings(hat <- predict(fit,data_new))
  y = data_new[as.character(formula[2])][,1]
  res = (y - hat);
  Y = res-mean(res)+mean(y)
  X = data_new[as.character(formula( terms.formula(formula)[1] ))[3]][,1]
  lvl = levels(X)
  
  # Refit and multiple testing
  refit = lm(Y~X-1,data_new)
  tk = TukeyHSD(aov(refit))
  tt = pairwise.t.test(Y,X,pool.sd=FALSE)
  
  # mean and sd
  mu = tapply(Y,X,mean)
  std = tapply(Y,X,sd)
  lwr = mu-std
  upr = mu+std
  COL = rainbow(nl)
  rnk = names(sort(mu,T))
  
  # t-test output
  ttX = matrix(NA,nl,nl,dimnames = list(rnk,rnk))
  for(i in lvl){for(j in lvl){
    if(i%in%rownames(tt$p.value) & j%in%colnames(tt$p.value)){
        ttX[i,j] = ttX[j,i] = tt$p.value[i,j] }}} 
  diag(ttX) = 1
  #
  ttrnk = rep(NA,nl); ttrnk[1] = 1
  for(i in 2:nl) ttrnk[i] = which(ttX[,i]>0.05)[1]
  for(i in 2:nl) if(ttrnk[i]>ttrnk[i-1]) ttrnk[i]=ttrnk[i-1]+1
  ttrnk = LETTERS[ttrnk]
  names(ttrnk) = rnk;
  
  # Tukey output
  tkX = matrix(NA,nl,nl,dimnames = list(rnk,rnk))
  for(i in lvl){for(j in lvl){
    if(paste(i,j,sep='-')%in%rownames(tk$X)){
    tkX[i,j] = tkX[j,i] = tk$X[which(paste(i,j,sep='-')==rownames(tk$X)),4]
    }}} 
  diag(tkX) = 1
  #
  tkrnk = rep(NA,nl)
  tkrnk[1] = 1
  for(i in 2:nl) tkrnk[i] = which(tkX[,i]>0.05)[1]
  for(i in 2:nl) if(tkrnk[i]>tkrnk[i-1]) tkrnk[i]=tkrnk[i-1]+1
  tkrnk = LETTERS[tkrnk]
  names(tkrnk) = rnk
  
  # bar plot
  if(Tukey){out = tkrnk}else{out = ttrnk}
  barplot(mu, ylim = c(0,max(upr)),space=F,main='Bar plot with SD and class',col=COL)
  for(i in 1:nl) lines(c(i-.5,i-.5),c(lwr[i],upr[i]),lwd=2)
  for(i in 1:nl) text(i-.35,mu[i],out[lvl][i],cex=1.25)
  
  # Return output
  return(out)
}
compiler::cmpfun(tlm)