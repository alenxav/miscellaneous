tlm = function(formula,data,Tukey=T,...){
  
  # Check the variable under evaluation
  key = gsub(' .+','',as.character(formula)[3])
  data[[key]] = factor(data[[key]])
  nl = length(levels(data[[key]]))
  
  # fit the model
  fnew = update(formula,~.-1)
  fit = lm(formula=fnew,data=data)
  stats = summary(fit)$coefficients
  stats = stats[grepl(key,rownames(stats)),]

  Y = predict(fit,data,terms=key,type='term')+fit$residuals
  data_new = fit$model
  data_new[['Y']] = Y
  X = data_new[key][,1]
  lvl = levels(X)
  
  # Refit and multiple testing
  refit = lm(Y~X-1,data_new)
  cat('overall p-value',c(anova(lm(Y~X-1,data_new))["Pr(>F)"])[[1]][1],'\n')
  tk = TukeyHSD(aov(refit))
  
  # Pool SD if it doesn't work otherwise
  tt = try(pairwise.t.test(Y,X,pool.sd=F),silent = T)
  if(class(tt)=="try-error")  tt = pairwise.t.test(Y,X,pool.sd=T)
  
  # mean and sd # Debugged with Rajat, 10/12/2023
  # mu = tapply(Y,X,mean)
  # std =  tapply(Y,X,sd)
  rownames(stats) = gsub(key,'',rownames(stats))
  mu = stats[,1]
  std = stats[,2]

  
  std[is.na(std)] = mean(std,na.rm = T)
  
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
  barplot(mu,..., ylim = c(min(lwr)*0.8,max(upr)*1.2),space=F,col=COL,xpd = F)
  for(i in 1:nl) lines(c(i-.5,i-.5),c(lwr[i],upr[i]),lwd=2)
  for(i in 1:nl) text(i-.2,(mu[i]*0.3+upr[i]*0.7),out[lvl][i],cex=1.25)
  
  # Return output
  return(out)
}
