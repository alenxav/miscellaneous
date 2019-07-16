
# ridge regression function
rr = function(y,gen){
  
  n = nrow(gen)
  p = ncol(gen)
  
  xx = apply(gen,2,crossprod)
  msx = lmb = sum(apply(gen,2,var))
  b = rep(0,p)
  mu = mean(y)
  e = y-mu
  
  for(i in 1:100){
    
    # Gauss-Seidel
    sapply(1:p,function(i){
      b0 = b[i]; 
      b[i] <<- (crossprod(e,gen[,i])[1,1]+xx[i]*b[i])/(xx[i]+lmb)
      e <<- e-gen[,i]*(b[i]-b0); return() })
    
    # Intercept
    mu0 = mean(e)
    mu = mu+mu0
    e = e-mu0
    
    # Variance components
    vy = crossprod(y-mu,y-mu)[1,1]/(n-1)
    ve = crossprod(y,e)[1,1]/(n-1)
    vb = (vy-ve)/msx
    lmb = ve/vb
    
  }
  
  fit = y-e
  out = list(mu=mu,b=b,hat=fit,vb=vb,ve=ve,h2=1-ve/vy)
  return(out)
  
}

compiler::cmpfun(rr)

# data(tpod,package='bWGR'); fit=rr(y,gen); plot(fit$b)