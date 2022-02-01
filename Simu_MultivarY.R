SimY = function(Z, k=10, h2=0.5,GC=0.5,  seed=123, unbalanced=FALSE){
  
  # Store inputs
  trueVal = c(h2=h2,GC=GC,seed=seed)
  n = nrow(Z)
  p = ncol(Z)
  
  # Pick a dataset
  set.seed(seed)
  
  # Genetic parameters
  h20 = h2
  h2 = rep(h20,k)
  
  # GC
  numGC = (k*(k-1))/2
  GC = rep(GC,numGC)
  G0 = diag(k)
  G0[lower.tri(G0)] = GC
  G0[upper.tri(G0)] = t(G0)[upper.tri(t(G0))]
  GC = mean(GC)
  
  # Sample effects
  alpha = 1/sum(apply(Z,2,var))
  trueVal['scaleG'] = alpha
  Vb = G0*alpha
  ev = eigen(Vb, symmetric = TRUE)
  UD = ev$vectors %*% diag(sqrt(ev$values))
  beta = matrix(rnorm(p * k), nrow = p)
  trueBeta = UD %*% t(beta)
  
  # True breeding values
  tbv = Z %*% t(trueBeta)
  
  # Residual variances and phenotypes
  ve = (1-h2)/h2;
  E = sapply(ve,function(ve) rnorm(n,0,sqrt(ve)))
  Y = 10 + tbv + E
  colnames(Y) = paste('y',1:k,sep='')
  
  # Unbalance
  if(unbalanced){
    Miss = sample(1:k,n,T)
    for(i in 1:k) Y[which(Miss!=i),i] = NA
  }
  
  return(list(Y=Y,tbv=tbv,settings=trueVal))
  
}