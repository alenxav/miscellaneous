# Extreme learning machine
ELM = function(Y,Z,HN1=800,HN2=400,R2=0.75){
  p = ncol(Z)
  w1 = matrix(rnorm(p*HN1,0,1/sqrt(p)),p,HN1)
  h1 = Z%*%w1
  h1[h1<0] = 0
  w2 = matrix(rnorm(HN1*HN2,0,1/sqrt(HN1)),HN1,HN2)
  h2 = h1%*%w2
  h2[h2<0] = 0
  h2 = cbind(1,h2)
  LHS = crossprod(h2)
  RHS = crossprod(h2,Y)
  lmb = mean(diag(LHS)[-1])*((1-R2)/R2)
  diag(LHS)[-1] = diag(LHS)[-1] + lmb
  w3 = solve(LHS,RHS)
  hat = h2 %*% w3
  out = list(hat=hat,w1=w1,w2=w2,w3=w3)
  class(out) = 'ELM'
  return(out)
}

# Function to run predictions
predict.ELM = function(object,newdata){
  h1 = newdata%*%object$w1
  h1[h1<0] = 0
  h2 = h1%*%object$w2
  h2[h2<0] = 0
  h2 = cbind(1,h2)
  hat = h2 %*% object$w3
  return(hat)
}

# Run example?
if(FALSE){
  
  # Get data and subset training set
  data(wheat,package = 'BGLR')
  set.seed(12345)
  set = sample(1:nrow(wheat.Y),nrow(wheat.Y)*0.8)
  
  # Fit ELM model with 80% of the data
  FIT = ELM(wheat.Y[set,], wheat.X[set,])
  PRED = predict(FIT, wheat.X[-set,])

  # Baseline: BayesB
  FIT2 = apply(wheat.Y[set,], 2, bWGR::emBB, gen=wheat.X[set,])
  FIT2 = sapply(FIT2, function(x) x$b )
  PRED2 = wheat.X[-set,] %*% FIT2
  
  # Compare
  PA = cor(PRED, wheat.Y[-set,])
  PA2 = cor(PRED2, wheat.Y[-set,])
  diag(PA)
  diag(PA2)
  
}
