require(NAM)
if(!exists('Y')){
  tmp = SoyNAM::BLUP(family=1:5)
  Y = tmp$Phen
  Z = CNT(tmp$Gen)
  rm(tmp)
  xx = apply(Z,2,crossprod)
  msx = mean(xx)
}

# Gauss-Seidel Gradients
gsg = function(Y,X,h2=1){
  xx = apply(X,2,crossprod)+1e-8
  lmb = mean(xx)*(1-h2)/h2
  #gs = function(Y) NAM::NOR(Y,X,lmb,xx,1)$b
  gs = function(Y){
    h = NAM::NOR(Y,X,lmb,xx,1)
    cat(h$v,'\n')
    return(h$b) } 
  w = apply(Y,2,gs)
  return(w)}

n = nrow(Z)
p = ncol(Z)

# Activation function
#ActFun = function(x){ 1/(1+exp(-x)) }  # Sigmoid
#ActFun = tanh  # Tanh
#ActFun = function(x){ x[x<0]=0; return(x) }  # ReLU
#ActFun = function(x){ x[x<0]=x[x<0]/100; return(x) }  # Leaky ReLU
ActFun = function(x){ x }  # linear

# Fit node
FN = function(Z,W,I){
  tmp = Z %*% W
  tmp = t(apply(tmp,1,function(x)x+I))
  tmp = ActFun(tmp)
  return(tmp)
}

# Number of nodes in hidden layer
HL1 = 50
HL2 = 20

# Some normalization function
nrm = function(X) apply(X,2,function(x) x/sqrt(c(crossprod(x))))

# Intercepts (I), Weights (W) and Fits (F)
I1 = rnorm(HL1)
W1 = nrm(matrix(runif(p*HL1),p,HL1))
H1 = FN(Z,W1,I1)
I2 = rnorm(HL2)
W2 = nrm(matrix(runif(HL1*HL2),HL1,HL2))
H2 = FN(H1,W2,I2)
I3 = mean(Y)
W3 = rnorm(HL2)
H3 = c(FN(H2,W3,I3))

# Check starting point
CHECK = function(){
  Hat = c(FN(FN(FN(Z,W1,I1),W2,I2),W3,I3))
  par(mfrow=c(1,1)); plot(Hat,Y,main=paste("COR =",round(cor(Hat,Y),4)))}
CHECK()

# Learning rate and L2 penalization
rates = c(0.25,0.5,1.0)
lmb = c(0.01,0.01,0.01)


########################### SECOND ATTEMPT

if(T){
  
  for(iter in 1:8){
    cat('\n ITERATION',iter,'\n')
    
    # Batch
    BatchSize = 0.5
    mb = sort(sample(1:n,n*BatchSize)) # mini batch
    
    # Fit model
    H1 = FN(Z[mb,],W1,I1)
    H2 = FN(H1,W2,I2)
    Hat = c(FN(H2,W3,I3))
    
    # Layer 3
    # Backprop
    dH3 = matrix(Y[mb]-Hat)
    #dW3 = t(H2)%*%dH3
    dW3 = gsg(dH3,H2)
    dI3 = mean(dH3)
    # Update parameter 
    W3 = W3 + c(dW3)*rates[2] #- (dW3+lmb[3])*rates[3]
    I3 = I3 - dI3*rates[3]
    CHECK()
    
    # Layer 2
    # Backprop
    dH2 = ActFun( dH3 %*% t(W3) )
    #dW2 = t(H1)%*%dH2
    dW2 = gsg(dH2,H1,h2=0.5)
    dI2 = colMeans(dH2)
    # Update parameter 
    W2 = W2 + dW2*rates[2] #- (dW2+lmb[2])*rates[2]
    I2 = I2 - dI2*rates[2]
    CHECK()
    
    # Layer 1
    # Backprop
    dH1 = ActFun( dH2 %*% t(W2) ) # y
    #dW1 = t(Z[mb,])%*%dH1 # rhs
    dW1 = gsg(dH1,Z[mb,],h2=0.5)
    dI1 = colMeans(dH1)
    # Update parameter 
    W1 = W1 + dW1*rates[1] #- (dW1+lmb[1])*rates[1]
    I1 = I1 - dI1*rates[1]
    CHECK()
    
    
  }
  
}


