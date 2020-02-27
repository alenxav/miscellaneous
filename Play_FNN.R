require(NAM)
if(!exists('Y')){
  tmp = SoyNAM::BLUP(family=1:5)
  Y = tmp$Phen
  Z = CNT(tmp$Gen)
}
rm(tmp)
xx = apply(Z,2,crossprod)
n = nrow(Z)
p = ncol(Z)

# Gradient descent
GD = function(y,Z,b=NULL,xx=NULL,h2=0.25,rate=0.1){
  if(is.null(xx)) xx = apply(Z,2,crossprod)
  if(is.null(b)) b = rep(0,ncol(Z))
  mu = mean(y)
  lmb = mean(xx)*(1-h2)/h2
  e = y-mu
  for(j in 1:ncol(Z)){
    b0 = b[j]
    b1 = (c(e%*%Z[,j])*rate+b0*xx[j])/(xx[j]+lmb+1e-8)
    b1 = c(b1)
    b[j] = b1
    e = e - c(Z[,j])*(b1-b0)}
  g = y-e-mu
  out = list(b=b,g=g,mu=mu)
}
compiler::cmpfun(GD)

# Activation function
#ActFun = function(x){ 1/(1+exp(-x)) }  # Sigmoid
ActFun = function(x){ x[x<0]=0; return(x) }  # ReLU
#ActFun = function(x){ (exp(x)-exp(-x))/(exp(x)+exp(-x)) }  # Tanh


# Fit node
FN = function(Z,W,I){
  tmp = Z %*% W
  tmp = t(apply(tmp,1,function(x)x+I))
  tmp = ActFun(tmp)
  return(tmp)
}

# Number of nodes in hidden layer
HL1 = 10
HL2 = 10

# Intercepts, Weights and Fits
I1 = rep(0,HL1)
I2 = rep(0,HL2)
I = 0
W1 = matrix(0,p,HL1)
W2 = matrix(0,HL1,HL2)
W3 = rep(0,HL2)
F1 = matrix(0,n,HL1)
F2 = matrix(0,n,HL2)
F3 = rep(0,n)

# Learning rate and L2 penalization
rates = c(0.1,0.8,1)
H2 = c(0.5,0.9,1)

# Fill layers
Y0 = Y
for(iter in 1:3){
  
  cat('\n ITERATION',iter,'\n')
  for(i in 1:HL1){
    
    cat('.')
    
    # Fill H1
    tmp = GD(Y0,Z,W1[,i],xx=xx,rate=rates[1],h2=H2[1])
    W1[,i] = tmp$b
    F1[,i] = tmp$g
    I1[i] = tmp$mu
    # plot(tmp$b)
    # plot(tmp$g,Y0)
    
    # Fill H2
    for(j in 1:HL2){
      
      tmp = GD(Y0,ActFun(F1),W2[,j],rate=rates[2],h2=H2[2])
      W2[,j] = tmp$b
      F2[,j] = tmp$g
      I2[j] = tmp$mu
      # plot(tmp$b)
      # plot(tmp$g,Y0)
      
      # Fill H3
      tmp = GD(Y,ActFun(F2),rate=rates[3],h2=H2[3])
      W3 = tmp$b
      F3 = tmp$g
      I3 = tmp$mu
      # plot(tmp$b)
      # plot(tmp$g,Y0)
      
      # Update target
      Y0 = Y - F3
      
    }
    
    
  }
  cat('R2',round(cor(F3,Y),2),'\n')
  
  cat('\n')
}

# Plot diagnostics
par(mfrow=c(2,1))

# Check fitness
Hat = F3+mu
plot(Hat,Y,main='Propagated')

# Double check fitting NN
hh = FN(FN(Z,W1,I1),W2,I2)%*%W3+I3
plot(hh,Y,main='Fitted')



