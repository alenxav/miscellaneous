if(!exists('Y')){
  require(NAM)
  tmp = SoyNAM::BLUP(family=1:12)
  fam = tmp$Fam
  # Train
  Y = tmp$Phen[fam<=8]
  Z = CNT(tmp$Gen[fam<=8,])
  # Test
  Y2 = tmp$Phen[fam>8]
  Z2 = CNT(tmp$Gen[fam>8,])
  rm(tmp)
  xx = apply(Z,2,crossprod)
  msx = mean(xx)
}

# Gauss-Seidel Gradients
gsg = function(Y,X,xx=NULL,h2=1){
  if(is.null(xx)) xx = apply(X,2,crossprod)+1e-8
  lmb = mean(xx)*(1-h2)/h2
  gs = function(Y) NAM::NOR(Y,X,lmb,xx,5)$b
  coefs = apply(Y,2,gs)
  return(coefs)}

n = nrow(Z)
p = ncol(Z)

# Activation function
#ActFun = function(x){ x }  # linear
#ActFun = function(x){ 1/(1+exp(-x)) }  # Sigmoid
#ActFun = tanh  # Tanh
ActFun = function(x){ x[x<0]=0; return(x) }  # ReLU
#ActFun = function(x){ x[x<0]=x[x<0]/100; return(x) }  # Leaky ReLU

# Fit node
FN = function(Z,W,I,AF=TRUE){
  tmp = Z %*% W
  for(k in 1:length(I)) tmp[,k] = tmp[,k]+I[k]
  if(AF) tmp = ActFun(tmp)
  return(tmp)
}

# Number of nodes in hidden layer
HL1 = 40
HL2 = 20

# Some normalization function
nrm = function(X) apply(X,2,function(x) x/sqrt(c(crossprod(x))))

# Set seed
set.seed(12345)

# Intercepts (I), Weights (W) and Fits (F)
I1 = rnorm(HL1)
W1 = nrm(matrix(rnorm(p*HL1,sd=2/sqrt(p)),p,HL1))
H1 = FN(Z,W1,I1)
I2 = rnorm(HL2)
W2 = nrm(matrix(rnorm(HL1*HL2,sd=2/sqrt(HL1)),HL1,HL2))
H2 = FN(H1,W2,I2)
I3 = mean(Y)
W3 = rnorm(HL2,sd=2/sqrt(HL2))
H3 = c(FN(H2,W3,I3))

# Check starting point
CHECK = function(a='check'){
  Hat = c(FN(FN(FN(Z,W1,I1),W2,I2),W3,I3,AF=FALSE))
  plot(Hat,Y,main=paste("COR =",round(cor(Hat,Y),4),a))}
CHECK()

# Store fitness and prediction
GOF = PA = 0

###########################

# Learning rate, L2 penalization and epochs
rate = 0.02
h2 = 0.5
epochs = 15

###########################

# Dinamics of learning rate
update_rate = TRUE
rate_update_factor = 1.15
max_rate = 0.5
perturb_y = 0.1

###########################

# RUN DNN
for(iter in 1:epochs){
  cat('ITERATION',iter,'\n')
  
  # Add some noise
  Ytmp = Y + rnorm(n,sd=sd(Y)*perturb_y)
  
  # Layer 1 backprop
  Hat = c(FN(FN(FN(Z,W1,I1),W2,I2),W3,I3,AF=FALSE))
  dH3 = matrix(Ytmp-Hat)
  dH2 = ActFun(dH3%*%t(W3))
  dH1 = ActFun(dH2%*%t(W2))
  dW1 = gsg(dH1,Z,xx,h2=h2)
  dI1 = colMeans(dH1)
  W1 = W1 + dW1*rate
  I1 = I1 + dI1*rate
  CHECK(' Layer 1')
  
  # Layer 2 backprop
  Hat = c(FN(FN(FN(Z,W1,I1),W2,I2),W3,I3,AF=FALSE))
  dH3 = matrix(Ytmp-Hat)
  dH2 = ActFun( dH3 %*% t(W3) )
  dW2 = gsg(dH2,H1,h2=h2)
  dI2 = colMeans(dH2)
  W2 = W2 + dW2*rate
  I2 = I2 + dI2*rate
  CHECK(' Layer 2')
  
  # Layer 3 backprop
  Hat = c(FN(FN(FN(Z,W1,I1),W2,I2),W3,I3,AF=FALSE))
  dH3 = matrix(Ytmp-Hat)
  dW3 = gsg(dH3,H2,h2=h2)
  dI3 = mean(dH3)
  W3 = W3 + c(dW3)*rate
  I3 = I3 + dI3*rate
  CHECK(' Layer 3')
  
  # Update rate
  if(update_rate) rate = min(max_rate,rate*rate_update_factor)
  
  # Store PA
  GOF[iter] = cor(Y,c(FN(FN(FN(Z,W1,I1),W2,I2),W3,I3,AF=FALSE)))
  PA[iter] = cor(Y2,c(FN(FN(FN(Z2,W1,I1),W2,I2),W3,I3,AF=FALSE)))
  cat('Rate',round(rate,4),'GOF/PA',round( c(GOF[iter],PA[iter]) ,4),'\n')
  
}

###########################

# Compute a baseline
rr = emML(Y,Z); rrhat=Z2%*%rr$b; parr = cor(Y2,rrhat); parr
bb = emBB(Y,Z); bbhat=Z2%*%bb$b; pabb = cor(Y2,bbhat); pabb
nnhat = c(FN(FN(FN(Z2,W1,I1),W2,I2),W3,I3,AF=FALSE))

# Checks prediction by epoch
par(mar=c(4,4,2,2))
plot(GOF,xlab='Iteration',type='l',main='Check performance of deep neural net',
     ylab='Correlation',lwd=3,ylim=c(0,1))
# Add baseline
lines(c(iter-2,iter),c(parr,parr),col=3,lwd=3)
lines(c(iter-2,iter),c(pabb,pabb),col=4,lwd=3)
lines(PA,col=2,lwd=3)
# Legend
legend('topleft',ncol=2,
       c('DNN fitness','DNN prediction','GBLUP prediction','BayesB prediction'),
       col=1:4,pch=20,bty='n')

# Within-family predictive ability
Q = data.frame(rrhat,bbhat,nnhat,y=Y2,fam=fam[fam>8])
Q =  by(data = Q,INDICES = Q[,5],FUN = function(q) cor(q[,1:4])[1:3,4] )
Q = sapply(Q, c)
Q