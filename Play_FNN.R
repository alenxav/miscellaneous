require(NAM)
if(!exists('Y')){
  tmp = SoyNAM::BLUP(family=1:5)
  Y = tmp$Phen
  Z = CNT(tmp$Gen)
  rm(tmp)
}

xx = apply(Z,2,crossprod)
msx = mean(xx)

n = nrow(Z)
p = ncol(Z)

# Activation function
#ActFun = function(x){ 1/(1+exp(-x)) }  # Sigmoid
#ActFun = tanh  # Tanh
ActFun = function(x){ x[x<0]=0; return(x) }  # ReLU
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

# Intercepts (I), Weights (W) and Fits (F)
I1 = rnorm(HL1,sd=0.1)
W1 = matrix(rnorm(p*HL1,sd=0.1),p,HL1)
H1 = FN(Z,W1,I1)
I2 = rnorm(HL2,sd=0.1)
W2 = matrix(rnorm(HL1*HL2,sd=0.1),HL1,HL2)
H2 = FN(H1,W2,I2)
I3 = mean(Y)
W3 = GD(Y-I3,H2)$b
H3 = c(FN(H2,W3,I3))

# Check starting point
hh = H3
plot(hh,Y,main='Fitted')
ee = c(crossprod(Y-hh))

# Learning rate and L2 penalization
rates = c(0.25,0.5,1.0)
lmb = c(msx,0.01,0.01)

########################### FIRST ATTEMPT

if(F){
  
  # Fill layers
  Y0 = Y
  for(iter in 1:3){
    
    cat('\n ITERATION',iter,'\n')
    for(i in 1:HL1){
      
      cat('.')
      
      # Fill H1
      Y0 = Y0 + c(ActFun(ActFun(Z%*%W1[,i]+I1[i])%*%W2[i,]+I2[i])%*%W3) 
      tmp = GD(Y0,Z,W1[,i],xx=xx,rate=rates[1],lmb=lmb[1])
      W1[,i] = tmp$b
      H1[,i] = tmp$g
      I1[i] = tmp$mu
      # plot(tmp$b)
      # plot(tmp$g,Y0)
      
      # Fill H2
      for(j in 1:HL2){
        
        Y00 = Y0 + c(ActFun(ActFun(Z%*%W1[,j]+I1[j])%*%W2[j,]+I2[j])%*%W3)-I3
        tmp = GD(Y00,ActFun(H1),W2[,j],rate=rates[2],lmb=lmb[2])
        W2[,j] = tmp$b
        H2[,j] = ActFun(tmp$g)
        I2[j] = tmp$mu
        # plot(tmp$b)
        # plot(tmp$g,Y0)
        
        # Fill H3
        tmp = GD(Y,ActFun(H2),rate=rates[3],lmb=lmb[3])
        W3 = tmp$b
        F3 = ActFun(tmp$g)
        I3 = tmp$mu
        # plot(tmp$b)
        # plot(tmp$g,Y0)
        
        # Update target
        Y0 = Y - F3
        
      }
      
      
    }
    cat('\n R',round(cor(F3,Y),2),'\n')
    
    cat('\n')
  }
  
  # Plot diagnostics
  par(mfrow=c(2,1))
  
  # Check fitness
  Hat = F3+I3
  plot(Hat,Y,main='Propagated')
  
  # Double check fitting NN
  hh = FN(FN(Z,W1,I1),W2,I2)%*%W3+I3
  plot(hh,Y,main='Fitted')
  
  
}

########################### SECOND ATTEMPT

par(mfrow=c(2,2))

if(T){
  
  for(iter in 1:8){
    cat('\n ITERATION',iter,'\n')
    
    # Fit model
    
    H1 = FN(Z,W1,I1)
    H2 = FN(H1,W2,I2)
    Hat = c(FN(H2,W3,I3))
    
    # Backprop
    
    # Layer 3
    dH3 = Y-Hat
    dW3 = t(H2)%*%dH3
    dI3 = mean(dH3)
    # Layer 2
    dH2 = ActFun( dH3 %*% t(W3) )
    dW2 = t(H1)%*%dH2
    dI2 = mean(dH2)
    # Layer 1
    dH1 = ActFun( dH2 %*% t(W2) )
    dW1 = t(Z)%*%dH1
    dI1 = mean(dH1)
    
    # Add L2 penalty
    dW3 = dW3+lmb[3]*W3
    dW2 = dW2+lmb[2]*W2
    dW1 = dW1+lmb[1]*W1
    
    # Update parameter 
    W3 = W3 - (dW3+lmb[3]*W3)*rates[3]
    W2 = W2 - (dW2+lmb[3]*W2)*rates[2]
    W1 = W1 - (dW1+lmb[3]*W1)*rates[1]
    I3 = I3 - dI3*rates[3]
    I2 = I2 - dI2*rates[2]
    I1 = I1 - dI1*rates[1]
    
    # Check fitness
    plot(Hat,Y,main=paste('Fitted',iter))
    noise = rnorm(n)
    cnoise = round(cor(rnorm(n),Y),2)
    cat('\n R =',round(cor(Hat+noise,Y),2)-cnoise,'\n')
    
  }
  
}


