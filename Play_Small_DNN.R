#ActFun = function(x){x[x<0]=0;return(x)} # RELU
ActFun = function(x){x[x<0]=x[x<0]*0.25;return(x)} # Leaky ReLU
#ActFun = function(x) x # Linear

# Data
require(bWGR)
data(tpod)
# Centralize
y = matrix(y-mean(y))
X = CNT(gen)
# Dimensions
n = nrow(X)
p = ncol(X)
k = ncol(y)
# Number of nodes
n1 = 20
n2 = 20
# Starting weights
b1 = matrix(rnorm(n1*p,0,1/p),p,n1)
b2 = matrix(rnorm(n1*n2,0,1/n1),n1,n2)
b3 = matrix(rnorm(n2,0,1/n2),n2,k)
# Iterations (backprop)
CNV1 = CNV2 = c()
for(i in 1:100){  
  cat('Iteration',i,'\n')
  # Fit hidden layers (H)
  H1 = ActFun(X%*%b1)
  H2 = ActFun(H1%*%b2)
  H3 = H2%*%b3
  # Gradients
  e3 = y-H3
  if(anyNA(e3)) e3[is.na(e3)]=0
  e2 = ActFun(e3 %*% t(b3))
  e1 = ActFun(e2 %*% t(b2))
  # Update coefficients
  b1 = b1 - ( 2*t(X)%*%(H1-e1)*0.0001  )
  b2 = b2 - ( 2*t(H1)%*%(H2-e2)*0.01  )
  b3 = b3 - ( 2*t(H2)%*%(H3-e3)*0.1  )
  # Print Residual Variance
  CNV1 = c(CNV1,mean(apply(e3,2,var,na.rm=T)))
  CNV2 = c(CNV2,mean(diag(cor(H3,y,use='p'))))
}


par(mfrow=c(1,2))
plot(CNV1,ylab='MSE',type='l')
plot(CNV2,ylab='COR',type='l')
