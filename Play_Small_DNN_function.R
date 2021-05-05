
dnn = function(y, X,
                nit=1000, batch=250,
                RELU=FALSE, Leak=0.1,
                dropout=0, Lambda=0.1,
                LrnRate = 1,
                Nodes1=4, Nodes2=4){
  # Normalization
  if(is.null(ncol(y))) y = matrix(y)
  muY = colMeans(y,na.rm=T)
  sdY = apply(y,2,sd,na.rm=T)
  y = apply(y,2,scale)
  # DNN functions
  ActFun = tanh
  if(RELU) ActFun = function(x){x[x<0]=x[x<0]*Leak;return(x)}
  DropOut = function(x,prc=dropout){x[sample(length(x),length(x)*prc)];return(x)}
  # Dimensions
  n = nrow(X)
  p = ncol(X)
  k = ncol(y)
  # Learning settings
  n1 = Nodes1; n2 = Nodes2;
  lmb = Lambda; rate = LrnRate/c(p,n1,n2)
  # Starting weights
  b1 = matrix(rnorm(n1*p,0,1/p),p,n1)

dnn = function(y, X,
                nit=1000, batch=250,
                RELU=FALSE, Leak=0.1,
                dropout=0, Lambda=0.1,
                LrnRate = 1,
                Nodes1=4, Nodes2=4){
  # Normalization
  if(is.null(ncol(y))) y = matrix(y)
  muY = colMeans(y,na.rm=T)
  sdY = apply(y,2,sd,na.rm=T)
  y = apply(y,2,scale)
  # DNN functions
  ActFun = tanh
  if(RELU) ActFun = function(x){x[x<0]=x[x<0]*Leak;return(x)}
  DropOut = function(x,prc=dropout){x[sample(length(x),length(x)*prc)];return(x)}
  # Dimensions
  n = nrow(X)
  p = ncol(X)
  k = ncol(y)
  # Learning settings
  n1 = Nodes1; n2 = Nodes2;
  lmb = Lambda; rate = LrnRate/c(p,n1,n2)
  # Starting weights
  E1 = b1 = matrix(rnorm(n1*p,0,1/p),p,n1)
  E2 = b2 = matrix(rnorm(n1*n2,0,1/n1),n1,n2)
  E3 = b3 = matrix(rnorm(n2,0,1/n2),n2,k)
  # Iterations (backprop)
  CNV1 = CNV2 = c()
  for(i in 1:nit){ 
    if((i+1)%%100==1) cat('Iteration',i,'\n')
    # Sample batch
    w = sample(n,batch,replace=T)
    y0 = y[w,]
    X0 = X[w,]
    # Fit hidden layers (H)
    H1 = ActFun(X0%*%b1)
    H2 = ActFun(H1%*%b2)
    H3 = H2%*%b3
    # Gradients
    e3 = y0-H3; if(anyNA(e3)) e3[is.na(e3)]=0
    e2 = ActFun(e3 %*% t(b3))
    e1 = ActFun(e2 %*% t(b2))
    # Update coefficients
    GRD1 = DropOut(t(X0)%*%(e1)-lmb*b1)*(2/n)*rate[1]
    GRD2 = DropOut(t(H1)%*%(e2)-lmb*b2)*(2/n)*rate[2]
    GRD3 = DropOut(t(H2)%*%(e3))*(2/n)*rate[3]
    b1 = b1 + GRD1 + E1*0.5
    b2 = b2 + GRD2 + E2*0.5
    b3 = b3 + GRD3 + E3*0.5
    # Store convergence
    CNV1 = c(CNV1,mean(apply(e3,2,var,na.rm=T)))
    CNV2 = c(CNV2,mean(diag(cor(H3,y0,use='p'))))
    # Store current gradient
    E1 = GRD1
    E2 = GRD2
    E3 = GRD3
  }
  out = list(af=ActFun,b1=b1,b2=b2,b3=b3,conv_MSE=CNV1,conv_GOF=CNV2,mu=muY,sd=sdY)
  class(out) = 'smalldnn'
  return(out)
}

predict.smalldnn = function(object,newdata){
  x = object
  h = x$af(x$af(newdata%*%x$b1)%*%x$b2)%*%x$b3
  for(i in 1:ncol(h)) h[,i] = x$mu[i] + x$sd[i]*h[,i]
  return(h)
}



#############################
# RUN EXAMPLE IN WHEAT DATA #
############################

RUN_EXAMPLE = T

if(RUN_EXAMPLE){
  
  # Run DNN
  data(wheat,package = 'BGLR')
  fit = dnn(wheat.Y,wheat.X)
  hat = predict(fit,wheat.X)
  cat('DNN fitness\n')
  print(diag(cor(hat,wheat.Y)))
  
  # Run BayesA
  ba = apply(wheat.Y,2,function(y) bWGR::emBA(y,wheat.X)$hat )
  cat('BayesA fitness\n')
  print(diag(cor(ba,wheat.Y)))
  
  # Compare DNN with BayesA
  cat('BayesA vs DNN\n')
  print(diag(cor(ba,hat)))
  
}
