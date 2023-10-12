
# Creating latent spaces
Autoencoder = function(Y,NumLatSpa=10,maxit=500,minit=20,LogConv=5,rate=0.1,decay=0.95,lmb=0.001,
                 PrintEvery=100,DropOut=0.5,ObsDropOut=FALSE,ADAM=FALSE,verb=TRUE){
  Y = X = apply(Y,2,scale)
  if(ObsDropOut){
    RunDropOut = function(x,prc=.5){x[sample(length(x),length(x)*prc)];return(x)};}else{
    RunDropOut = function(x,prc=.5){x[,sample(ncol(x),round(ncol(x))*prc)];return(x)};}
  n = nrow(Y); X[is.na(Y)] = 0; k = ncol(Y); learning = 2/n*rate; if(ADAM){GG1 = GG0 = 0}; 
  V = matrix(rnorm(NumLatSpa*k)/k,k,NumLatSpa)*0.01;
  A = matrix(rnorm(NumLatSpa*k)/NumLatSpa,NumLatSpa,k)*0.01;
  hat = X %*% V %*% A;  E2 = (Y - hat);
  MSE0 = mean(apply(E2,2,var,na.rm=T))
  for(i in 1:maxit){
    H = X %*% V
    if(anyNA(E2)) E2[is.na(E2)]=0
    E1 = E2 %*% t(A)
    if(DropOut>0){
      G1 = RunDropOut((t(H)%*%E2+A*lmb)*learning,DropOut)
      G0 = RunDropOut((t(X)%*%E1+V*lmb)*learning,DropOut)
    }else{
      G1 = ((t(H)%*%E2+A*lmb)*learning)
      G0 = ((t(X)%*%E1+V*lmb)*learning)
    }
    A = A + G1; V = V + G0
    rate = rate*decay;
    if(ADAM){ A = A+GG1; V = V+GG0; GG1 = G1*0.5; GG0 = G0*0.5}
    hat = X %*% V %*% A; E2 = (Y - hat);
    MSE = mean(apply(E2,2,var,na.rm=T))
    cnv = -log10(abs(MSE-MSE0))
    if(!i%%PrintEvery & verb) cat('maxit',i,'| R2 =',1-round(MSE,4),'| cnv =',round(cnv,2),'\n')
    if(i>minit & MSE>MSE0){cat('Peaked R2 at',i,'\n'); break()}; V0=V; A0=A; MSE0=MSE;
    if(i>minit & cnv>LogConv){cat('Converged at',i,'\n');break()} }
  Hsd = apply(H,2,sd); H = apply(H,2,scale);
  for(i in 1:NumLatSpa){ V0[,i]=V0[,i]/Hsd[i]; A0[i,]=A0[i,]*Hsd[i];} 
  out = list(Encoder=V0,Decoder=A0,LatentSpace=X%*%V0,R2=1-MSE0,NumbOfIterations=i)}

# Abstraction function
MegaLmm = function(Y,Z,rerun=0,...){
  RR = function(y,X) bWGR::MRR3(matrix(y,ncol=1),X,verbose=F,NLfactor=1)
  UV = function(y,alg=RR){w=which(!is.na(y));return(c(alg(y[w],Z[w,])$b))}
  # Get FA
  cat('Fit FA ')
  FitFA = Autoencoder(Y,...)
  FA = FitFA$LatentSpace %*% FitFA$Decoder
  # Put FA in the phenotypic scale
  Scale_LatSp = sapply(1:ncol(Y), function(i){ w = which(!is.na(Y[,i])); return(qr.solve(cbind(1,FA[w,i]),Y[w,i]))}  )
  for(i in 1:ncol(Y)) FA[,i] = FA[,i]* Scale_LatSp[2,i]+Scale_LatSp[1,i]
  # Solve F
  cat('Fit F model\n')
  Fmodel = apply(FitFA$LatentSpace,2,UV); 
  # Final estimates
  cat('Fit R model\n')
  Beta_R = apply(Y-FA,2,UV)
  Hat = FA + Z%*%Beta_R
  Beta = t(t(Fmodel %*% FitFA$Decoder) * Scale_LatSp[2,]) + Beta_R
  A = t(t(FitFA$Decoder)*Scale_LatSp[2,])
  # Re-iterate
  if(rerun>0){
    for(i in 1:rerun){
      cat('Rerun -',i,'\n')
      FitFA = Autoencoder(Hat,...)
      FA = FitFA$LatentSpace %*% FitFA$Decoder
      # Put FA in the phenotypic scale
      Scale_LatSp = sapply(1:ncol(Y), function(i){ w = which(!is.na(Y[,i])); return(qr.solve(cbind(1,FA[w,i]),Y[w,i]))}  )
      for(i in 1:ncol(Y)) FA[,i] = FA[,i]* Scale_LatSp[2,i]+Scale_LatSp[1,i]
      # Solve F
      cat('Fit F model\n')
      Fmodel = apply(FitFA$LatentSpace,2,UV); 
      # Final estimates
      cat('Fit R model\n')
      Beta_R = apply(Y-FA,2,UV)
      Hat = FA + Z%*%Beta_R
      Beta = t(t(Fmodel %*% FitFA$Decoder) * Scale_LatSp[2,]) + Beta_R
      A = t(t(FitFA$Decoder)*Scale_LatSp[2,])
    }
  }
  # Output
  mu=colMeans(Y,na.rm=T)
  return(list(FA=FA,A=A,hat=Hat,mu=mu,b=Beta,gebv=Z%*%Beta))
}

