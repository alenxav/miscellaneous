# Multi-Trait Model Via Expectation-Maximization
MEM = function(Y,Z,nit=200,cnv=1e-10){
  
  # Template
  if(F){
    require(bWGR)
    data(tpod)
    Y = cbind(y1=y,y2=y)
    Y[fam==1,1]=NA
    Y[fam==2,2]=NA
    Z = gen
    rm(y,chr,fam,gen)
    nit=200
  }
  
  # Obs
  tr = function(x) sum(diag(x))
  require(Matrix)
  W = !is.na(Y)
  n = colSums(W)
  q = ncol(Z)
  k = ncol(Y)
  
  # Variances
  I = Matrix::diag(q)
  Vy = var(Y,use='p')
  Vy[is.na(Vy)] = 0
  MSx = mean(apply(Z,2,crossprod))
  Sa = Se = (Vy*0.5)
  Sa = Sa/MSx
  Se = diag(diag(Se))
  Traces = matrix(0,k,k)
  
  # Baseline for LHS
  Z = apply(Z,2,function(x) x-mean(x) )
  ZpZ0 = lapply(data.frame(W),function(x) crossprod(Matrix(Z[x,])))
  
  ZpZ = ZpZ0
  for(z in 1:k) ZpZ[[z]] = ZpZ0[[z]]*Se[z,z]
  ZpZ = bdiag(ZpZ)
  
  # Baseline for RHS
  Y0 = Y
  Mu = colMeans(Y0,na.rm = T)
  Y = apply(Y,2,function(x) x-mean(x,na.rm=T) )
  y = c(Y)
  y = y[!is.na(y)]
  Y[is.na(Y)] = 0
  Zy = c(crossprod(Z,Y))
  
  # Loop
  CONV = rep(nit)
  for(i in 1:nit){
    
    Sa0 = Sa
    Se0 = Se
    cat(i,'\n')
    
    ZpZ = ZpZ0
    YY = Y
    for(z in 1:k){
      ZpZ[[z]] = ZpZ0[[z]]/Se[z,z]
      YY[,z] = Y[,z]/Se[z,z]
    } 
    ZpZ = bdiag(ZpZ)
    Zy = c(crossprod(Z,YY))
    
    lmb = solve(Sa)
    Lmb = kronecker(lmb,I)
    
    iC = solve(ZpZ+Lmb)
    u = crossprod(iC,Zy)
    U = matrix(u,ncol=k)
    
    for(a in 1:k){
      for(b in 1:k){
        a0 = (q*(a-1)+1):(a*q)
        b0 = (q*(b-1)+1):(b*q)
        Traces[a,b] = tr(iC[a0,b0])
      }
    }
    
    G = Z%*%U
    E = (Y-G)*W
    Se = diag(diag(crossprod(E))/(n-1))
    Sa = (crossprod(U)+Traces)/q
    
    print(round(cov2cor(Sa),4))
    cat('\n')
    
    #conv = crossprod(c(c(Sa-Sa0),c(Se-Se0)))
    conv = sum(abs(c(c(Sa-Sa0),c(Se-Se0))))
    CONV[i] = conv
    
    if( conv<(cnv) ) break()
    
  }
  
  # Output
  CONV = log10(CONV[!is.na(CONV)])
  HAT = G;
  H2 = 1-diag(Se)/diag(Vy)
  for(i in 1:k) HAT[,i]=G[,i]+Mu[i]
  return(list(mu=Mu,b=U,hat=HAT,GC=cov2cor(Sa),h2=H2,Va=Sa,Ve=Se,conv=CONV))
  
}

# Test
if(F){
  require(bWGR)
  data(tpod)
  Y = cbind(y1=y,y2=y)
  Y[fam==1,1]=NA
  Y[fam==2,2]=NA
  Z = gen
  rm(y,chr,fam,gen)
  nit=200
  fit = MEM(Y,Z)
  h = mrr(Y,Z)
  c(h2=fit$h2,gc=fit$GC[1,2]);plot(fit$conv)
  c(h2=h$h2,gc=h$GC[1,2])
}
