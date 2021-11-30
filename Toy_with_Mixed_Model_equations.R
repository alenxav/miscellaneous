require(bWGR)
data(tpod)
tr = function(x) sum(diag(x))
Y = sqrt(y); rm(y)

X = model.matrix(~factor(fam)-1)
Z = apply(gen,2,function(x) x-mean(x))/sqrt(ncol(gen))
h0 = NAM::reml(Y,X,Z)

b = h0$Fixed[,1]
u = h0$EBV
g = Z %*% u

n = nrow(Z)
q = ncol(Z)

vu = h0$VC$Vg
ve = h0$VC$Ve
I = diag(n)
G = diag(q)*vu
iG = solve(G)
R = I*ve
iR = solve(I*ve)
V = Z%*%G%*%t(Z) + R
iV = solve(V)  
Sigma = Matrix::bdiag(diag(0,2),iG)


W = cbind(X,Z)
C = t(W) %*% iR %*% W + Sigma
iC = solve(C)

C22 = iC[-c(1:2),-c(1:2)]
C22 = as.matrix(C22)

S = I - X %*% solve( t(X)%*%X ) %*% t(X)
H = W %*% iC %*% t(W) %*% iR

P = iV - iV %*% X %*% solve( t(X)%*%iV%*%X ) %*% t(X) %*% iV
plot(I-H,P*ve)
plot(iR%*%(I-H),P)

P2 = S%*%iV%*%S


YY = tcrossprod(S%*%Y)
diag(YY) = diag(YY)+0.0001
iYY = solve(YY)

VarG = G %*% t(Z) %*% iV %*% Z %*% G
VarG2 = G %*% t(Z) %*% iYY %*% Z %*% G

plot(diag(VarG2),diag(VarG))


u2 = G %*% t(Z) %*% P %*% Y
plot(u2,u)

yHat = X%*%b + Z%*%u
plot(P%*%Y*ve,Y-yHat)


# check the use of iV X
vx = iV %*% X
vx2 = (iR - iR %*% Z %*% solve( t(Z) %*% iR %*% Z + iG ) %*% t(Z) %*% iR ) %*% X
vx3 = iR %*% X - iR %*% Z %*% solve( t(Z) %*% iR %*% Z + iG ) %*% t(Z) %*% iR  %*% X

# check VC

# MIVQUE
Vi = tcrossprod(Z)
ss = t(Y) %*% P %*% Vi %*% P %*% Y
df = tr( P %*% Vi )# %*% P %*% Vi )
ss/df


# MIVQUE for both VCs
solve(a=matrix(c(
  tr( P %*% Vi %*% P %*% Vi ),
  tr( P %*% Vi %*% P %*% I ),
  tr( P %*% I %*% P %*% Vi ),
  tr( P %*% I %*% P %*% I )),2,2),
      b=c(t(Y) %*% P %*% Vi %*% P %*% Y,
             t(Y) %*% P %*%  I %*% P %*% Y))


# EM
ss = crossprod(u) + tr(C22)
df = q
ss/df
vu

e = Y - yHat
ss = crossprod(Y,e)
df = tr(S)
ss/df
ve


# Schaeffer's
ss = t(Y) %*% S %*% Z %*% u
df = tr( S %*% Vi )
ss/df
vu
ss = t(Y) %*% S %*% e
df = tr( S )
ss/df
ve

# AI
SecDer = t(Y) %*% P %*% Vi %*% P %*% Vi %*% P %*% Y
FirDer = -0.5 * ( q/vu - crossprod(u)/(vu^2) - tr(C22)/(vu^2) )
vu + FirDer / SecDer
vu
SecDer = t(Y) %*% P %*% I %*% P %*% I %*% P %*% Y
FirDer = -0.5 * ( tr(S)/ve - crossprod(e)/(ve^2) - 1/ve * (q-tr(C22)/(vu))  )
ve + FirDer / SecDer
ve

# Direct calculation
sqrt(crossprod(u)/tr(P%*%Vi))
vu
sqrt(crossprod(e)/tr(P%*%I))
ve

# Half derivation
e = Y - yHat
t(e/ve) %*% Z %*% u / tr(P%*%Vi)
vu
e = Y - yHat
t(e/ve) %*% e / tr(P%*%I)
ve

# No-P No-C solvers for vu
          
# Odd solver 1
sqrt(crossprod(u)/crossprod(crossprod(Z,e))*(ve^2))
# Odd solver 2
sqrt(crossprod(u)/(crossprod(Z%*%u/vu,e/ve)))
# Odd solver 3
crossprod(u)/(crossprod(Z%*%u,e/ve))

          
          
