require(NAM)
data(tpod)
Z = CNT(gen)
rm(gen,fam,chr)
y = y-mean(y)
h = emML(y,Z) # Baseline = EMREML
lmb = h$Ve/h$Vb
b = rep(0,ncol(Z))
Zb = c(Z%*%b)
e = y-Zb
a = 1/ncol(Z)
plot(Zb,y)

for(i in 1:50){
  e = y-Zb
  d = c(-2*crossprod(e,Z)+2*lmb*b)/length(y)
  b = b - d*a
  Zb = c(Z%*%b)
  cat('Cor',round(cor(Zb,h$hat),4),'\n')
}

plot(b)
