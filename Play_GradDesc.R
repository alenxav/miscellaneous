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

par(mfrow=c(3,3),mar=c(1,1,2,1),cex=0.8)
for(i in 1:45){
  e = y-Zb
  d = c(-2*crossprod(e,Z)+2*lmb*b)/length(y)
  b = b - d*a
  Zb = c(Z%*%b)
  cat('Cor',round(cor(Zb,h$hat),4),'\n')
  if(i %% 5 == 0){
    p = paste('Cor',round(cor(Zb,h$hat),8),'\n')
    plot(b,h$b,pch=20,main=paste('Iteration',i),xaxt='n',yaxt='n')
    legend('topleft',legend = p,bty='n')
  }
}
