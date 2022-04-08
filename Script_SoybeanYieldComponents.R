
## THIS IS THE SCRIPT UTILIZED TO ANALYZE SOYBEAN YIELD COMPONENTS
## THE PAPER IS PUBLISHED HERE https://doi.org/10.1534/g3.119.400896

# Load packages
require(lme4)
require(NAM)

# Get data
data(soyin,package = 'SoyNAM')
data = rbind(data.line.in,data.check.in)
gen = gen.in

# QC phenotypes
data$yield[data$yield>95] = NA
data$yield[data$yield<15] = NA
data$Pods[data$Pods>70] = NA
data$Nodes[data$Nodes>20] = NA
data$Pods.Nodes[data$Pods.Nodes>5] = NA
data$Pods.Nodes[data$Pods.Nodes<1] = NA

# Spatial covariate
data = data.frame(data,spPod=SPC(data$Pods,data$BLOCK,data$ROW,data$COL,2,2) )
data = data.frame(data,spNodes=SPC(data$Nodes,data$BLOCK,data$ROW,data$COL,2,2) )
data = data.frame(data,spPN=SPC(data$Pods.Nodes,data$BLOCK,data$ROW,data$COL,2,2) )
data = data.frame(data,spYLD=SPC(data$yield,data$BLOCK,data$ROW,data$COL,2,2) )

# QC genotypic files on MAF
gen = snpQC(gen, MAF=0.05)

# Family and SNPs/Chromosome
fam = as.numeric(gsub("...$|DS1.-", "", rownames(gen)))
chr = c(table(substr(colnames(gen),3,4)))
cc = (as.numeric(substr(colnames(gen),3,4))%%2+1)*2

# Fit mixed model
fit1 = lmer(formula = Pods~spPod+(1|environ)+(1|strain), data=data) 
fit2 = lmer(formula = Nodes~spNodes+(1|environ)+(1|strain), data=data) 
fit3 = lmer(formula = Pods.Nodes~spPN+(1|environ)+(1|strain), data=data) 
fit4 = lmer(formula = yield~spPN+(1|environ)+(1|strain), data=data) 

# Get genetic values from mixed model
blup1 = ranef(fit1)$strain[rownames(gen),1]+mean(data$Pods,na.rm=T) # Pods
blup2 = ranef(fit2)$strain[rownames(gen),1]+mean(data$Nodes,na.rm=T) # Nodes
blup3 = ranef(fit3)$strain[rownames(gen),1]+mean(data$Pods.Nodes,na.rm=T) # PN
blup4 = ranef(fit4)$strain[rownames(gen),1]+mean(data$yield,na.rm=T) # Yield

# Broad-sense heritability
H1 = VarCorr(fit1)$strain[1,1]/(VarCorr(fit1)$strain[1,1]+(getME(fit1,"devcomp")$cmp['sigmaREML']^2)/2.3)
H2 = VarCorr(fit2)$strain[1,1]/(VarCorr(fit2)$strain[1,1]+(getME(fit2,"devcomp")$cmp['sigmaREML']^2)/2.3)
H3 = VarCorr(fit3)$strain[1,1]/(VarCorr(fit3)$strain[1,1]+(getME(fit3,"devcomp")$cmp['sigmaREML']^2)/2.3)
H4 = VarCorr(fit4)$strain[1,1]/(VarCorr(fit4)$strain[1,1]+(getME(fit3,"devcomp")$cmp['sigmaREML']^2)/2.3)

########################
# ASSOCIATION ANALYSIS #
########################

# QC for repeated individuals
QCing = cleanREP(y=cbind(blup1,blup2,blup3,blup4),gen,fam,thr=0.95)
QCing$y = IMP(QCing$y)
blup1 = QCing$y[,1]
blup2 = QCing$y[,2]
blup3 = QCing$y[,3]
blup4 = QCing$y[,3]
Y = QCing$y
gen = QCing$gen
fam = QCing$fam

# GWAS 1
ex = eigX(gen,fam)
am1 = gwas2(blup1,gen,fam,chr,EIG=ex)
am2 = gwas2(blup2,gen,fam,chr,EIG=ex)
am3 = gwas2(blup3,gen,fam,chr,EIG=ex)
rm(ex)

# GWAS 2
bb1 = BGLR::BGLR(blup1,nIter = 20000, burnIn=1000, thin =1, ETA=list(list(X=gen,model='BayesB')), verbose=F); Pbb1=(-log(1-bb1$ETA[[1]]$d))
bb2 = BGLR::BGLR(blup2,nIter = 20000, burnIn=1000, thin =1, ETA=list(list(X=gen,model='BayesB')), verbose=F); Pbb2=(-log(1-bb2$ETA[[1]]$d))
bb3 = BGLR::BGLR(blup3,nIter = 20000, burnIn=1000, thin =1, ETA=list(list(X=gen,model='BayesB')), verbose=F); Pbb3=(-log(1-bb3$ETA[[1]]$d))

# GWAS 3
rf1 = ranger::ranger(y~., data=data.frame(y=blup1,gen), importance='impurity', num.trees = 10000)$variable.importance
rf2 = ranger::ranger(y~., data=data.frame(y=blup2,gen), importance='impurity', num.trees = 10000)$variable.importance
rf3 = ranger::ranger(y~., data=data.frame(y=blup3,gen), importance='impurity', num.trees = 10000)$variable.importance

# Permutaion threshold for RF
thr1 = c(); for(i in 1:1000){cat(i,' ');thr1 = c(thr1, max(ranger::ranger(y~.,data=data.frame(y=sample(blup1),gen), importance='impurity')$variable.importance) ); cat(round(mean(thr1),2),'\n')}; Thr1 = quantile(thr1,.95); 
thr2 = c(); for(i in 1:1000){cat(i,' ');thr2 = c(thr2, max(ranger::ranger(y~.,data=data.frame(y=sample(blup2),gen), importance='impurity')$variable.importance) ); cat(round(mean(thr2),2),'\n')}; Thr2 = quantile(thr2,.95); 
thr3 = c(); for(i in 1:1000){cat(i,' ');thr3 = c(thr3, max(ranger::ranger(y~.,data=data.frame(y=sample(blup3),gen), importance='impurity')$variable.importance) ); cat(round(mean(thr3),2),'\n')}; Thr3 = quantile(thr3,.95); 

#####################
# CROSS VALIDATIONS #
#####################

# CV1 - Leave family out
CV1.1 = bWGR::emCV(blup1,gen,llo=fam,avg=FALSE)
CV1.2 = bWGR::emCV(blup2,gen,llo=fam,avg=FALSE)
CV1.3 = bWGR::emCV(blup3,gen,llo=fam,avg=FALSE)

# CV2 - Across family k=5, 20x
CV2.1 = bWGR::emCV(blup1,gen,5,20,avg=FALSE)
CV2.2 = bWGR::emCV(blup2,gen,5,20,avg=FALSE)
CV2.3 = bWGR::emCV(blup3,gen,5,20,avg=FALSE)

# CV3 - Within family k=5, 20x
require(snow)
cl <- makeCluster(15, type = "SOCK")
CV3.1 = parSapply(X=unique(fam),y=blup1,FUN=function(i,y,gen,fam){w=which(fam==i);return(bWGR::emCV(y[w],gen[w,],5,20)[c("emRR","emEN","emBL","emDE","emBA","emBB","emBC","emML","emMX")])},gen=gen,fam=fam,cl=cl)
CV3.2 = parSapply(X=unique(fam),y=blup2,FUN=function(i,y,gen,fam){w=which(fam==i);return(bWGR::emCV(y[w],gen[w,],5,20)[c("emRR","emEN","emBL","emDE","emBA","emBB","emBC","emML","emMX")])},gen=gen,fam=fam,cl=cl)
CV3.3 = parSapply(X=unique(fam),y=blup3,FUN=function(i,y,gen,fam){w=which(fam==i);return(bWGR::emCV(y[w],gen[w,],5,20)[c("emRR","emEN","emBL","emDE","emBA","emBB","emBC","emML","emMX")])},gen=gen,fam=fam,cl=cl)
stopCluster(cl)

# mcCV1 - Leave family out
mcCV1.1 = bWGR::mcmcCV(blup1,gen,llo=fam,avg=FALSE)
mcCV1.2 = bWGR::mcmcCV(blup2,gen,llo=fam,avg=FALSE)
mcCV1.3 = bWGR::mcmcCV(blup3,gen,llo=fam,avg=FALSE)

# mcCV2 - Across family k=5, 20x
mcCV2.1 = bWGR::mcmcCV(blup1,gen,5,20,avg=FALSE)
mcCV2.2 = bWGR::mcmcCV(blup2,gen,5,20,avg=FALSE)
mcCV2.3 = bWGR::mcmcCV(blup3,gen,5,20,avg=FALSE)

# mcCV3 - Within family k=5, 20x
require(snow)
cl <- makeCluster(15, type = "SOCK")
mcCV3.1 = parSapply(X=unique(fam),y=blup1,FUN=function(i,y,gen,fam){w=which(fam==i);return(bWGR::mcmcCV(y[w],gen[w,],5,20)[c("BayesA","BayesB","BayesC","BayesL","BayesCpi","BayesDpi","BayesRR")])},gen=gen,fam=fam,cl=cl)
mcCV3.2 = parSapply(X=unique(fam),y=blup2,FUN=function(i,y,gen,fam){w=which(fam==i);return(bWGR::mcmcCV(y[w],gen[w,],5,20)[c("BayesA","BayesB","BayesC","BayesL","BayesCpi","BayesDpi","BayesRR")])},gen=gen,fam=fam,cl=cl)
mcCV3.3 = parSapply(X=unique(fam),y=blup3,FUN=function(i,y,gen,fam){w=which(fam==i);return(bWGR::mcmcCV(y[w],gen[w,],5,20)[c("BayesA","BayesB","BayesC","BayesL","BayesCpi","BayesDpi","BayesRR")])},gen=gen,fam=fam,cl=cl)
stopCluster(cl)

###################
# MANHATTAN PLOTS #
###################

par(mfcol=c(3,3),mar=c(4,4,1,1))
plot(Pbb1,pch=20,col=cc,ylab='-log(p)',xaxt='n',xlab=''); abline(h=1.3,lty=2)
plot(rf1,pch=20,col=cc,ylab='Importance',xaxt='n',xlab=''); abline(h=Thr1,lty=2)
plot(am1,pch=20,alpha = 0.05/4240);
plot(Pbb2,pch=20,col=cc,ylab='-log(p)',xaxt='n',xlab=''); abline(h=1.3,lty=2)
plot(rf2,pch=20,col=cc,ylab='Importance',xaxt='n',xlab=''); abline(h=Thr2,lty=2)
plot(am2,pch=20,alpha = 0.05/4240);
plot(Pbb3,pch=20,col=cc,ylab='-log(p)',xaxt='n',xlab=''); abline(h=1.3,lty=2)
plot(rf3,pch=20,col=cc,ylab='Importance',xaxt='n',xlab=''); abline(h=Thr3,lty=2)
plot(am3,pch=20,alpha = 0.05/4240);

# r max
pltS=function(X) boxplot(t(apply(X,1,function(x)(x-min(x))/(max(x)-min(x)))))
rS=function(X) sort(colMeans(t(apply(X,1,function(x)(x-min(x))/(max(x)-min(x))))),T)
pltS(CV1.1)
rS(CV1.1)

###############################
# VARIANCE COMPONENT ANALYSIS #
###############################

vc = c()
require(BGLR)
for(i in unique(fam)){
  
  w = which(fam==i)
  cat('Fam',i,length(w),'\n')
  
  yA = blup1[w] # pods
  yB = blup2[w] # nodes
  yC = blup3[w] # pn
  G = CNT(gen[w,])
  KA = tcrossprod(G)
  KAA = KA*KA;
  KA = KA/mean(diag(KA))
  KAA = KAA/mean(diag(KAA))
  
  fA = BGLR(yA,ETA=list(A=list(K=KA,model='RKHS'),AA=list(K=KAA,model='RKHS')),verbose = F,nIter = 5000)
  fB = BGLR(yB,ETA=list(A=list(K=KA,model='RKHS'),AA=list(K=KAA,model='RKHS')),verbose = F,nIter = 5000)
  fC = BGLR(yC,ETA=list(A=list(K=KA,model='RKHS'),AA=list(K=KAA,model='RKHS')),verbose = F,nIter = 5000)
  
  Pod = c(Va=fA$ETA$A$varU,Vaa=fA$ETA$AA$varU,Ve=fA$varE); Pod = Pod/sum(Pod)
  Node = c(Va=fB$ETA$A$varU,Vaa=fB$ETA$AA$varU,Ve=fB$varE); Node = Node/sum(Node)
  PN = c(Va=fC$ETA$A$varU,Vaa=fC$ETA$AA$varU,Ve=fC$varE); PN = PN/sum(PN)
  
  vc0 = c(Pod=Pod,Node=Node,PN=PN)
  vc = rbind(vc,vc0)
  
}
rownames(vc) = gsub('1\\.|099','',paste('Fam',unique(fam)/100+1.00099,sep=''))
write.csv(vc,'VarComp.csv',quote=F)
barplot(t(vc[order(vc[,3]),1:3]),las=2,col=2:4)

# PLOT
par(mfrow=c(1,3),mar=c(4,4.5,2,1))
barplot(t(vc[order(vc[,3]),1:3]),las=2,col=c(1,8,0),border=0);legend('topleft','Pod number',bty='n',cex=1.5)
barplot(t(vc[order(vc[,6]),4:6]),las=2,col=c(1,8,0),border=0);legend('topleft','Node number',bty='n',cex=1.5)
barplot(t(vc[order(vc[,9]),7:9]),las=2,col=c(1,8,0),border=0);legend('topleft','Pods per node',bty='n',cex=1.5)
legend('topright',c('Additive variance','Epistatic variance','Residual variance'),bty='n',cex=2,col=c(1,8,1),pch=c(16,16,1))

# Summary
range(vc0$Pod.Vaa+vc0$Pod.Va)
range(vc0$Node.Vaa+vc0$Node.Va)
range(vc0$PN.Vaa+vc0$PN.Va)

range(vc0$Pod.Va)
range(vc0$Node.Va)
range(vc0$PN.Va)

range(vc0$Pod.Vaa)
range(vc0$Node.Vaa)
range(vc0$PN.Vaa)

mean(vc0$Pod.Va)
mean(vc0$Node.Va)
mean(vc0$PN.Va)

mean(vc0$Pod.Vaa)
mean(vc0$Node.Vaa)
mean(vc0$PN.Vaa)

####################
# PHENOTYPIC TREND #
####################

# Scatter plot
if()
library(psych)
dt = data[,c(16,19,18,20)]
names(dt) = c('Yield','Number of pods','Number of nodes','Pods per node')
pairs.panels(dt, 
             method = "pearson",
             hist.col = 0,
             density = TRUE,
             ellipses = TRUE)

# Box plots
par(mfrow=c(3,1),mar=c(4,4.5,1,1))
boxplot(data$Pods~data$family,ylab='Count',xlab='Family',pch=20);legend('topleft','Number of pods',bty='n')
boxplot(data$Nodes~data$family,ylab='Count',xlab='Family',pch=20);legend('topleft','Number of nodes',bty='n')
boxplot(data$Pods.Nodes~data$family,ylab='Ratio',xlab='Family',pch=20);legend('topleft','Pods per node',bty='n')

# Genetic correlations
G = GRM(gen)
mv = bWGR::mkr(Y,gen)
mv$GC

