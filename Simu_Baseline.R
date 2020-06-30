library(AlphaSimR)
library(ggplot2)
library(bWGR)

# Set (make it reproducible)
set.seed(12345)

# Selection intesity
Intensity = 0.2
GS_model = BayesA

# Founder settings
founderPop <- runMacs(nInd = 100,
                      nChr = 20,
                      segSites = 3000,
                      inbred = TRUE,
                      species = "GENERIC")

# Simulation parameters
SP = SimParam$new(founderPop)
SP$addSnpChip(nSnpPerChr=50)
SP$addTraitAE(nQtlPerChr=20,mean=50,var=30,relAA=0.25)
SP$setVarE(h2=0.25)

# Create a pop
pop = newPop(founderPop, simParam=SP)
pop = randCross(pop, nCrosses=20, nProgeny = 10, simParam=SP)
pop = makeDH(pop, nDH=1, simParam=SP)

# Object to store population means (track gain)
genMean = c()
genVar = c()

# Loop across generations
for(i in 1:20){
  
  # Store population mean and var
  genMean = c(genMean, meanG(pop))
  genVar = c(genVar, varG(pop))

  # Print generation
  cat('Generation',i,'\n')
  cat(' --Mu =',round(meanG(pop),2),'\n')
  cat(' --Var =',round(varG(pop),2),'\n')
  
  #######
  # SSD #
  #######
  
  # F1
  if(i==1){ # No GS in the first round
    pop = selectCross(pop, nInd=pop@nInd, nCrosses=50, use="pheno", nProgeny=1, simParam=SP)
  }else{ # Genomic Selection
    pop = selectCross(pop, nInd=round(pop@nInd*Intensity), use="ebv", nCrosses=50, nProgeny=1, simParam=SP)
  }
  # F2
  pop = self(pop, nProgeny=10, simParam=SP)
  # F2.3
  pop = self(pop, nProgeny=1, simParam=SP)
  # F2.4
  pop = self(pop, nProgeny=1, simParam=SP)
  
  #############
  # Genotypes #
  #############
  
  # get genotypes
  gen = pullSnpGeno(pop, simParam = SP)
  
  # fit built-in WGP model
  ans = RRBLUP(pop, simParam=SP)
  
  # fit using bWGR package
  fit = GS_model(pop@pheno[,1],gen)

  # replacing original parameters by bWGR's
  ans@bv[[1]]@addEff = fit$b
  ans@bv[[1]]@intercept = fit$mu-mean(pop@pheno[,1])
  ans@Vu[1,1] = (var(pop@pheno[,1])-fit$ve)/sum(apply(gen,2,var))
  ans@Ve[1,1] = fit$ve
  
  # Add EBVs to the pop
  pop = setEBV(pop, ans, simParam=SP)

}

# CHECK THE POPULATION DEVELOPMENT AFTER 20 CYCLES
par(mfrow = c(1,2))
plot(genMean,main="Population Mean",xlab='Generation',ylab='Yield',col=2,type='o',pch=20)
plot(genVar,main="Genetic variance",xlab='Generation',ylab='Variance',col=3,type='o',pch=20)

