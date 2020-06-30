library(AlphaSimR)

# Selection intesity
Intensity = 0.10

# Founder settings
founderPop <- runMacs(nInd = 100,
                      nChr = 20,
                      segSites = 3000,
                      inbred = TRUE,
                      species = "GENERIC")

# Simulation parameters
SP = SimParam$new(founderPop)
SP$addSnpChip(nSnpPerChr=50)
SP$addTraitAE(nQtlPerChr=20,mean=50,var=30,relAA=0.5)
SP$setVarE(varE=30) #, h2=0.25)

# Settings of runs
Number_of_runs = 20
Number_of_generations = 40
Plot_at_the_end_of_each_run = TRUE

# Replicate loop
RES = list()
for(run in 1:Number_of_runs){
  
  # Print current run
  cat('\n\n RUN',run,'\n\n')
  
  # Reset pop
  pop = newPop(founderPop, simParam=SP)
  pop = randCross(pop, nCrosses=20, nProgeny = 10, simParam=SP)
  pop = makeDH(pop, nDH=1, simParam=SP)
  
  # Store parameters
  genMean = c()
  genVar = c()
  H2 = c()
  
  # Loop across generations
  for(i in 1:Number_of_generations){
    
    # Set (make it reproducible)
    set.seed(i)
    
    # Print generation
    cat('Generation',i,'\n')
    
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
    
    # Add EBVs to the pop
    pop = setEBV(pop, ans, simParam=SP)
    
    # Store population mean and var
    genMean = c(genMean, meanG(pop))
    genVar = c(genVar, varG(pop))
    h2_tmp = varG(pop)/varP(pop)
    H2 = c(H2, h2_tmp)
    cat(' --Mu =',round(meanG(pop),2),'\n')
    cat(' --Var =',round(varG(pop),2),'\n')
    cat(' --H2 =',round(h2_tmp,2),'\n')
    
  }
  
  # Check run
  if(Plot_at_the_end_of_each_run){
    par(mfrow = c(1,3),mar=c(4,2,2,1),cex=1)
    plot(genMean,main="Population Mean",xlab='Generation',ylab='',col=2,type='o',pch=20)
    plot(genVar,main="Genetic variance",xlab='Generation',ylab='',col=3,type='o',pch=20)
    plot(H2,main="Heritability",xlab='Generation',ylab='',col=4,type='o',pch=20)
  }
  
  # Store run
  RES[[run]] = list(Mu=genMean,GV=genVar,H2=H2)
  
}


#######################
# Final summary graph #
######################

par(mfrow = c(1,3),mar=c(4,4,2,1),cex=1)
MU = t(sapply(RES,function(x)x$Mu))
boxplot(MU,main="Population Mean",xlab='Generation',ylab='Phenotypic scale',col='pink',border='darkred',type='o',pch=20)
GV = t(sapply(RES,function(x)x$GV))
boxplot(GV,main="Genetic variance",xlab='Generation',ylab='Variance',col='lightgreen',border='darkgreen',type='o',pch=20)
H2 = t(sapply(RES,function(x)x$H2))
boxplot(H2,main="Heritability",xlab='Generation',ylab='Heritability',col='lightblue',border='darkblue',type='o',pch=20)
