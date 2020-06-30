
# Simulator
library(AlphaSimR)

# Prediction models
library(bWGR)
RF = function(y,gen) list(hat=ranger::ranger(y~.,data.frame(y=y,gen),mtry=ncol(gen)*0.5)$predictions)
RKHS = function(y,gen){ K = GAU(gen); diag(K)=diag(K)+0.00001; E = eigen(K,symmetric = T); fit = emML(y,E$vectors,E$values); return(fit)} 
XGB = function(y,gen){require(xgboost); X = as(data.matrix(gen), "dgCMatrix"); fit0 = xgboost(data=X,label=y,params=list(subsample=0.25),nrounds=20,objective="reg:squarederror"); return(list(hat=predict(fit0,X)))}
phenotypic = function(y,gen) return(list(hat=y))
#Rcpp::sourceCpp('dnn_chol_norm.cpp'); DNN = function(y,gen) FNN(as.matrix(y),gen) # source('kerasNN.R')

# Final list of models
GS_model = list(
  #'DeepNet'=DNN,
  'Pheno'=phenotypic,
  'GBLUP'=emML,
  'BayesB'=emBB,
  'FLM'=emDE,
  'XGBoost'=XGB,
  'RForest'=RF,
  'RKHS'=RKHS)

# Founder settings
founderPop <- runMacs(nInd = 100,
                      nChr = 20,
                      segSites = 750,
                      inbred = TRUE,
                      species = "GENERIC")

# Simulation parameters
SP = SimParam$new(founderPop)
SP$addSnpChip(nSnpPerChr=100)
SP$addTraitAEG(nQtlPerChr=250,varEnv=10,
               mean=100,var=30,
               varGxE=25,relAA=1)
SP$setVarE(varE=50)

# Settings of runs
Number_of_runs = 10
Number_of_families = 40
Number_of_generations = 40 
Family_size = 10
FamSelIntensity = 0.5
SelIntensity = c(0.1,0.2,0.3)


# Store here
RES = c()

# Loop
for(run in 1:Number_of_runs){
  
  # Set (make it reproducible)
  set.seed(run)
  
  for(Intensity in SelIntensity){
    for(model in names(GS_model)){
      for(Strategy in c('AF','WF','WBF')){
        
        # Print current run
        cat('\n\n',model,Intensity,run,'\n\n')
        
        # Reset pop
        pop = newPop(founderPop, simParam=SP)
        pop = randCross(pop, nCrosses=20, nProgeny = Family_size, simParam=SP)
        pop = makeDH(pop, nDH=1, simParam=SP)
        
        # Loop across generations
        for(i in 1:Number_of_generations){
          
          # Print generation
          cat('Generation',i,'\n')
          
          #######
          # SSD #
          #######
          
          # F1
          if(i==1){ # No GS in the first round
            pop = selectCross(pop, nInd=pop@nInd, nCrosses=Number_of_families, use="pheno",nProgeny=1, simParam=SP)
          }else{
            
            # Across-family
            if(Strategy=='AF'){
              pop = selectCross(pop,nInd=round(pop@nInd*Intensity), use = 'ebv',
                                nCrosses=Number_of_families, nProgeny=1, simParam=SP)
            }
            
            # Within-family
            if(Strategy=='WF'){
              pop = selectWithinFam(pop, nInd=round(Intensity*Family_size),use='ebv')
              pop = selectCross(pop, nInd=pop@nInd, nCrosses=Number_of_families, nProgeny=1, simParam=SP)
            }
            
            # Within-family
            if(Strategy=='WBF'){
              pop = selectFam(pop,nFam=Number_of_families*FamSelIntensity,use="ebv",simParam=SP)
              pop = selectWithinFam(pop, nInd=round(Intensity*Family_size/FamSelIntensity),use='ebv')
              pop = selectCross(pop, nInd=pop@nInd, nCrosses=Number_of_families, nProgeny=1, simParam=SP)
            }
            
            
          }
          
          # F2
          pop = self(pop, nProgeny=Family_size, simParam=SP)
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
          if(!exists('ans')) ans = RRBLUP(pop, simParam=SP)
          
          # fit using models
          if(length(unique(gen%*%rnorm(ncol(gen))))<=5){
            fit = list(hat=rnorm(length(fit$hat)))
          }else{
            fit = GS_model[[model]](pop@pheno[,1],gen)
            if(anyNA(fit$hat)) fit$hat=rnorm(length(fit$hat))
          }
          
          # Add EBVs to the pop
          pop = setEBV(pop, ans, simParam=SP)
          pop@ebv = as.matrix(fit$hat,ncol=1)
          
          #################
          # Store outcome #
          #################
          
          # Store results
          h2_tmp = varG(pop)/varP(pop)
          out = c(
            'Strategy'=Strategy,
            'SelInt'=Intensity,
            'Model'=model,
            'Run'=run,
            'Generation'=i,
            'Accuracy'=cor(pop@gv,pop@ebv),
            'PopMean'=meanG(pop),
            'GenVar'=varG(pop)[1,1],
            'H2'=h2_tmp[1,1])
          RES = rbind(RES,out)
          print(out)
          
        }
      }
      write.csv(RES,'AX_BigLoop_Results.csv',quote=F,row.names=F)
    }
  }
  
  # Analyze
  if(run>1) source('plot.R')
  
}


