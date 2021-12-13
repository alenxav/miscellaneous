
# GWAS for multiple populations based on the MAGIC model (missing data friendly)
ExtFlexGwas = function(y,gen,fam,rescale=TRUE){
  tag = 0
  require(lme4)
  cat('Running marker: ')
  sma = apply(gen,2,function(x,y,fam,rescale){
    # Print marker under evaluation
    tag <<- tag+1
    cat(paste(tag,'. ',sep=''))
    TmpDta = data.frame(y=y,f=factor(fam),x=x)
    lvl = levels(TmpDta$f)
    TmpDta = droplevels.data.frame(TmpDta[rowMeans(!is.na(TmpDta))==1,])
    if(rescale){TmpDta$y = scale(TmpDta$y); TmpDta$x = scale(TmpDta$x)}
    Vy = c(var(TmpDta$y))
    # Null model
    fit0 = lmer(y~(1|f),TmpDta)
    ll0 = logLik(fit0)
    ## Model 1 - Within-family effect
    fit1 = suppressMessages(lmer(y~(1|f)+(1|f):x,TmpDta))
    eff1 = ranef(fit1)$f[,2]
    names(eff1) = rownames(ranef(fit1)$f)
    eff1 = eff1[lvl]
    ll1 = logLik(fit1)
    LRT1 = ll1-ll0
    PVAL1 = -log10(1-pchisq(LRT1,1))
    ## Model 2 - Across-family effect
    fit2 = suppressMessages(lmer(y~x+(1|f),TmpDta))
    eff2 = fit2@beta[2]
    ll2 = logLik(fit2)
    LRT2 = ll2-ll0
    PVAL2 = -log10(1-pchisq(LRT2,1))
    ## Coeff of determination
    R2 = 1-c(fit0@devcomp$cmp['sigmaREML'],
             fit1@devcomp$cmp['sigmaREML'],
             fit2@devcomp$cmp['sigmaREML'])/Vy
    names(R2)=paste0('R2.model',0:2)
    ## Output
    NumObs = nrow(TmpDta)
    out = c(NumObs,P1=PVAL1,P2=PVAL2,R2,Fxd=eff2,eff1)
    return(out)},fam=fam,y=y,rescale=rescale)
  cat('\n')
  rownames(sma) = c('NumObs','PVAL1','PVAL2','NHR2','WFR2','AFR2','EffOverall',
                    paste('Eff',sort(unique(fam)),sep=''))
  sma[is.na(sma)] = 0
  return(data.frame(t(sma)))
}


