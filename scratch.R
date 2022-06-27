# scratch.R
# https://github.com/uhkniazi/BRC_CirrhosisBiomarker_McPhailM/blob/85df3c8ef110e2ae8f15847d2e43b7c8ba048deb/02_variableSelectionLNEG.R#L121
# https://github.com/uhkniazi/BRC_CSF_Proteomics_Jamie_PID_15/blob/bf8c7acad352667a474af7aa7e5018fc546e70b5/03_deAnalysis.R#L111
# https://github.com/uhkniazi/BRC_FOG1KO_John_PID_25/blob/8f9a034eceb6399760db58759f45321338bff737/10_deAnalysis.R#L168
# https://github.com/uhkniazi/BRC_Keloid/blob/57032328576d24948a611f5c6fe1ee7193cd5b05/Keloid_main/de_4_contrasts.R#L146

##################################################
library(UpSetR)
mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")

df = read.csv(file.choose(), header=T, stringsAsFactors = F)
upset(df, sets = colnames(df)[2:5], sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = NULL)

stanDso = rstan::stan_model(file='tResponse1RandomEffectsMultipleScales.stan')


### setup data object for stan and run model
## subset the data to get the second level of nested parameters
## this is done to avoid loops in the stan script to map the scale parameters
## of each ind/gene to the respective set of coefficients for jitters
d = dfData[!duplicated(dfData$Coef), ]
#d2 = dfData[!duplicated(dfData$Coef.2), ]

lStanData = list(Ntotal=nrow(dfData), 
                 Nclusters1=nlevels(dfData$Coef),
                 #Nclusters2=nlevels(dfData$Coef.2),
                 NScaleBatches1 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 #NScaleBatches2 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NgroupMap1=as.numeric(dfData$Coef),
                 #NgroupMap2=as.numeric(dfData$Coef.2),
                 NBatchMap1=as.numeric(d$ind), # this is where we use the second level mapping
                 #NBatchMap2=as.numeric(d2$ind), # this is where we use the second level mapping
                 Nnu=nlevels(dfData$ind),
                 NsigmaPop=nlevels(dfData$ind),
                 NnuMap=as.numeric(dfData$ind),
                 NsigmaPopMap=as.numeric(dfData$ind),
                 y=dfData$values, 
                 intercept = mean(dfData$values), intercept_sd= sd(dfData$values)*3)

initf = function(chain_id = 1) {
  list(sigmaRan1 = rep(1, times=lStanData$NScaleBatches1),
       #sigmaRan2= rep(0.1, times=lStanData$NScaleBatches2),
       rGroupsJitter1 = rep(0, times=lStanData$Nclusters1))
  #rGroupsJitter2_scaled = rep(0, times=lStanData$Nclusters2),
  #phi_scaled=rep(15, times=lStanData$Nphi))
}


ptm = proc.time()

fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2,
                    pars=c('sigmaRan1',
                           #'sigmaRan2',
                           'nu',
                           #'mu',
                           'rGroupsJitter1',
                           #'rGroupsJitter2',
                           'betas',
                           'sigmaPop'
                    ),
                    cores=2)#, init=initf)#, control=list(adapt_delta=0.99, max_treedepth = 11), init=initf)
#save(fit.stan, file='results/fit.stan.nb_3Mar.rds')
ptm.end = proc.time()
print(fit.stan)


m = model.matrix(values ~ 1 + age + Coef.1 + Coef.2 + Coef.3, data = dfData)
c = coef(fit.lm)
c[2] = 2
dfData.old = dfData
dfData$values = m %*% c + rnorm(10)
