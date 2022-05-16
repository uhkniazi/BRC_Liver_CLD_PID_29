# File: 02_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 13/5/2022
# Desc: use the GEO data matrix and metadata to find DEGs

source('header.R')

library(RMySQL)
##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

dbGetQuery(db, paste('describe Sample;'))
q = paste0('select Sample.* from Sample where Sample.idData = 52')
dfSample = dbGetQuery(db, q)
dim(dfSample)
head(dfSample)
# close connection after getting data
dbDisconnect(db)

## load the downloaded data and extract matrix
library(GEOquery)
gse =  getGEO(filename = 'dataExternal/GSE84954_series_matrix.txt.gz')
mData = exprs(gse)
dim(mData)
rn = dfSample$description
rn = strsplit(rn, ';')
rn = sapply(rn, function(x) return(x[2]))
identical(colnames(mData), rn)
rownames(dfSample) = rn

## select one tissue
dfSample$description[1]
table(dfSample$group2)
i = which(dfSample$group2 == 'liver')

dfSample = dfSample[i,]
mData = mData[,i]
identical(colnames(mData), rownames(dfSample))

### set up data for modelling with stan
mData.orig = mData
mData = t(as.matrix(mData.orig[1,]))
dim(mData)

plot(density(mData))
range(mData)

## prepare covariates for model matrix
dfData = data.frame(t(mData))
dim(dfData); dim(dfSample);
dfData = stack(dfData)
dim(dfData)
head(dfData)
table(dfSample$group1); table(dfSample$group2); table(dfSample$group3)
f = factor(dfSample$group1)
levels(f)

dfData$fTreatment = f
dfData$Coef = factor(dfData$fTreatment:dfData$ind)
dim(dfData)
dfData = droplevels.data.frame(dfData)
#dfData = dfData[order(dfData$Coef, dfData$Coef.adj1, dfData$Coef.adj2), ]
str(dfData)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

fit.stan = sampling(stanDso, data=lStanData, iter=500, chains=2,
                    pars=c('sigmaRan1',
                           #'sigmaRan2',
                           'nu',
                           'mu',
                           'rGroupsJitter1',
                           #'rGroupsJitter2',
                           'betas',
                           'sigmaPop'
                    ),
                    cores=2)#, control=list(adapt_delta=0.99, max_treedepth = 11), init=initf)
#save(fit.stan, file='results/fit.stan.nb_3Mar.rds')
ptm.end = proc.time()