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

## load additional metadata 
df = read.csv(file.choose(), header=T)
identical(dfSample$id, df$id)
dfSample = df
rm(df)

## load the downloaded data and extract matrix
library(GEOquery)
gse =  getGEO(filename = 'dataExternal/GSE84954_series_matrix.txt.gz')
## load the annotation
library(org.Hs.eg.db)
x = org.Hs.eg.db
gpl = getGEO(annotation(gse))
dfAnnotation = gpl@dataTable@table
table(dfAnnotation$ID %in% as.numeric(rownames(gse)))
## select the probes that have a gene id
i = which(dfAnnotation$GB_ACC %in% '')
length(i)
dfAnnotation = dfAnnotation[-i,]
dim(dfAnnotation)
df = select(x, keys = dfAnnotation$GB_ACC, columns = 'ENTREZID', keytype = 'REFSEQ')
## some probes have no corresponding REFSEQ to ENTREZ entry, remove those
table(is.na(df$ENTREZID))
df = na.omit(df)
dim(df)
## some genes are annotated multiple times having multiple probes
table(duplicated(df$REFSEQ))
## don't remove these duplicated probes for this analysis
i = dfAnnotation$GB_ACC %in% df$REFSEQ
table(i)
dfAnnotation = dfAnnotation[i,]
## remove the non-matching probes from the expression matrix
i = rownames(gse) %in% as.character(dfAnnotation$ID)
table(i)
gse = gse[i,]
dim(gse)

mData = exprs(gse)
dim(mData)
plot(density(mData))
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
#i = which(rownames(mData.orig) %in% c('16727967', '16681593', '16972339'))
cvSel = '17006392'#scan(what=character())
i = which(rownames(mData.orig) %in% cvSel)
length(i)
mData = (as.matrix(mData.orig[i,]))
dim(mData)
dim(dfSample)

plot(density(mData))
range(mData)

## prepare covariates for model matrix
dfData = data.frame((mData))
dim(dfData); dim(dfSample);
dfData = stack(dfData)
dim(dfData)
head(dfData)
table(dfSample$group1); table(dfSample$group2); table(dfSample$group3)
f = factor(dfSample$group1)
levels(f)

dfData$fTreatment = f
dfData$fGender = factor(dfSample$Gender)
dfData$fEthnicity = factor(dfSample$Ethnicity)
dfData$age = (dfSample$Age)
dfData$Coef.1 = factor(dfData$fTreatment:dfData$ind)
dfData$Coef.2 = factor(dfData$fGender:dfData$ind)
dfData$Coef.3 = factor(dfData$fEthnicity:dfData$ind)
dfData$Coef.4 = factor(dfData$ind)
dim(dfData)
dfData = droplevels.data.frame(dfData)
#dfData = dfData[order(dfData$Coef, dfData$Coef.adj1, dfData$Coef.adj2), ]
str(dfData)

# these db packages interfere with stan so unload 
detach('package:org.Hs.eg.db', unload=T)
## setup the stan model
library(Rcpp)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='tResponse3RandomEffects1SlopeMultipleScales.stan')


### setup data object for stan and run model
## subset the data to get the second level of nested parameters
## this is done to avoid loops in the stan script to map the scale parameters
## of each ind/gene to the respective set of coefficients for jitters
d = dfData[!duplicated(dfData$Coef.1), ]
d2 = dfData[!duplicated(dfData$Coef.2), ]
d3 = dfData[!duplicated(dfData$Coef.3), ]
d4 = dfData[!duplicated(dfData$Coef.4), ]

lStanData = list(Ntotal=nrow(dfData), 
                 Nclusters1=nlevels(dfData$Coef.1),
                 Nclusters2=nlevels(dfData$Coef.2),
                 Nclusters3=nlevels(dfData$Coef.3),
                 Nclusters4=nlevels(dfData$Coef.4),
                 NScaleBatches1 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NScaleBatches2 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NScaleBatches3 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NScaleBatches4 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NgroupMap1=as.numeric(dfData$Coef.1),
                 NgroupMap2=as.numeric(dfData$Coef.2),
                 NgroupMap3=as.numeric(dfData$Coef.3),
                 NgroupMap4=as.numeric(dfData$Coef.4),
                 NBatchMap1=as.numeric(d$ind), # this is where we use the second level mapping
                 NBatchMap2=as.numeric(d2$ind), # this is where we use the second level mapping
                 NBatchMap3=as.numeric(d3$ind), # this is where we use the second level mapping
                 NBatchMap4=as.numeric(d4$ind), # this is where we use the second level mapping
                 Nnu=nlevels(dfData$ind),
                 NsigmaPop=nlevels(dfData$ind),
                 NnuMap=as.numeric(dfData$ind),
                 NsigmaPopMap=as.numeric(dfData$ind),
                 y=as.numeric(dfData$values),
                 X = as.numeric(dfData$age),
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
                           'sigmaRan2',
                           'sigmaRan3',
                           'nu',
                           #'mu',
                           'rGroupsJitter1',
                           'rGroupsJitter2',
                           'rGroupsJitter3',
                           'rGroupsSlopes',
                           'betas',
                           'slope',
                           'sigmaPop'
                    ),
                    cores=2)#, init=initf)#, control=list(adapt_delta=0.99, max_treedepth = 11), init=initf)
#save(fit.stan, file='results/fit.stan.nb_3Mar.rds')
ptm.end = proc.time()
print(fit.stan)
## compare with lme4
library(lme4)
fit.lme = lmer(values ~ 1 + age + (1 | Coef.1) + (1 | Coef.2) + (1 | Coef.3), data=dfData)
summary(fit.lme)

## extract results
## get the coefficient of interest - Coef.1 (treatment) in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# # ## get the intercept at population level
iIntercept = as.numeric(extract(fit.stan)$betas)
# ## add the intercept to each random effect variable, to get the full coefficient
mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
# the split is done below on : symbol, but factor name has a : symbol due
# to creation of interaction earlier, do some acrobatics to sort that issue
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
# d$`1` = d$`1`:d$`2`
# d = d[,-4]
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
str(d)
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='CN', deflection='CLD-BA') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
#dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

