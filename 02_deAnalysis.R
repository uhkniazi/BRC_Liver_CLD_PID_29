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
mData = mData.orig[1:3,]
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
