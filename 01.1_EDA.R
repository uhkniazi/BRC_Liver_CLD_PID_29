# File: 01.1_EDA.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 24/5/2022
# Desc: use the metadata and matrix for EDA

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
rn = as.character(dfSample$description)
rn = strsplit(rn, ';')
rn = sapply(rn, function(x) return(x[2]))
identical(colnames(mData), rn)
rownames(dfSample) = rn

# mCounts = mData
# dim(mCounts)
# dfSample$fAge = cut(dfSample$Age, 4, include.lowest = T, labels = 1:4)

str(dfSample)

xtabs( ~ dfSample$group1 + dfSample$Gender)
xtabs( ~ dfSample$group1 + dfSample$Ethnicity)
xtabs( ~ dfSample$group1 + dfSample$fAge)
xtabs( ~ dfSample$Ethnicity + dfSample$Gender)

iSub = which(dfSample$group1 %in% 'CLD-BA')

mData = mData[,iSub]
dfSample = dfSample[iSub,]

# iSub = which(dfSample$group2 %in% 'muscle')
# 
# mData = mData[,-iSub]
# dfSample = dfSample[-iSub,]
dfSample = droplevels.data.frame(dfSample)
str(dfSample)

dim(mData)
iGapdh = '16747338'
iSrsf4 = '16684222'
iSel = '16946825'
iOthers = rownames(mData)[1:5000]
mData.sub = mData[c(iGapdh, iSrsf4, iSel, iOthers), ]
dim(mData.sub)
t = dfSample$group2

# x = mData['16946825',]
# h = colMeans(mData[c(iSrsf4,iGapdh),])
# df = data.frame(cbind(x, h))
# df = stack(df)
# df$tissue = dfSample$group2
# df$pid = dfSample$group3
library(lme4)
# df$f = df$ind:df$tissue
# f = lmer(values ~ 1 + (1 | f), data=df)
# summary(f)
# coef(f)
# d = as.data.frame(VarCorr(f))

v = rep(NA, times=5003)
for (i in 1:5003){
  f = lmer(mData.sub[i,] ~ 1 + (1 | t))
  d = as.data.frame(VarCorr(f))
  v[i] = d[d$grp == 't', 'sdcor']
}
names(v) = rownames(mData.sub)
sort(round(v, 3), decreasing = T)[1:30]

df = dfAnnotation[dfAnnotation$ID %in% rownames(mData.sub), ]
dim(df)
i = match(names(v), df$ID)
df = df[i,]
identical(names(v), as.character(df$ID))
dfExport = data.frame(ID=as.character(names(v)), SD=round(v, 3), matchedID=as.character(df$ID), REFSEQ=as.character(df$GB_ACC), stringsAsFactors = F)
df = select(x, keys = dfExport$REFSEQ, columns = 'ENTREZID', keytype = 'REFSEQ')
dfExport$ENTREZID = df$ENTREZID
write.csv(dfExport, file='temp/varianceCheck.csv')

mCounts = sweep(log(mData), 2, log(mData[iSrsf4,]), '/')
range(mCounts)
## some EDA diagnostic plots on the data matrix

library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

colnames(mCounts) = as.character(dfSample$group1)
oDiag.1 = CDiagnosticPlots(mCounts, 'cld-ba')

# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
str(dfSample)
fBatch = factor(dfSample$group2)
levels(fBatch)

# choose a different grouping variable
# summary(dfMeta$CGAdelivery)
# fBatch = cut(dfMeta$CGAdelivery, 4, include.lowest = T)
# table(fBatch)
## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F
oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)

boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.5)
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.5)
plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.5, cex.main=1)
plot.PCA(oDiag.1, fBatch, cex.main=1, legend.pos = 'topright')#, csLabels = as.character(dfMeta$fGroups))
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)

## extract the PCA components and model the variation
######## modelling of PCA components to assign sources of variance to covariates in the design
par(p.old)
plot(oDiag.1@lData$PCA$sdev)
# use the first 3 principal components
mPC = oDiag.1@lData$PCA$x[,1:2]

## try a linear mixed effect model to account for varince
library(lme4)
# prepare data for input
dfData = data.frame(mPC)
dfData = stack(dfData)
str(dfData)
dfData$values = as.numeric(scale(dfData$values))

library(lattice)
densityplot(~ values, data=dfData)
densityplot(~ values | ind, data=dfData, scales=list(relation='free'))

# add covariates of interest to the data frame
dfData$fTreatment = factor(dfSample$group1[iSub])
dfData$fGender = dfSample$Gender[iSub]
dfData$fEthnicity = dfSample$Ethnicity[iSub]
dfData$fAge = dfSample$fAge[iSub]
dfData$age = dfSample$Age[iSub]
str(dfData)

densityplot(~ values | ind, groups=fTreatment, data=dfData, auto.key = list(columns=3), scales=list(relation='free'))
densityplot(~ values | ind, groups=fGender, data=dfData, auto.key = list(columns=3), scales=list(relation='free'))
densityplot(~ values | ind, groups=fEthnicity, data=dfData, auto.key = list(columns=3), scales=list(relation='free'))
densityplot(~ values | ind, groups=fAge, data=dfData, auto.key = list(columns=4))

# format data for modelling, i.e. create coefficients to estimate
str(dfData)
dfData$Coef.1 = factor(dfData$fTreatment:dfData$ind)
dfData$Coef.2 = factor(dfData$fGender:dfData$ind)
dfData$Coef.3 = factor(dfData$fEthnicity:dfData$ind)
dfData$Coef.4 = factor(dfData$fAge:dfData$ind)
str(dfData)

fit.lme1 = lmer(values ~ 1  + (1 | Coef.1), data=dfData)
fit.lme2 = lmer(values ~ 1 + (1 | Coef.1) + (1 | Coef.2), data=dfData)
fit.lme3 = lmer(values ~ 1  + (1 | Coef.1) + (1 | Coef.3), data=dfData)
fit.lme4 = lmer(values ~ 1 + age + (1 | Coef.1) + (1 | Coef.2) + (1 | Coef.3), data=dfData)
fit.lme5 = lmer(values ~ 1  + age + (1 | Coef.1), data=dfData)

summary(fit.lme1)


anova(fit.lme1, fit.lme5)

# summary(fit.lme1)
# summary(fit.lme2)
# 
# plot((fitted(fit.lme2)), resid(fit.lme2), pch=20, cex=0.7)
# lines(lowess((fitted(fit.lme2)), resid(fit.lme2)), col=2)
# hist(dfData$values, prob=T)
# lines(density(fitted(fit.lme2)))

# these db packages interfere with stan so unload 
detach('package:org.Hs.eg.db', unload=T)
## fit model with stan with various model sizes
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(rethinking)

stanDso = rstan::stan_model(file='tResponsePartialPooling.stan')

######## models of various sizes using stan
dfData = df
dfData$values = log(dfData$values)
str(dfData)
m1 = model.matrix(values ~ f - 1, data=dfData)
m2 = model.matrix(values ~ pid - 1, data=dfData)
m3 = model.matrix(values ~ Coef.3 - 1, data=dfData)
#m4 = model.matrix(values ~ age - 1, data=dfData)

m = cbind(m1, m2)#, m4)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$f)),
                                              rep(2, times=nlevels(dfData$pid))),
                                              # rep(3, times=nlevels(dfData$Coef.3))),
                                              #rep(1, times=1)),
                 y=dfData$values)

fit.stan.4 = sampling(stanDso, data=lStanData, iter=2000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
                                                                           'nu', 'mu', 'log_lik'),
                      cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
print(fit.stan.4, c('populationMean', 'sigmaPop', 'sigmaRan', 'nu', 'betas'), digits=3)

traceplot(fit.stan.4, 'populationMean')
traceplot(fit.stan.4, 'sigmaPop')
traceplot(fit.stan.4, 'sigmaRan')

## extract coefficients of interest to plot
mCoef = extract(fit.stan.4)$betas
dim(mCoef)
dim(m)
colnames(mCoef) = colnames(m)
colnames(mCoef)
muscle = mCoef[,'fx:muscle'] - mCoef[,'fh:muscle']
liver = mCoef[,'fx:liver'] - mCoef[,'fh:liver']
fat = mCoef[,'fx:fat'] - mCoef[,'fh:fat']
ml = mCoef[,'fx:muscle'] - mCoef[,'fx:liver']


mCor = (cor(mCoef[,grep('PC2', colnames(mCoef))]))
dim(mCor)
colnames(mCor) = gsub('Coef.', '', colnames(mCor))
rownames(mCor) = gsub('Coef.', '', rownames(mCor))
heatmap(abs(mCor), Rowv = NA, Colv = NA, symm = T, scale = 'none', cexRow = 0.8, cexCol = 0.8)

i = which(colnames(mCoef) %in% c("Coef.1CN:PC1", "Coef.1CN:PC2", "Coef.1CN:PC3"))
i2 = which(colnames(mCoef) %in% c("Coef.3Middle_East:PC1", "Coef.3Middle_East:PC2", "Coef.3Middle_East:PC3"))
colnames(mCoef)[c(i, i2)]
pairs(mCoef[,c(i, i2)], pch=1, cex=0.2)
round(cor(mCoef[,c(i, i2)]), 2)

# # remove outliers
# mPlot = scale(mCoef[,c(i, i2)])
# mPlot[abs(mPlot) > 2.5] = NA
# mPlot = apply(mPlot, 2, function(x) sample(x, 300, replace = T))
# dim(mPlot)
# pairs(mPlot, pch=1, cex=0.2)
# 
# #i = grep('PC1', colnames(mCoef))
# mPlot = scale(mCoef)
# mPlot = apply(mPlot, 2, function(x) sample(x, 500, replace = F))
# dim(mPlot)
# mPlot[abs(mPlot) > 4] = NA
# colnames(mPlot) = gsub('Coef.|:PC1', '', colnames(mPlot))
# pairs(mPlot[,c(5, 6, 7, 9, 10, 11)], pch=1, cex=0.2, col='grey')
# pairs(mPlot[,c(3, 6, 7, 9, 10, 11)], pch=1, cex=0.2, col='grey')
# 
# #pairs(mCoef[,c(i)], pch=1, cex=0.2)

m = extract(fit.stan.4)$sigmaRan
dim(m)
colnames(m) = c('T', 'G', 'E')
pairs(log(m), pch=1, cex=0.2)


### model with only treatment and gender
m = cbind(m1, m2)
lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$Coef.1)),
                                              rep(2, times=nlevels(dfData$Coef.2))),
                 y=dfData$values)

fit.stan.2 = sampling(stanDso, data=lStanData, iter=2000, chains=2, pars=c('betas', 'populationMean', 'sigmaPop', 'sigmaRan',
                                                                           'nu', 'mu', 'log_lik'),
                      cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
print(fit.stan.2, c('populationMean', 'sigmaPop', 'sigmaRan', 'nu', 'betas'), digits=3)

traceplot(fit.stan.2, 'populationMean')
traceplot(fit.stan.2, 'sigmaPop')
traceplot(fit.stan.2, 'sigmaRan')

## extract coefficients of interest to plot
mCoef = extract(fit.stan.2)$betas
dim(mCoef)
dim(m)
colnames(mCoef) = colnames(m)
mCor = (cor(mCoef[,grep('PC2', colnames(mCoef))]))
dim(mCor)
colnames(mCor) = gsub('Coef.', '', colnames(mCor))
rownames(mCor) = gsub('Coef.', '', rownames(mCor))
heatmap(abs(mCor), Rowv = NA, Colv = NA, symm = T, scale = 'none', cexRow = 0.8, cexCol = 0.8)

pairs(mCoef[,c('Coef.1CN:PC1', 'Coef.1CLD-BC:PC1', 'Coef.1CLD-BA:PC1', 'Coef.2F:PC1', 'Coef.2M:PC1' )], pch=1, cex=0.2)

## some model scores and comparisons
compare(fit.stan.4, fit.stan.2)
compare(fit.stan.4, fit.stan.2, func = LOO)
plot(compare(fit.stan.4, fit.stan.2))

############### new simulated data
###############
### generate some posterior predictive data
## generate random samples from alternative t-distribution parameterization
## see https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
## follow the algorithm in section 14.3 page 363 in Gelman 2013
simulateOne = function(mu, sigma, nu){
  yrep = rt_ls(length(mu), nu, mu,  sigma)
  return(yrep)
}

## sample n values, numerous times
mDraws.sim = matrix(NA, nrow = nrow(dfData), ncol=300)
l = extract(fit.stan.2)
for (i in 1:300){
  p = sample(1:nrow(l$mu), 1)
  mDraws.sim[,i] = simulateOne(l$mu[p,], 
                               l$sigmaPop[p],
                               l$nu[p])
}

dim(mDraws.sim)
plot(density(dfData$values), main='posterior predictive density plots, model 4')
apply(mDraws.sim, 2, function(x) lines(density(x), lwd=0.5, col='lightgrey'))
lines(density(dfData$values))

## plot residuals
plot(dfData$values - colMeans(l$mu) ~ colMeans(l$mu))
lines(lowess(colMeans(l$mu), dfData$values - colMeans(l$mu)))
apply(l$mu[sample(1:nrow(l$mu), 100),], 1, function(x) {
  lines(lowess(x, dfData$values - x), lwd=0.5, col=2)
})

## plot the original PCA and replicated data
plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=c(1:5)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], main='PCA Components - original and simulated',
     xlab='PC1', ylab='PC2')
points(rowMeans(mDraws.sim)[dfData$ind == 'PC1'], rowMeans(mDraws.sim)[dfData$ind == 'PC2'],
       col=c(1:5)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], pch='1')

plot(dfData$values[dfData$ind == 'PC1'], dfData$values[dfData$ind == 'PC2'], 
     col=c(1:5)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], main='PCA Components - original and model 2',
     xlab='PC1', ylab='PC2', xlim=c(-3, 3), ylim=c(-3, 3), pch=15)

apply(mDraws.sim, 2, function(x) {
  points(x[dfData$ind == 'PC1'], x[dfData$ind == 'PC2'],
         col=c(1:5)[as.numeric(dfData$fTreatment[dfData$ind == 'PC1'])], pch=20)
})


m = cbind(extract(fit.stan.4)$sigmaRan, extract(fit.stan.4)$sigmaPop) 
dim(m)
m = log(m)
colnames(m) = c('Treatment', 'Gender', 'Ethnicity', 'Residual')
pairs(m, pch=20, cex=0.5, col='grey')

df = stack(data.frame(m[,-4]))
histogram(~ values | ind, data=df, xlab='Log SD', scales=list(relation='free'))

## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}
## define some test quantities to measure the lack of fit
## define a test quantity T(y, theta)
## variance
T1_var = function(Y) return(var(Y))

## min quantity
T1_min = function(Y){
  return(min(Y))
} 

## max quantity
T1_max = function(Y){
  return(max(Y))
} 

## mean quantity
T1_mean = function(Y){
  return(mean(Y))
} 

## mChecks
ivResp = dfData$values
mChecks = matrix(NA, nrow=4, ncol=1)
rownames(mChecks) = c('Variance', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('model 1')

t1 = apply(mDraws.sim, 2, T1_var)
mChecks['Variance', 1] = getPValue(t1, var(ivResp))

## testing for outlier detection i.e. the minimum value show in the histograms earlier
t1 = apply(mDraws.sim, 2, T1_min)
t2 = T1_min(ivResp)
mChecks['Min',1] = getPValue(t1, t2)

## maximum value
t1 = apply(mDraws.sim, 2, T1_max)
t2 = T1_max(ivResp)
mChecks['Max', 1] = getPValue(t1, t2)

## mean value
t1 = apply(mDraws.sim, 2, T1_mean)
t2 = T1_mean(ivResp)
mChecks['Mean', 1] = getPValue(t1, t2)

mChecks
