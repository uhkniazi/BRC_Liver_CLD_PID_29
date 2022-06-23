# File: 04_gseaAllContrasts.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: gene set enrichment analysis for the datasets
# Date: 23/6/2022


## set variables and source libraries
## libraries to load
library(gage)

lFiles = list.files('results/', pattern='DEAnalysis*', full.names = T, ignore.case = T)

ldfData = lapply(lFiles, function(x) as.data.frame(read.csv(x, header=T, row.names=1, stringsAsFactors = F)))
names(ldfData) = lFiles
sapply(ldfData, nrow)

# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

sapply(ldfData, function(df) identical(rownames(df), rn))

# ## map the ensemble ids to enterez ids
# ## some genes will be missing or duplicated
# # do some acrobatics to match both the tables and names
# library(org.Hs.eg.db)
# df = select(org.Hs.eg.db, keys = rn, columns = 'ENTREZID', keytype = 'ENSEMBL')
# df = na.omit(df)
# table(duplicated(df$ENSEMBL))
# table(duplicated(df$ENTREZID))
# # choose unique ids
# i = match(rn, df$ENSEMBL)
# df = df[i,]
# df = na.omit(df)
# i = match(df$ENSEMBL, rn)
# rn = rn[i]
# identical(rn, df$ENSEMBL)
# dim(df)
# dim(na.omit(df))

# ## match the tables and names
# ldfData = lapply(ldfData, function(df2){
#   df2 = df2[rn,]
#   df2$ENTREZID = df$ENTREZID
#   return(df2)
# })
# 
# # sanity checks
# rn = rownames(ldfData[[1]])
# head(rn)
# rm(df)
# sapply(ldfData, function(df) identical(rownames(df), rn))

cvTitle = gsub('results//DEAnalysis_', '', names(ldfData))
cvTitle = gsub('.xls', '', cvTitle)

## load msig db data
oMsigGS.c2 = readList('dataExternal/msigdb/c2.all.v7.5.1.entrez.gmt')
oMsigGS.c3 = readList('dataExternal/msigdb/c3.all.v7.5.1.entrez.gmt')
oMsigGS.c5 = readList('dataExternal/msigdb/c5.all.v7.5.1.entrez.gmt')
oMsigGS.c7 = readList('dataExternal/msigdb/c7.all.v7.5.1.entrez.gmt')
oMsigGS.c8 = readList('dataExternal/msigdb/c8.all.v7.5.1.entrez.gmt')
oMsigGS.hm = readList('dataExternal/msigdb/h.all.v7.5.1.entrez.gmt')


## choose a contrast to work with loop through
for (i in 1:length(ldfData)){
  dfContrast = ldfData[[i]]
  # for a contrats of choice create the list
  iContFc = dfContrast$logFC
  ## add enterez ids
  names(iContFc) = as.character(dfContrast$ENTREZID)
  head(iContFc)
  head(dfContrast)
  oGage = gage(iContFc, oMsigGS.c2)
  
  dfGreater = data.frame(oGage$greater)
  #str(dfGreater)
  #i = which(dfGreater$p.val < 0.01)
  #rownames(dfGreater[i,])
  
  dfLess = data.frame(oGage$less)
  #str(dfLess)
  #i = which(dfLess$p.val < 0.01)
  #rownames(dfLess[i,])
  
  write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_upregulated_pathways_mSigDb_c2.xls', sep=''))
  write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_downregulated_pathways_mSigDb_c2.xls', sep=''))
  
  ## c3
  oGage = gage(iContFc, oMsigGS.c3)
  
  dfGreater = data.frame(oGage$greater)
  #str(dfGreater)
  #i = which(dfGreater$p.val < 0.01)
  #rownames(dfGreater[i,])
  
  dfLess = data.frame(oGage$less)
  #str(dfLess)
  #i = which(dfLess$p.val < 0.01)
  #rownames(dfLess[i,])
  
  write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_upregulated_pathways_mSigDb_c3.xls', sep=''))
  write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_downregulated_pathways_mSigDb_c3.xls', sep=''))
  
  
  ## c5
  oGage = gage(iContFc, oMsigGS.c5)
  
  dfGreater = data.frame(oGage$greater)
  dfLess = data.frame(oGage$less)
  
  write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_upregulated_pathways_mSigDb_c5.xls', sep=''))
  write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_downregulated_pathways_mSigDb_c5.xls', sep=''))

  ## c7
  oGage = gage(iContFc, oMsigGS.c7)
  
  dfGreater = data.frame(oGage$greater)
  dfLess = data.frame(oGage$less)
  
  write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_upregulated_pathways_mSigDb_c7.xls', sep=''))
  write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_downregulated_pathways_mSigDb_c7.xls', sep=''))
  
  ## c8
  oGage = gage(iContFc, oMsigGS.c8)
  
  dfGreater = data.frame(oGage$greater)
  dfLess = data.frame(oGage$less)
  
  write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_upregulated_pathways_mSigDb_c8.xls', sep=''))
  write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_downregulated_pathways_mSigDb_c8.xls', sep=''))
  
  ## hm
  oGage = gage(iContFc, oMsigGS.hm)
  
  dfGreater = data.frame(oGage$greater)
  dfLess = data.frame(oGage$less)
  
  write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_upregulated_pathways_mSigDb_hm.xls', sep=''))
  write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file=paste('results/gsea/', cvTitle[i], '_downregulated_pathways_mSigDb_hm.xls', sep=''))
  
}
