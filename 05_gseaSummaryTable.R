# File: 05_gseaSummaryTable.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: merge the gsea results for all contrasts in one table
# Date: 23/6/2022


lFiles = list.files('results/gsea/', pattern='*pathways_mSigDb_c2.xls', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('results/gsea//(.+VS.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('results/gsea//(.+VS.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c2 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c2)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c2)[o], ncol = 2, byrow = T)
colnames(mMerged.c2)[o]
mMerged.c2 = mMerged.c2[,o]

# remove na sections
dim(mMerged.c2)
mMerged.c2 = na.omit(mMerged.c2)
dim(mMerged.c2)
head(mMerged.c2)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c2.bin = getBinaryMatrix(mMerged.c2)

## group this matrix into combinations
mMerged.c2.bin.grp = mMerged.c2.bin
set.seed(123)
dm = dist(mMerged.c2.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c2.bin.grp = cbind(mMerged.c2.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c2.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c2.bin)
dfMerged.c2 = data.frame(round(mMerged.c2, 3), sig.pvals, groups, DB='C2')
str(dfMerged.c2)
head(dfMerged.c2)
tail(dfMerged.c2)


############# C3
lFiles = list.files('results/gsea/', pattern='*pathways_mSigDb_c3.xls', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('results/gsea//(.+VS.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('results/gsea//(.+VS.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c3 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c3)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c3)[o], ncol = 2, byrow = T)
colnames(mMerged.c3)[o]
mMerged.c3 = mMerged.c3[,o]

# remove na sections
dim(mMerged.c3)
mMerged.c3 = na.omit(mMerged.c3)
dim(mMerged.c3)
head(mMerged.c3)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c3.bin = getBinaryMatrix(mMerged.c3)

## group this matrix into combinations
mMerged.c3.bin.grp = mMerged.c3.bin
set.seed(123)
dm = dist(mMerged.c3.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c3.bin.grp = cbind(mMerged.c3.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c3.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c3.bin)
dfMerged.c3 = data.frame(round(mMerged.c3, 3), sig.pvals, groups, DB='C3')
str(dfMerged.c3)
head(dfMerged.c3)
tail(dfMerged.c3)
#############

##################### C5
lFiles = list.files('results/gsea/', pattern='*pathways_mSigDb_c5.xls', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('results/gsea//(.+VS.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('results/gsea//(.+VS.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c5 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c5)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c5)[o], ncol = 2, byrow = T)
colnames(mMerged.c5)[o]
mMerged.c5 = mMerged.c5[,o]

# remove na sections
dim(mMerged.c5)
mMerged.c5 = na.omit(mMerged.c5)
dim(mMerged.c5)
head(mMerged.c5)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c5.bin = getBinaryMatrix(mMerged.c5)

## group this matrix into combinations
mMerged.c5.bin.grp = mMerged.c5.bin
set.seed(123)
dm = dist(mMerged.c5.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c5.bin.grp = cbind(mMerged.c5.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c5.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c5.bin)
dfMerged.c5 = data.frame(round(mMerged.c5, 3), sig.pvals, groups, DB='C5')
str(dfMerged.c5)
head(dfMerged.c5)
tail(dfMerged.c5)
#####################

#################### C7
lFiles = list.files('results/gsea/', pattern='*pathways_mSigDb_c7.xls', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('results/gsea//(.+VS.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('results/gsea//(.+VS.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c7 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c7)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c7)[o], ncol = 2, byrow = T)
colnames(mMerged.c7)[o]
mMerged.c7 = mMerged.c7[,o]

# remove na sections
dim(mMerged.c7)
mMerged.c7 = na.omit(mMerged.c7)
dim(mMerged.c7)
head(mMerged.c7)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c7.bin = getBinaryMatrix(mMerged.c7)

## group this matrix into combinations
mMerged.c7.bin.grp = mMerged.c7.bin
set.seed(123)
dm = dist(mMerged.c7.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c7.bin.grp = cbind(mMerged.c7.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c7.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c7.bin)
dfMerged.c7 = data.frame(round(mMerged.c7, 3), sig.pvals, groups, DB='C7')
str(dfMerged.c7)
head(dfMerged.c7)
tail(dfMerged.c7)
####################

#################### C8
lFiles = list.files('results/gsea/', pattern='*pathways_mSigDb_c8.xls', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('results/gsea//(.+VS.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('results/gsea//(.+VS.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c8 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c8)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.c8)[o], ncol = 2, byrow = T)
colnames(mMerged.c8)[o]
mMerged.c8 = mMerged.c8[,o]

# remove na sections
dim(mMerged.c8)
mMerged.c8 = na.omit(mMerged.c8)
dim(mMerged.c8)
head(mMerged.c8)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c8.bin = getBinaryMatrix(mMerged.c8)

## group this matrix into combinations
mMerged.c8.bin.grp = mMerged.c8.bin
set.seed(123)
dm = dist(mMerged.c8.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c8.bin.grp = cbind(mMerged.c8.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c8.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c8.bin)
dfMerged.c8 = data.frame(round(mMerged.c8, 3), sig.pvals, groups, DB='C8')
str(dfMerged.c8)
head(dfMerged.c8)
tail(dfMerged.c8)
####################

##################### HM
lFiles = list.files('results/gsea/', pattern='*pathways_mSigDb_hm.xls', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('results/gsea//(.+VS.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('results/gsea//(.+VS.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.hm = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.hm)
o = c(1, 3, 2, 4)
## sanity check
matrix(colnames(mMerged.hm)[o], ncol = 2, byrow = T)
colnames(mMerged.hm)[o]
mMerged.hm = mMerged.hm[,o]

# remove na sections
dim(mMerged.hm)
mMerged.hm = na.omit(mMerged.hm)
dim(mMerged.hm)
head(mMerged.hm)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.hm.bin = getBinaryMatrix(mMerged.hm)

## group this matrix into combinations
mMerged.hm.bin.grp = mMerged.hm.bin
set.seed(123)
dm = dist(mMerged.hm.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.hm.bin.grp = cbind(mMerged.hm.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.hm.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.hm.bin)
dfMerged.hm = data.frame(round(mMerged.hm, 3), sig.pvals, groups, DB='HM')
str(dfMerged.hm)
head(dfMerged.hm)
tail(dfMerged.hm)
#####################

ncol(dfMerged.c2)
length(Reduce(intersect, list(colnames(dfMerged.c2), colnames(dfMerged.c3), colnames(dfMerged.c5),
          colnames(dfMerged.c7), colnames(dfMerged.c8), colnames(dfMerged.hm))))

write.csv(dfMerged.c2, file='results/gsea/gsea_msigdb_c2_merged.xls')
write.csv(dfMerged.c3, file='results/gsea/gsea_msigdb_c3_merged.xls')
write.csv(dfMerged.c5, file='results/gsea/gsea_msigdb_c5_merged.xls')
write.csv(dfMerged.c7, file='results/gsea/gsea_msigdb_c7_merged.xls')
write.csv(dfMerged.c8, file='results/gsea/gsea_msigdb_c8_merged.xls')
write.csv(dfMerged.hm, file='results/gsea/gsea_msigdb_hm_merged.xls')

## merge together into one dataframe
# drop the group with most zeros
table(dfMerged.c2$groups)
t = rowSums(mMerged.c2.bin)
table(t, dfMerged.c2$groups)
dfMerged.c2.sub = dfMerged.c2[dfMerged.c2$groups != 3,]

table(dfMerged.c3$groups)
dfMerged.c3.sub = dfMerged.c3[dfMerged.c3$groups != 3,]

table(dfMerged.c5$groups)
dfMerged.c5.sub = dfMerged.c5[dfMerged.c5$groups != 4,]

table(dfMerged.c7$groups)
dfMerged.c7.sub = dfMerged.c7[dfMerged.c7$groups != 4,]

table(dfMerged.c8$groups)
dfMerged.c8.sub = dfMerged.c8[dfMerged.c8$groups != 4,]

table(dfMerged.hm$groups)
t = rowSums(mMerged.hm.bin)
table(t, dfMerged.hm$groups)
dfMerged.hm.sub = dfMerged.hm[dfMerged.hm$groups != 3,]

dfMerged = rbind(dfMerged.c2.sub, dfMerged.c3.sub,
                 dfMerged.c5.sub, dfMerged.c7.sub,
                 dfMerged.c8.sub, dfMerged.hm.sub)
dfMerged = droplevels.data.frame(dfMerged)
dim(dfMerged)
str(dfMerged)
write.csv(dfMerged, file='results/gsea_msigdb_significant_merged.xls')

### heatmaps
### just for a quick visual check, do not use for results
df = dfMerged
head(df)
dim(df)
mMat = as.matrix(df[,c(1:4)])
head(mMat)
mMat = -1*log(mMat+1e-16)
g1 = df[,'groups']
g1 = factor(as.character(g1))
levels(g1)
g2 = df[,'DB']
g2 = factor(as.character(g2))
levels(g2)

ann = data.frame(DB=g2, DBG=g2:g1)#, Group=g1 )
str(ann)
#ann = data.frame(Group=g1 )
range(mMat)
quantile(as.vector(mMat), 0:20/20)
mMat[mMat > 10] = 10 
#mMat[mMat > 8] = 8

library(NMF)
library(RColorBrewer)

i = which(g2 == 'C5')
table(g2)
aheatmap(mMat, annRow = ann$DB, scale = 'none', Rowv = order(g2), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

pdf('dataID50/results/gsea_msigdb_significant_merged.pdf')
# aheatmap(mMat, annRow = NA, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
#          col=c('white', brewer.pal(9, 'YlOrRd')))
aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

dev.off(dev.cur())

############### choose a subset of the data for heatmaps

pdf('results/figures/gsea_msigdb_significant_merged_kegg.pdf')
i = grep('KEGG', rownames(mMat))
m = mMat[i,]
rownames(m) = gsub('KEGG_', '', rownames(m))
aheatmap(m, annRow = NA, scale = 'none', Rowv = NA, Colv=NA, cexRow=1, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))
dev.off(dev.cur())

pdf('results/figures/gsea_msigdb_significant_merged_C8.pdf')
#i = grep('GOBP', rownames(mMat))
i = which(g2 == 'C2')
m = mMat[i,]
dim(m)
#rownames(m) = gsub('REACTOME_', '', rownames(m))
aheatmap(m, annRow = NA, scale = 'none', Rowv = NA, Colv=NA, cexRow=1, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))
dev.off(dev.cur())
