# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 12/5/2022


## set variables and source libraries
source('header.R')

library(GEOquery)
## open the soft format and raw data
gse =  getGEO(filename = 'dataExternal/GSE84954_series_matrix.txt.gz')

# get the samples from the expression set object
dfSamples = pData(gse)
identical(rownames(dfSamples), colnames(exprs(gse)))
dim(dfSamples)

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe MetaFile;'))$Field[-1]

## grouping factors in the sample data frame
fDisease = gse$`disease:ch1`
fSubject = gse$`subjectid:ch1`
fTissue = gse$`tissue:ch1`

table(fDisease)
d = gsub('chronic liver disease', 'CLD', fDisease)
d = gsub('Crigler-Najjar', 'CN', d)
fDisease = d

table(fSubject)
d = gsub('(^\\d+)', 'S\\1', fSubject)
fSubject = d

table(fTissue)
fTissue[fTissue == 'muscle (rectus abdominis)'] = 'muscle' 
fTissue[fTissue == 'subcutaneous fat'] = 'fat' 

xtabs( ~ fDisease + fSubject + fTissue)

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfSamples$title,
                       location='downloaded from GEO GSE84954_series_matrix.txt.gz',
                       description= paste('GEO Sample', as.character(rownames(dfSamples)),
                                          'group1 is Disease',
                                          'group2 is Tissue',
                                          'group3 is Subject ID',
                                          sep=';'),
                       group1 = fDisease,
                       group2= fTissue,
                       group3= fSubject)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 52;'))

dbDisconnect(db)
