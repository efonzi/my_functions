#--------------------------------------
# DEPENDENCIES
#--------------------------------------
source("https://www.bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)


#--------------------------------------
# create objects for annotation
#--------------------------------------
a = org.Hs.egSYMBOL2EG
symbol2entrez = as.list(a[mappedkeys(a)])

a = org.Hs.egCHR
entrez2chr= as.list(a[mappedkeys(a)])


#--------------------------------------
# WD and files preparation
#--------------------------------------
wd = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/aml/1612_chromosome_annotation"
SCRIPTPATH = "/Users/Eugenio/Desktop/Lavoro/Seragnoli/scripts"
setwd(wd)

tableName="chr_genes_table.txt"
table_pref = strsplit(tableName, "[.]")[[1]][1]
system(paste0("sh ", SCRIPTPATH, "/1612_chr_annot_aml.sh < ", wd, "/", tableName, " > ", wd, "/", table_pref, "_ed.txt"))
tableName = paste0(table_pref, "_ed.txt")

fileName = "gain_dup.txt"
#fileName = "loss_del.txt"
file_pref = strsplit(fileName, "[.]")[[1]][1]
system(paste0("dos2unix -c Mac ", fileName))

#-------------------------------------------------------------------------------------------------------
# TAKES GENE SYMBOLS LIST AND ANNOTATES THE CHR NUMS VIA THE ENTREZ ID, BY MEANS OF org.Hs.eg.db -
# ENTRIES WITH DOUBLE ENTREZ ID, NULL ENTRIES AND ENTRIES WITH DOUBLE CHR ARE SAVED IN SEPARATE VECTORS
#-------------------------------------------------------------------------------------------------------
df = read.table(fileName, header = T, stringsAsFactors = F)
df$chr = rep(".", dim(df)[1])

nulls = vector()
doubleEntr = vector()
XY = vector()

for(i in 1:dim(df)[1]){
  
  entrez = symbol2entrez[[as.character(df[i,1])]]
  if(is.null(entrez) == TRUE){
    df[i,2] = "NULL"
    nulls = append(nulls, as.character(df[i,1]))
  }else{
    if(length(entrez) > 1){
      df[i,2] = "2entr"
      doubleEntr = append(doubleEntr, as.character(df[i,1]))
    }else{
      entrez2 = as.character(entrez2chr[[entrez]])
      if(length(entrez2) > 1){
        df[i,2] = entrez2[1]
        XY = append(XY, as.character(df[i,1]))
      }else{
        df[i,2] = entrez2
      }
    }
  }
}
table(df$chr)
nulls = sort(nulls)
doubleEntr = sort(doubleEntr)
XY = sort(XY)


#-------------------------------------------------------------------------------------------------------
# PARSES DATA FRAME OF GENE SYMBOLS CREATED IN LAST STEP FOR CHR=NULL AND CORRECTS THEM ACCORDING TO THE 
# INFO IN "TABLE" - THEN WRITES TO FILE
#-------------------------------------------------------------------------------------------------------

table = unique(read.table(tableName, header = T, sep="\t", stringsAsFactors = F))
sel = which(is.na(table[,2]))
table = table[-sel,]


for(i in 1:dim(df)[1]){
  if(df$chr[i] == "NULL"){
    for(j in 1:dim(table)[1]){
      if(df[,1][i] == table$GENE[j]){
        if(table$CHR[j] == 23){df$chr[i] = "X"}
        if(table$CHR[j] == 24){df$chr[i] = "Y"}
        else{df$chr[i] = table$CHR[j]}
      }
    }
  }
}
write.table(df, file=paste0(file_pref, "_annot.txt"), quote=F, row.names=F, col.names=T, sep="\t")

#-----------------------------------------------------------------
# CREATES LIST OF GENE SYMBOLS WITH CHR=NULL AND WRITES IT TO FILE
#-----------------------------------------------------------------

nulls = vector()
for(i in 1:dim(df)[1]){
  if(df$chr[i] == "NULL"){
    nulls = append(nulls, df[,1][i])
  }
}

write.table(nulls, file=paste0(file_pref, "_nulls.txt"), quote=F, row.names=F, col.names=F, sep="\t")

table(df$chr)
