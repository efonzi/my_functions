
library(org.Hs.eg.db)
library(reactome.db)
library(GO.db)
library(KEGG.db)
library(biomaRt)

# REACTOME PA - CLUSTER PROFILER
# Reactome total number of pathways (gene sets) = 
# Reactome total number of genes (universe) = 6748
# KEGG total number of Pathways (gene sets) = 
# KEGG total number of genes (universe) = 7216
# GOBP total number of Pathways (gene sets) = 
# GOBP total number of genes (universe) = 16672

# get OrgDb object for H.Sapiens
hub <- AnnotationHub()
hsapiens = query(hub, c("OrgDb", "Sapiens"))

######## use bioMart to get all genes in a GO term #######
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations, create a Mart object
# convert goterm to goid
goid = mapIds(GO.db, keys=goterm, column="GOID", keytype="TERM") # get GOid from GOterm

#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),filters = 'go_id', values = 'GO:0007507', mart = ensembl)
# "attributes" will be the columns of the resulting table
# "hgnc_symbol" indicates the gene symbols


# annotate chromosome number from gene symbol
a = org.Hs.egSYMBOL2EG
symbol2entrez = as.list(a[mappedkeys(a)])

a = org.Hs.egCHR
entrez2chr= as.list(a[mappedkeys(a)])


rownames(df[df[,2] == ".",])
df[df[,1] == "MEMO1",]
df[3055,1]
symbol2entrez[["CDCA8"]]
entrez2chr[["51072"]]

