##############################################################################################################################
#################################################### FUNCTIONS ###############################################################
##############################################################################################################################

#------------------------------------------------- FOR ARRAY DATA PREPROCESSING ---------------------------------------------------------------

# load a PATIENTxGENE table from file, edit sample names, filter out undesired genes and patients and write to file (the old file will be overwritten)
filterOutGenesAndPatientsChromothripsis <- function(geneTable, sampleID){
  
  sel = grep("LOC", colnames(geneTable))
  if("CLOCK" %in% colnames(geneTable)){
    sel = sel[sel!=grep("CLOCK", colnames(geneTable))]
    if(length(sel)>0) geneTable<-geneTable[, -sel]} 
  else{if(length(sel)>0) geneTable<-geneTable[, -sel]}
  sel = grep("LINC", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("^OR[0-9]", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("MUC", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("KRT", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("RYR", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("CSMD", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("NRXN", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("CNTN", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("[.]AS", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("MIR", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("ANK", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("SNORD", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("[.]IT", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("SCARNA", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  
  sel = which(geneTable$Sample %in% sampleID$Sample)
  if(length(sel)>0) geneTable<-geneTable[sel,]
  
  return(geneTable)
  
}


# load a PATIENTxGENE table from file, filter out undesired genes and write to file (the old file will be overwritten)
filterOutGenesAnna <- function(geneTable, sampleID){
  
  sel = grep("^OR[0-9]", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("MUC", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("KRT", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("-AS", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("TTT", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  sel = grep("-IT", colnames(geneTable))
  if(length(sel)>0) geneTable<-geneTable[, -sel]
  
  # # edit gene names ending with a dot followed by a number
  # sel = grep("[.]", colnames(geneTable))
  # if(length(sel)>0) colnames(geneTable)[sel] = sapply(strsplit(colnames(geneTable)[sel], "[.]"), function(x) x[[1]][1])
  
  return(geneTable)
  
}

# accept two PATIENTxGENE tables, adapt the samples (rows) and genes (columns) of both to table1, merge their content (summing up the count values), write merged table to file
# table1 and table2 do not have to contain the same patients(rows) and genes(columns)
# "event" must be a string with the combination of the two events (es. "gain-dup", "loss-del")
mergePatientGeneTables <- function(table1, table2){
  
  # adds samples (rows) unique of table1 to table2
  sel<-which(table1$Sample %in% table2$Sample)
  toadd<-table1[-sel, 'Sample'] #select the samples in table1 that are not contained in table2
  mtoadd<-data.frame(toadd, matrix(rep(0, length(toadd)*(dim(table2)[2]-1)), ncol=(dim(table2)[2]-1)))  #make dataframe of 0 with rows like toadd and columns like table2
  colnames(mtoadd)<-colnames(table2)
  table2<-rbind(table2,mtoadd)
  
  # adds samples (rows) unique of table2 to table1
  sel<-which(table2$Sample %in% table1$Sample)
  toadd<-table2[-sel, 'Sample'] #select the samples in table2 that are not contained in table1
  mtoadd<-data.frame(toadd, matrix(rep(0, length(toadd)*(dim(table1)[2]-1)), ncol=(dim(table1)[2]-1)))  #make dataframe of 0 with rows like "toadd" and columns like table1
  colnames(mtoadd)<-colnames(table1)
  table1<-rbind(table1,mtoadd)
  
  # merge tables (sum up values in common columns and adds the unique columns in table2 to table1)
  for(i in 2:dim(table2)[2]){
    
    sel<-colnames(table2)[i]
    if(sel %in% colnames(table1)){
      k<-which(colnames(table1)==sel)
      table1[,k]<-table1[,k]+table2[,i]
    }
    else {
      table1<-cbind(table1, table2[,i])
      names(table1)[dim(table1)][2]<-names(table2)[i]
    }
  }
  
  return(table1)
  
} 


#------------------------------------------------- FOR WES DATA PREPROCESSING ---------------------------------------------------------------

# take "wesTable", convert it to PATIENTxGENE table and fill each cell with counts of how 
# many times each genes is altered in each patient
# return PATIENTxGENE table
convertWesTableToPatientGeneTable <- function(wesTable){
  
  # from "wesTable" create PATIENTxGENE matrix ("geneTable") and fill it with zeroes
  patients = sort(unique(wesTable$Sample))
  genes = sort(unique(wesTable$GENE))
  mat = matrix(rep(0, length(patients)*length(genes)), ncol=(length(genes)))
  colnames(mat) = genes
  geneTable = data.frame(Sample = patients, mat, stringsAsFactors = F)
  
  # loop through genes(columns) and patients(rows) of "geneTable" and fill the cells according to how many times each gene is altered in each patient
  for(i in 1:dim(geneTable)[2]){
    gen = colnames(geneTable)[i]
    subset = wesTable[wesTable$GENE == gen,]
    for(j in 1:dim(geneTable)[1]){
      pat = geneTable$Sample[j]
      if(pat %in% subset$Sample){
        pat_tab = table(subset$Sample)
        geneTable[j,i] = as.numeric(pat_tab[pat])}}
  } 
  
  return(geneTable)
  
}


#------------------------------------------------- FOR GENE ENRICHMENT - FISHER -----------------------------------------------------------------

# load file containing sample classification and split it in lists of Cases and Ctrls
separateCaseCtrl <- function(criterion, case, ctrl, sampleID){
  
  # "criterion" must be the column title of the desired sample classification in the file "wesnames.tsv"/"chipnames.tsv"
  # "cases" and "ctrl" can be both integers or strings
  # "platform" must be either "wes" or "array"
  
  criterion = which(colnames(sampleID) == criterion)
  
  selCase<-sampleID[which(sampleID[,criterion]==case & sampleID$Sample!=""),]
  selCtrl<-sampleID[which(sampleID[,criterion]==ctrl & sampleID$Sample!=""),]
  selNA<-sampleID[which(sampleID[,criterion]=="na" & sampleID$Sample!=""),]
  
  return(list(selCase,selCtrl,selNA))
  
}

# take a PATIENTxGENE table, convert values to zeroes and ones and
# perform fisher's test on genes between cases and ctrls
geneCaseCtrlFisher <- function(geneTable, selCase, selCtrl){
  
  # "selCase" and "selCtrl" must be data frames containing a column "Sample" with the sample names of CASES and CTRLS, respectively
  
  # split table in Cases and Ctrls (by row)
  geneCase = geneTable[which(geneTable$Sample %in% selCase$Sample),]
  geneCtrl = geneTable[which(geneTable$Sample %in% selCtrl$Sample),]
  
  # create Fmatrix, with event counts for each gene in cases and ctrls
  caseCounts<-apply(geneCase[,2:dim(geneCase)[2]], 2, sum)
  ctrlCounts<-apply(geneCtrl[,2:dim(geneCtrl)[2]], 2, sum)
  Fmatrix<-data.frame(caseCounts, ctrlCounts)
  
  # create Rmatrix and execute two.sided fisher's test with FDR p-adjustment
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
  colnames(Rmatrix)<-c("GENE", "%Case", "%Ctrl", "PVAL", "QVAL")
  sumCase<-dim(geneCase)[1]
  sumCtrl<-dim(geneCtrl)[1]
  for(i in 1:dim(Fmatrix)[1]){
    if(i%%1000 == 0){cat("gene #", i, "\n")}
    posCase<-Fmatrix[i, 'caseCounts']
    posCtrl<-Fmatrix[i, 'ctrlCounts']
    negCase<-sumCase-posCase
    negCtrl<-sumCtrl-posCtrl
    
    x<-matrix(c(posCase, posCtrl, negCase, negCtrl), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=20)
      
      Rmatrix[i, 2]<-posCase/sumCase*100
      Rmatrix[i, 3]<-posCtrl/sumCtrl*100
      Rmatrix[i, 4]<-pval
    } 
  }
  Rmatrix$GENE<-rownames(Fmatrix)
  Rmatrix$QVAL<-round(p.adjust(Rmatrix$PVAL, method="BH"), digits=20)
  Rmatrix<-Rmatrix[order(Rmatrix$QVAL),]
  
  return(Rmatrix)
}


#------------------------------------------------------ FOR PATHWAY ENRICHMENT ------------------------------------------------------------------

##### OVER-REPRESENTATION TEST

# loop over patients and perform over-representation test for REACTOME, KEGG, GOBP pathways (based on hypergeometric distribution)
# requires a PATIENTxGENE matrix of zeroes and ones "geneTable"
hypergeometricTestOnMultipleIndividuals <- function(geneTable){
  
  reactomeALL=NULL
  keggALL=NULL
  gobpALL=NULL
  
  # loop over table rows (patients)
  for(r in 1:dim(geneTable)[1]){
    
    # get name of the patient
    patient = geneTable[r,1]
    
    # select altered genes for the patient
    genes = colnames(geneTable)[which(geneTable[r,] == 1)]
    
    cat("     \n")
    cat("PATIENT", r, "of", dim(geneTable)[1], " - ", length(genes), "events", " - ", folderName, "\n")
    
    if(length(genes)>0){
      
      cat("pathway over-representation test \n")
      
      # transform gene symbols to entrez ids
      ids = as.vector(unlist(mget(genes,envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
      ids = ids[complete.cases(ids)]
      
      # perform over-representation test for reactome, add sample name, add case/ctrl classification number and
      # concatenate to "reactomeALL"
      if(length(ids)>0){
        if(TRUE %in% (ids %in% keys(reactome.db, keytype="ENTREZID"))){ # to check whether the ids map in reactome
          reactome = try(enrichPathway(ids, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1))
          if(!inherits(reactome, "NULL")){
            reactome = reactome@result
            reactome$Sample = patient
            reactomeALL = rbind(reactomeALL, reactome)
            cat("reactome --> DONE \n")
          }else{cat("SKIP reactome \n")}
        }else{cat("SKIP reactome \n")}}
      
      # perform over-representation test for kegg, add sample name, add case/ctrl classification number and
      # concatenate to "keggALL"
      if(length(ids)>0){
        kegg = try(enrichKEGG(ids, organism="hsa", pAdjustMethod = "BH",  pvalueCutoff = 1, qvalueCutoff = 1))
        if(!inherits(kegg, "NULL")){
          kegg = kegg@result
          kegg$Sample = patient
          keggALL = rbind(keggALL, kegg)
          cat("kegg --> DONE \n")
        }else{cat("SKIP kegg \n")}}
      
      # perform over-representation test for gobp, add sample name, add case/ctrl classification number and
      # concatenate to "gobpALL"
      # try() and inherits() allow to decide whether to skip the rest of the loop in case enrichGO() finds an error and returns NULL
      if(length(ids)>0){
        gobp = try(enrichGO(gene=ids, OrgDb=org.Hs.eg.db, keytype = "ENTREZID", ont = "BP", pAdjustMethod = "BH",  pvalueCutoff = 1, qvalueCutoff = 1))
        if(!inherits(gobp, "NULL")){
          gobp = gobp@result
          gobp$Sample = patient
          gobpALL = rbind(gobpALL, gobp)
          cat("gobp --> DONE \n")
        }else{cat("SKIP gobp\n")}}
      
    }
  }
  
  return(list(reactomeALL, keggALL, gobpALL))
  
}


##### FOR LOGISTIC REGRESSION 

# create PATIENTxPATHWAY matrix and fill it with adjusted p-values from the over-representation tests performed on each patient
# "databaseALL" must contain the concatenated output of over-representation tests performed on several individuals with reactomePA and clusterProfiler packages
# "databaseALL" must have a "ID" column with pathway ID and "Sample" column with patient name
# "geneTable" is a PATIENTxGENE matrix with sample names in its first column
# "sampleID" is a CASE/CTRL classification table for the samples and "criterion" is the type of classification to be used (which must be a column name in "sampleID" table)
createPatientPathwayMatrixFromHypergeometricTestOutput <- function(databaseALL, geneTable, sampleID, criterion){
  
  cat("pre-processing for PATIENTxPATHWAY matrix\n")
  
  t = unique(databaseALL[,c(1,2)])
  if(TRUE %in% duplicated(t$Description)){
    dup = t[which(duplicated(t$Description)),]
    databaseALL[databaseALL$Description == dup$Description, ] = dup$ID}
  databaseALL = unique(databaseALL)
  
  # select all unique pathways
  pathways = sort(unique(databaseALL$ID))
  
  # create PATIENTxPATHWAY matrix with:
  # sample names as first column
  # all zeroes in second column (called "case_ctrl")
  # NA in as many columns as pathways
  table = data.frame(Sample = geneTable$Sample, case_ctrl = rep(0, dim(geneTable)[1]), matrix(data=NA, nrow=dim(geneTable)[1], ncol=length(pathways)))
  
  # assign pathways as column names and samples as rownames
  colnames(table)[-(1:2)] = pathways
  rownames(table) = table$Sample
  
  # remove patient not included in the study
  table = table[sampleID$Sample,]
  
  # loop over PATIENTxPATHWAY matrix rows and fill "case_ctrl" column with the classification number of each sample
  # the info is gotten from "sampleID" table according to the "criterion" specified 
  cat("filling case_ctrl column \n")
  for(r in 1:dim(table)[1]){table$case_ctrl[r] = sampleID[sampleID$Sample == table$Sample[r], criterion]}
  
  # exclude patients with NA in "case_ctrl"
  table = table[table$case_ctrl != "na",]
  table$case_ctrl = as.numeric(table$case_ctrl)
  
  # loop over "databaseALL" rows, get adjusted p-value at that row and write it in corresponding cell of PATIENTxPATHWAY matrix
  for(r in 1:dim(databaseALL)[1]){
    if(r%%1000 == 0){cat("filling adjusted p-value - ", r, "/", dim(databaseALL)[1],"\n")}
    if(databaseALL$Sample[r] %in% rownames(table)){table[databaseALL$Sample[r], as.character(databaseALL$ID[r])] = databaseALL$p.adjust[r]}}
  
  # substitute all NAs with mean of values in a certain column
  # remove all columns with only NAs
  to_remove = NULL
  for(c in 2:dim(table)[2]){
    if(!FALSE %in% is.na(table[,c])){to_remove = append(to_remove, c)
    }else{
      if(c%%100 == 0){cat("change NAs with column's mean - ", c, "/", dim(table)[2],"\n")}
      table[is.na(table[,c]),c] = mean(table[,c], na.rm=T)}}
  if(length(to_remove)>0){table = table[,-to_remove]}
  
  return(table)
  
}


# change ctrls and cases to 0 and 1 ("case_ctrl" column)
# perform logistic regression for "case_ctrl" against all pathways (one by one) with following command:
# glm(case_ctrl ~ pathway, family = binomial(link = 'logit'), data=pathwayTable)
# values in "pathway" are the adjusted p-values output by over-representation hypergeometric test
# "pathwayTable" is a PATIENTxPATHWAY matrix with CASE/CTRL classification as second column and 
# adjusted p-values for the other columns
# "databaseALL" contains the concatenated results of over-representation test on all samples
logisticRegressionOnPatientPathwayMatrix <- function(pathwayTable, databaseALL, case, ctrl){
  
  # CHANGE CASE_CTRL (the OUTCOME) TO 0 AND 1, USE SAMPLE NAMES AS ROWNAMES, REMOVE FIRST COLUMN
  cat("pre-processing for LR\n")
  pathwayTable = pathwayTable[(pathwayTable$case_ctrl == case) | (pathwayTable$case_ctrl == ctrl),]
  pathwayTable[pathwayTable$case_ctrl == ctrl, "case_ctrl"] = 0
  pathwayTable[pathwayTable$case_ctrl == case, "case_ctrl"] = 1
  rownames(pathwayTable) = pathwayTable[,1]
  pathwayTable = pathwayTable[,-1]
  
  # perform Logistic Regression with "case_ctrl" as outcome and each pathway's adjusted p-values as predictor
  # save results in LRmatrix (matrix for results of Logistic Regression)
  
  
  CI=NULL #<--------------------------------------------------------------------------------------- for CI
  LRmatrix=NULL
  for(c in 2:dim(pathwayTable)[2]){
    if(c%%100 == 0){cat("fitting LR models - ", c, "/", dim(pathwayTable)[2], "\n")}
    df = pathwayTable[,c(1,c)]
    fit <- glm(case_ctrl ~ ., family = binomial(link = 'logit'), data=df)
    results = as.data.frame(coef(summary(fit)))
    t = try(confint(fit, level=0.999999))
    if(!inherits(t, "try-error")){CI = rbind(CI, as.data.frame(confint(fit, level=0.999999))[-1,])} #<-------------------------------------------- for CI
    LRmatrix = rbind(LRmatrix, results[-1,])}
  
  # adjust p-values according to BH method, sort pathwayTable according to q-values
  colnames(LRmatrix)[4] = "PVAL"
  LRmatrix$QVAL<-round(p.adjust(LRmatrix$PVAL, method="BH"), digits=20)
  # add columns for CI in LRmatrix
  LRmatrix = cbind(LRmatrix, one=rep(NA,dim(LRmatrix)[1]), two=rep(NA,dim(LRmatrix)[1]))
  colnames(LRmatrix)[(dim(LRmatrix)[2]-1):dim(LRmatrix)[2]] = c("CI_low", "CI_up")
  # copy CI values in LRmatrix, only for pathways common to both tables
  sel = which(rownames(LRmatrix) %in% rownames(CI))
  if(length(sel)>0){LRmatrix[sel, c("CI_low", "CI_up")] = CI[rownames(LRmatrix)[sel], c(1, 2)]}
  
  # add column and compute if CI comprise 1
  LRmatrix$one_in_CI_99.9999 = rep(NA, dim(LRmatrix)[1])
  for(r in 1:dim(LRmatrix)[1]){
    if(!TRUE %in% is.na(LRmatrix[r,c("CI_low", "CI_up")])){
      if(LRmatrix[r,"CI_low"]<1 & LRmatrix[r,"CI_up"]>1){LRmatrix[r,"one_in_CI_99.9999"] = "yes"
      }else{LRmatrix[r,"one_in_CI_99.9999"] = "no"}}}
  
  # sort table according to adjusted p-value
  LRmatrix<-LRmatrix[order(LRmatrix$QVAL),]
  
  # remove "`" from pathway IDs, which remains after logistic regression in REACTOME and GOBP
  if(TRUE %in% grepl("`",rownames(LRmatrix))){rownames(LRmatrix)<-unlist(lapply(strsplit(rownames(LRmatrix), "`"), "[[", 2))}
  
  # create pathwayTable with pathway IDs and pathway name "pathId2name" 
  pathId2name = unique(databaseALL[,c(1,2)])
  rownames(pathId2name) = pathId2name$ID
  
  # add columns to LRmatrix and fill with pathway IDs and assign pathway names
  # using the info in "pathId2name" pathwayTable
  cat("converting pathway IDs to pathway names\n")
  LRmatrix = cbind(pathID=rownames(LRmatrix), pathName=rep(".", dim(LRmatrix)[1]), LRmatrix, stringsAsFactors=F)
  for(r in 1:dim(LRmatrix)[1]){LRmatrix$pathName[r] = pathId2name[LRmatrix$pathID[r],2]}
  
  return(LRmatrix)
  
}


####### FOR FISHER BASED ON OVER-REPRESENTATION TEST RESULTS ---> TO TEST!!!!

# create PATIENTxPATHWAY matrix and fill it with zero or one depending on the significance of over-representation test for a certain pathway in a certain patient
# "databaseALL" must contain the concatenated output of over-representation tests performed on several individuals with reactomePA and clusterProfiler packages
# "databaseALL" must have a "ID" column with pathway ID and "Sample" column with patient name
# "geneTable" is a PATIENTxGENE matrix with sample names in its first column
# "sampleID" is a CASE/CTRL classification table for the samples and "criterion" is the type of classification to be used (which must be a column name in "sampleID" table)
createPatientPathwayMatrixFromHypergeometricTestOutputBinary <- function(databaseALL, geneTable, sampleID, criterion, cutoff){
  
  cat("pre-processing for PATIENTxPATHWAY matrix")
  
  t = unique(databaseALL[,c(1,2)])
  if(TRUE %in% duplicated(t$Description)){
    dup = t[which(duplicated(t$Description)),]
    databaseALL[databaseALL$Description == dup$Description, ] = dup$ID}
  databaseALL = unique(databaseALL)
  
  # select all unique pathways
  pathways = sort(unique(databaseALL$ID))
  
  # create PATIENTxPATHWAY matrix with:
  # sample names as first column
  # all zeroes in second column (called "case_ctrl") and other columns (for the pathways)
  table = data.frame(Sample = geneTable$Sample, case_ctrl = rep(0, dim(geneTable)[1]), matrix(data=0, nrow=dim(geneTable)[1], ncol=length(pathways)))
  
  # assign pathway names to columns
  colnames(table)[-(1:2)] = pathways
  
  # loop over PATIENTxPATHWAY matrix rows and fill "case_ctrl" column with the classification number of each sample
  # the info is taken from "sampleID" table according to the "criterion" specified 
  for(r in 1:dim(table)[1]){
    cat("filling case_ctrl column \n")
    table$case_ctrl[r] = sampleID[sampleID$Sample == table$Sample[r], criterion]}
  
  # exclude patients with NA in "case_ctrl"
  table = table[table$case_ctrl != "na",]
  table$case_ctrl = as.numeric(table$case_ctrl)
  
  # loop over "databaseALL" rows and get adjusted p-value at that row and fill corresponding cell of PATIENTxPATHWAY matrix with one,
  # if adjusted p-value is < "cutoff"(specified in arguments)
  for(r in 1:dim(databaseALL)[1]){
    if(r%%1000 == 0){cat("filling adjusted p-value - ", r, "/", dim(databaseALL)[1],"\n")}
    if(as.numeric(databaseALL$p.adjust[r]) < cutoff){table[table$Sample == databaseALL$Sample[r], as.character(databaseALL$ID[r])] = 1}}
  
  return(table)
  
}

# take a PATIENTxPATHWAY matrix from previous script and perform fisher's test on pathways between cases and ctrls
pathwayCaseCtrlFisher_TOtest <- function(pathwayTable, databaseALL, case, ctrl){
  
  cat("\n pre-processing for FISHER\n")
  
  # SET SAMPLE NAMES AS ROWNAMES
  rownames(pathwayTable) = pathwayTable[,1]
  
  # split table in Cases and Ctrls (by row), remove first and second columns (sample names and case/ctrl)
  pathCase = pathwayTable[pathwayTable$case_ctrl == case,]
  pathCase = pathCase[,-(1:2)]
  pathCtrl = pathwayTable[pathwayTable$case_ctrl == ctrl,]
  pathCtrl = pathCtrl[,-(1:2)]
  
  # create Fmatrix, with event counts for each gene in cases and ctrls
  caseCounts<-apply(pathCase[,1:dim(pathCase)[2]], 2, sum)
  ctrlCounts<-apply(pathCtrl[,1:dim(pathCtrl)[2]], 2, sum)
  Fmatrix<-data.frame(caseCounts, ctrlCounts)
  
  # create Rmatrix and execute two.sided fisher's test with FDR p-adjustment
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=5))
  colnames(Rmatrix)<-c("PATHWAY", "%Case", "%Ctrl", "PVAL", "QVAL")
  sumCase<-dim(pathCase)[1]
  sumCtrl<-dim(pathCtrl)[1]
  for(i in 1:dim(Fmatrix)[1]){
    if(i%%1000 == 0){cat("pathway #", i, "\n")}
    posCase<-Fmatrix[i, 'caseCounts']
    posCtrl<-Fmatrix[i, 'ctrlCounts']
    negCase<-sumCase-posCase
    negCtrl<-sumCtrl-posCtrl
    
    x<-matrix(c(posCase, posCtrl, negCase, negCtrl), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=20)
      
      Rmatrix[i, 2]<-posCase/sumCase*100
      Rmatrix[i, 3]<-posCtrl/sumCtrl*100
      Rmatrix[i, 4]<-pval
    } 
  }
  Rmatrix$PATHWAY<-rownames(Fmatrix)
  Rmatrix$QVAL<-round(p.adjust(Rmatrix$PVAL, method="BH"), digits=20)
  Rmatrix<-Rmatrix[order(Rmatrix$QVAL),]
  
  
  # create pathwayTable with pathway IDs and pathway name "pathId2name"
  pathId2name = unique(databaseALL[,c(1,2)])
  rownames(pathId2name) = pathId2name$ID
  
  # add columns to Rmatrix and fill with pathway IDs and assign pathway names
  # using the info in "pathId2name" pathwayTable
  cat("converting pathway IDs to pathway names\n")
  Rmatrix = cbind(pathID=Rmatrix[,1], pathName=rep(".", dim(Rmatrix)[1]), Rmatrix[,2:5], stringsAsFactors=F)
  for(r in 1:dim(Rmatrix)[1]){Rmatrix$pathName[r] = pathId2name[Rmatrix$pathID[r],2]}
  
  return(Rmatrix)
  
}



####### FOR FISHER BASED ON PATHWAYS ALTERED IN EACH PATIENT

# fill 1 in PATIENTxPATHWAY matrix whenever at least a gene in a pathway is altered in a patient
# take a PATIENTxPATHWAY matrix of zeroes and a pathwayID:entrezID table
fillPatientPathwayTable <- function(pathIds, patients, path2entr){
  
  # create a PATIENTxPATHWAY matrix of zeroes
  pathTable = data.frame(matrix(rep(0, length(patients)*length(pathIds)), ncol=length(pathIds)))
  colnames(pathTable) = pathIds
  rownames(pathTable) = patients
  
  # loop over pathways (columns)
  for(c in 1:dim(pathTable)[2]){
    
    if(c%%10==0){cat("filling data for pathway", c, "/", dim(pathTable)[2], "\n")}
    # get pathway[c] ID
    pathway = colnames(pathTable)[c]
    
    # get entrez ids contained in pathway[c]
    ids = path2entr[path2entr[,1] == pathway, "ENTREZID"]
    
    # convert entrez ids to symbols
    genes = as.vector(unlist(mget(ids,envir=org.Hs.egSYMBOL, ifnotfound=NA)))
    genes = genes[complete.cases(genes)]
    
    # remove genes that are not contained in PATIENTxGENE matrix
    genes = genes[which(genes %in% colnames(geneTable))]
    
    # subset PATIENTxGENE matrix for genes contained in pathway[c]
    if(length(genes)>0){
      df = as.data.frame(geneTable[,genes], row.names = patients)
      if(length(genes)==1){colnames(df) = genes} # deal with cases of only one gene
      
      # loop over patients (rows)
      for(r in 1:dim(pathTable)[1]){
        # get patient[r] name
        patient = rownames(pathTable)[r]
        # if patient[r] has at least one event in at least a gene in pathway[c] ----> fill 1 in PATIENTxPATHWAY matrix
        if(sum(df[patient,]) > 0){pathTable[patient,pathway] = 1}}
    }
  }
  
  # add patient names as first column
  pathTable = cbind(Sample=rownames(pathTable), pathTable)
  
  return(pathTable)
}  


# take a PATIENTxPATHWAY MATRIX zeroes and ones
# perform fisher's test on pathways between cases and ctrls
pathwayCaseCtrlFisher <- function(pathTable, selCase, selCtrl){
  
  # "selCase" and "selCtrl" must be data frames containing a column "Sample" with the sample names of CASES and CTRLS, respectively
  
  # remove first column (patient names)
  pathTable = pathTable[,-1]
  
  # split table in Cases and Ctrls (by row)
  cat("splitting matrix in cases/ctrls\n")
  pathCase = pathTable[which(rownames(pathTable) %in% selCase$Sample),]
  pathCtrl = pathTable[which(rownames(pathTable) %in% selCtrl$Sample),]
  
  # create Fmatrix, with event counts for each path in cases and ctrls
  cat("creating Fmatrix with event counts\n")
  caseCounts<-apply(pathCase[,1:dim(pathCase)[2]], 2, sum)
  ctrlCounts<-apply(pathCtrl[,1:dim(pathCtrl)[2]], 2, sum)
  Fmatrix<-data.frame(caseCounts, ctrlCounts)
  
  # create Rmatrix and execute two.sided fisher's test with FDR p-adjustment
  Rmatrix<-as.data.frame(matrix(rep(NA), dim(Fmatrix)[1], ncol=6))
  colnames(Rmatrix)<-c("PATHID", "PATHNAME", "%Case", "%Ctrl", "PVAL", "QVAL")
  sumCase<-dim(pathCase)[1]
  sumCtrl<-dim(pathCtrl)[1]
  for(i in 1:dim(Fmatrix)[1]){
    if(i%%100 == 0){cat("fisher on pathway", i, "/", dim(Fmatrix)[1], "\n")}
    posCase<-Fmatrix[i, 'caseCounts']
    posCtrl<-Fmatrix[i, 'ctrlCounts']
    negCase<-sumCase-posCase
    negCtrl<-sumCtrl-posCtrl
    
    x<-matrix(c(posCase, posCtrl, negCase, negCtrl), ncol=2, byrow=T)
    if(sum(!is.na(x))==4){
      
      pval<-round(fisher.test(x,alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)$p.value, digits=20)
      
      Rmatrix[i, "%Case"]<-posCase/sumCase*100
      Rmatrix[i, "%Ctrl"]<-posCtrl/sumCtrl*100
      Rmatrix[i, "PVAL"]<-pval
    } 
  }
  Rmatrix$PATHID<-rownames(Fmatrix)
  Rmatrix$QVAL<-round(p.adjust(Rmatrix$PVAL, method="BH"), digits=20)
  Rmatrix<-Rmatrix[order(Rmatrix$QVAL),]
  
  return(Rmatrix)
}



#------------------------------------------------------ FOR FINAL DATA ARRANGEMENT -------------------------------------------

# take Rmatrix/LRmatrix and a pathwayID:symbol table, find genes contained in each pathway and
# repeat each row as many times as genes are contained in each pathway
# column name for pathway id must be "pathID"
# column name for gene symbol must be "symbol"
getGeneSymbolsAndCreateNewRmatrix <- function(resTable, path2symbol){
  
  resTable2=NULL
  
  # loop over resTable rows (pathways)
  for(r in 1:dim(resTable)[1]){
    
    if(r%%10==0){cat("pathway", r, "/", dim(resTable)[1], "\n")}
    # get pathway[r]
    pathway = as.character(resTable$pathID[r])
    # get genes in pathway[r]
    genes = path2symbol[path2symbol$pathID==pathway,"symbol"]
    genes = genes[complete.cases(genes)]
    
    if(length(genes)>0){
      # repeat row [r] in resTable for N=length(genes) times and append everything to resTable2
      df=NULL
      for(l in 1:length(genes)){df = rbind(df, resTable[r,])}
      df$symbol = genes
      resTable2 = rbind(resTable2, df)
    }
    
  }
  
  return(resTable2)
  
}


# take the output of logistic regression on pathways and calculate frequency and percentage of events in cases and ctrls for each gene
# "LRmatrix2" is the modified output of logistic regression (LRmatrix); it has a row for each gene in each pathway
geneFrequencyAndPercentageInCasesAndCtrls <- function(LRmatrix2, geneTable, sampleID, criterion, case, ctrl){
  
  # IN CASE OF LOGISTIC REGRESSION RESULTS TABLE, remove columns 4 and 5 (which contain logistic regression coefficients)
  if("Estimate" %in% colnames(LRmatrix2)){LRmatrix2 = LRmatrix2[,-(4:5)]}
  
  # add gene frequency and gene percentage columns (for cases and ctrl)
  LRmatrix2$frequency_case = rep(0, dim(LRmatrix2)[1])
  LRmatrix2$frequency_ctrl = rep(0, dim(LRmatrix2)[1])
  LRmatrix2$percentage_case = rep(0, dim(LRmatrix2)[1])
  LRmatrix2$percentage_ctrl = rep(0, dim(LRmatrix2)[1])
  
  # split PATIENTxGENE matrix in cases and ctrl PATIENTxGENE matrices
  selCase<-sampleID[which(sampleID[,criterion]==case),"Sample"]
  caseTable = geneTable[which(geneTable$Sample %in% selCase),]
  selCtrl<-sampleID[which(sampleID[,criterion]==ctrl),"Sample"]
  ctrlTable = geneTable[which(geneTable$Sample %in% selCtrl),]
  
  # loop over LRmatrix rows (pathway:gene pair)
  for(r in 1:dim(LRmatrix2)[1]){
    if(r%%100==0){cat("count events for gene", r, "/", dim(LRmatrix2)[1], "\n")}
    gene = LRmatrix2$symbol[r]
    # if gene is altered in cases, calculate frequency and percentage of that event in cases and write in LRmatrix
    if(gene %in% colnames(caseTable)){
      caseCounts = sum(caseTable[1:dim(caseTable)[1], gene])
      LRmatrix2$frequency_case[r] = caseCounts
      LRmatrix2$percentage_case[r] = (caseCounts/dim(caseTable)[1])*100}
    # if gene is altered in ctrls, calculate frequency and percentage of that event in ctrls and write in LRmatrix
    if(gene %in% colnames(ctrlTable)){
      ctrlCounts = sum(ctrlTable[1:dim(ctrlTable)[1], gene])
      LRmatrix2$frequency_ctrl[r] = ctrlCounts
      LRmatrix2$percentage_ctrl[r] = (ctrlCounts/dim(ctrlTable)[1])*100}}
  
  # remove genes with no events
  LRmatrix2 = LRmatrix2[which(LRmatrix2$frequency_ctrl>0 | LRmatrix2$frequency_case>0),]
  
  return(LRmatrix2)
  
}


#### FOR ANNA (control if previous two are the same)

# for a list of genes, compute percentage of an event in cases and ctrls
# remove genes with no data
# for remaining genes, compute ratio and difference between percentage in cases and ctrls
# "resTable" must have a column with genes titled "symbol"
# "eventMatrix" is the PATIENTxGENE matrix for that event
computePercentageInCaseAndCtrlForEvent <- function(resTable, eventMatrix, selCase, selCtrl){
  
  # get gene list
  genes = resTable$symbol
  
  # subset "eventMatrix" for cases and ctrl
  mat_case = eventMatrix[selCase,]
  mat_ctrl = eventMatrix[selCtrl,]
  
  # get total number of cases and ctrls(n. of rows)
  totCase = dim(mat_case)[1]
  totCtrl = dim(mat_ctrl)[1]
  
  # if a gene has events (i.e. is contained in eventMatrix) compute the percentage of its events in cases and ctrls
  resTable$percentage_case = sapply(genes, function(x) if(x %in% colnames(mat_case)){sum(mat_case[, x])/totCase*100})
  resTable$percentage_ctrl = sapply(genes, function(x) if(x %in% colnames(mat_ctrl)){sum(mat_ctrl[, x])/totCtrl*100})
  
  # get index of genes that were not altered in cases (genes not appearing in eventMatrix)
  sel = which(sapply(resTable$percentage_case, function(x) is.null(x)))
  # remove selected genes
  if(length(sel)>0){resTable = resTable[-sel,]}
  
  # as new columns are saved as lists, convert to "numeric" (this would not be possible if NULLs were not removed in previous step)
  resTable$percentage_case = as.numeric(resTable$percentage_case)
  resTable$percentage_ctrl = as.numeric(resTable$percentage_ctrl)
  
  # compute ratio and difference between percentage in cases and ctrls
  resTable$ratio_case_ctrl = resTable$percentage_case/resTable$percentage_ctrl
  resTable$diff_case_ctrl = resTable$percentage_case - resTable$percentage_ctrl
  
  return(resTable)
  
}

