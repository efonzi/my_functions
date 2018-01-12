
##### OVER-REPRESENTATION TEST

hypergeometricTestOnMultipleIndividuals_01 <- function(sampleGene){

  # loop over samples and perform over-representation test for REACTOME, KEGG, GOBP pathways (based on hypergeometric distribution)
  # requires a SAMPLExGENE matrix of zeroes and ones "sampleGene"

  reactomeALL=NULL
  keggALL=NULL
  gobpALL=NULL
  
  # loop over table rows (samples)
  for(r in 1:dim(sampleGene)[1]){
    
    # get name of the sample
    sample = sampleGene[r,1]
    
    # select altered genes for the sample
    genes = colnames(sampleGene)[which(sampleGene[r,] == 1)]
    
    cat("     \n")
    cat("SAMPLE", r, "of", dim(sampleGene)[1], " - ", length(genes), "events", "\n")
    
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
            reactome$Sample = sample
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
          kegg$Sample = sample
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
          gobp$Sample = sample
          gobpALL = rbind(gobpALL, gobp)
          cat("gobp --> DONE \n")
        }else{cat("SKIP gobp\n")}}
      
    }
  }
  
  return(list(reactomeALL, keggALL, gobpALL))
  
}



##### FOR LOGISTIC REGRESSION 

createSamplePathwayMatrixFromHypergeometricTestOutput <- function(over_representation, sampleGene, classification, criterion){
  
  # create SAMPLExPATHWAY matrix and fill it with adjusted p-values from the over-representation tests performed on each sample
  # "over_representation" must contain the concatenated output of over-representation tests performed on several individuals with reactomePA and clusterProfiler packages
  # "over_representation" must have a "ID" column with pathway ID and "Sample" column with sample name
  # "sampleGene" is a SAMPLExGENE matrix with sample names in its first column
  # "classification" is a CASE/CTRL classification table for the samples and "criterion" is the type of classification to be used (which must be a column name in "classification" table)
  
  cat("pre-processing for SAMPLExPATHWAY matrix\n")
 
  #############################################################
  ##### SOMETIMES, THE SAME PATHNAME CORRESPOND TO MULTIPLE PATHIDS...
   
  # get unique DF of pathID/pathName
  t = unique(over_representation[,c('ID', 'Description')])
  
  # if any pathNames are duplicated (different pathID with more than one pathName)
  if(TRUE %in% duplicated(t$Description)){
    
    # get duplicated pathID/pathName
    dup = t[which(duplicated(t$Description)),]
    
    # get indeces of duplicated pathIDs in OVER-REPRESENTATION table
    sel = which(over_representation$ID %in% dup$ID)
    
    # remove duplicated
    over_representation = over_representation[-sel,]
  }  

  
  #############################################################
    
  # select all unique pathways
  pathways = sort(unique(over_representation$ID))
  
  # create SAMPLExPATHWAY matrix with:
  # sample names as first column
  # all zeroes in second column (called "case_ctrl")
  # NA in as many columns as pathways
  table = data.frame(Sample = sampleGene$Sample, case_ctrl = rep(0, dim(sampleGene)[1]), matrix(data=NA, nrow=dim(sampleGene)[1], ncol=length(pathways)))
  
  # assign pathways as column names and samples as rownames
  colnames(table)[-(1:2)] = pathways
  rownames(table) = table$Sample
  
  # remove sample not included in the study
  table = table[classification$Sample,]
  
  # loop over SAMPLExPATHWAY matrix rows and fill "case_ctrl" column with the classification number of each sample
  # the info is gotten from "classification" table according to the "criterion" specified 
  cat("filling case_ctrl column \n")
  for(r in 1:dim(table)[1]){table$case_ctrl[r] = classification[classification$Sample == table$Sample[r], criterion]}
  
  # exclude samples with NA in "case_ctrl"
  table = table[table$case_ctrl != "na",]
#  table$case_ctrl = as.numeric(table$case_ctrl)
  
  # loop over "over_representation" rows, get adjusted p-value at that row and write it in corresponding cell of SAMPLExPATHWAY matrix
  for(r in 1:dim(over_representation)[1]){
    if(r%%1000 == 0){cat("filling adjusted p-value - ", r, "/", dim(over_representation)[1],"\n")}
    if(over_representation$Sample[r] %in% rownames(table)){table[over_representation$Sample[r], as.character(over_representation$ID[r])] = over_representation$p.adjust[r]}}
  
  # substitute all NAs with mean of values in a certain column
  # remove all columns with only NAs
  to_remove = NULL
  for(c in 3:dim(table)[2]){
    if(!FALSE %in% is.na(table[,c])){to_remove = append(to_remove, c)
    }else{
      if(c%%100 == 0){cat("change NAs with column's mean - ", c, "/", dim(table)[2],"\n")}
      table[is.na(table[,c]),c] = mean(table[,c], na.rm=T)}}
  if(length(to_remove)>0){table = table[,-to_remove]}
  
  return(table)
  
}



logisticRegressionOnSamplePathwayMatrix <- function(samplePathway, over_representation, case, ctrl){

  # change ctrls and cases to 0 and 1 ("case_ctrl" column)
  # perform logistic regression for "case_ctrl" against all pathways (one by one) with following command:
  # glm(case_ctrl ~ pathway, family = binomial(link = 'logit'), data=samplePathway)
  # values in "pathway" are the adjusted p-values output by over-representation hypergeometric test
  # "samplePathway" is a SAMPLExPATHWAY matrix with CASE/CTRL classification as second column and 
  # adjusted p-values for the other columns
  # "over_representation" contains the concatenated results of over-representation test on all samples

  print(samplePathway$case_ctrl)
    
  # CHANGE CASE_CTRL (the OUTCOME) TO 0 AND 1, USE SAMPLE NAMES AS ROWNAMES, REMOVE FIRST COLUMN
  cat("pre-processing for LR\n")
  samplePathway = samplePathway[(samplePathway$case_ctrl == case) | (samplePathway$case_ctrl == ctrl),]
#  samplePathway[samplePathway$case_ctrl == ctrl, "case_ctrl"] = 0
#  samplePathway[samplePathway$case_ctrl == case, "case_ctrl"] = 1
  samplePathway$case_ctrl = sapply(samplePathway$case_ctrl, function(x) if(x == case){x = 1}else{x = 0})
  rownames(samplePathway) = samplePathway[,1]
  samplePathway = samplePathway[,-1]
  
  # perform Logistic Regression with "case_ctrl" as outcome and each pathway's adjusted p-values as predictor
  # save results in LRmatrix (matrix for results of Logistic Regression)
  
  print(samplePathway$case_ctrl)
  
  CI=NULL #<--------------------------------------------------------------------------------------- for CI
  LRmatrix=NULL
  for(c in 2:dim(samplePathway)[2]){
    if(c%%100 == 0){cat("fitting LR models - ", c, "/", dim(samplePathway)[2], "\n")}
    df = samplePathway[,c(1,c)]
    fit <- glm(case_ctrl ~ ., family = binomial(link = 'logit'), data=df)
    results = as.data.frame(coef(summary(fit)))
    t = try(confint(fit, level=0.999999))
    if(!inherits(t, "try-error")){CI = rbind(CI, as.data.frame(confint(fit, level=0.999999))[-1,])} #<-------------------------------------------- for CI
    LRmatrix = rbind(LRmatrix, results[-1,])}
  
  # adjust p-values according to BH method, sort samplePathway according to q-values
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
  
  # create samplePathway with pathway IDs and pathway name "pathId2name" 
  pathId2name = unique(over_representation[,c(1,2)])
  rownames(pathId2name) = pathId2name$ID
  
  # add columns to LRmatrix and fill with pathway IDs and assign pathway names
  # using the info in "pathId2name" samplePathway
  cat("converting pathway IDs to pathway names\n")
  LRmatrix = cbind(pathID=rownames(LRmatrix), pathName=rep(".", dim(LRmatrix)[1]), LRmatrix, stringsAsFactors=F)
  for(r in 1:dim(LRmatrix)[1]){LRmatrix$pathName[r] = pathId2name[LRmatrix$pathID[r],2]}
  
  return(LRmatrix)
  
}


