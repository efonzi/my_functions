
# coding: utf-8


def session_info():

    ### PRINT SESSION INFO TO STDOUT
    ### TO BE USED INSIDE PYTHON SCRIPTS    
    
    import IPython
    import subprocess as sp
    
    for i in IPython.sys_info().splitlines():
        print(i)
        
    for i in sp.run("pip3 freeze", shell=True, stdout=sp.PIPE).stdout.decode().splitlines():
        print(i)    





    
def overRepresentationFromGeneSets_02(genes, geneSets, symbols=None):
    '''Perform over-representation test on a set of 'genes', given a 'geneSets' background.
    See overRepresentationFromGeneSets_01 for details.
    
    CHANGES TO PREVIOUS VERSION
    - 'N' is counted only among the genes that belong to the gene set background
    - Name of column 'hits_count' in output was changed to 'hits'
    - Added column 'pathSize' in output with total number of genes contained in pathway (n)
    - Added column 'total_genes' in output with total number of altered genes
    - Added column 'tested_genes' in output with total number of tested genes (N)
    - Added column 'background' in output with total number of background genes (M)
    
    CREATED
    171004
    '''
    
    import numpy as np
    import pandas as pd
    from scipy.stats import hypergeom as hypg
    from statsmodels.sandbox.stats.multicomp import multipletests
    
    
    # in case 'symbols' argument was not provided, derive it from 'geneSets'
    if not symbols:
        symbols = list(set([s for p in geneSets for s in p['symbols']]))
    
    # get 'M'
    M = len(symbols)
    
    # get 'N'
    N = len([g for g in genes if g in symbols])
    
    # loop over pathways in 'geneSets' and check if they contain any of 'genes'
    # if so, save hit genes in a nested list (pathways with no hits will have an empty list)
    hits = [[g for g in genes if g in p['symbols']] for p in geneSets]
    
    # create another list with the counts of the hits per pathway
    hits_count = [len(h) for h in hits]
    
    # get number of genes contained in each pathway
    pathSize = [len(p['symbols']) for p in geneSets]
    
    # create DF with all pathways from 'geneSets' and the columns: 'hits', pathID', 'pathName', 'x'
    df = pd.DataFrame({'pathID': [p['pathID'] for p in geneSets], 
             'pathName': [p['pathName'] for p in geneSets], 
             'hits': hits, 
             'hits_count': hits_count,
             'pathSize': pathSize,
             'total_genes': [str(len(genes))] * len(geneSets),
             'tested_genes': [str(N)] * len(geneSets),
             'background': [str(M)] * len(geneSets),
             'pval': ['.'] * len(geneSets)})
    
    # subset DF for pathways that have at least one hit (x>0)
    df = df[[h > 0 for h in df.hits_count]]
    # reset DF indeces
    df = df.reset_index()
    # remove 'index' column, containing the old indeces
    df = df.drop('index', 1)
    
    # loop over pathways and perform OVER-REPRESENTATION test
    for i in range(len(df)):
        
        # get current pathID
        pat = df.loc[i, 'pathID']
        
        # get 'n' for current pathID
#        n = len([g for p in geneSets if p['pathID'] == pat for g in p['symbols']])
        n = len([p['symbols'] for p in geneSets if p['pathID'] == pat][0])
        
        # compute hypergeometric pval and save to 'pval' column
        df.loc[i, 'pval'] = hypg.sf(df.loc[i, 'hits_count']-1, M, n, N)

    # adjust pvals (Benjamini-Hochberg)
    df['fdr'] = multipletests(df['pval'], alpha=0.05, method='fdr_bh')[1]

    # sort by adjusted pvals
    df = df.sort_values('pval')
    
    return(df)



def getCytobandFromNCBI_01(geneSymbol):
    '''Using Bio.Entrez library in Biopython, send a request to NCBI website and extract the info of the cytoband were a gene symbol is located.
    
    It works in 2 steps:
      1. a query with the gene symbol to get the corresponding Entrez ID (if there is more than one ID, the first will be chosen
      2. a query with the Entrez ID to get the cytoband
    
    INPUT
    geneSymbol: gene symbol (string)
    
    OUTPUT
    the cytoband where the gene is located (if no info is available, an empty string is returned)
    
    CAVEAT
    The following code must be run before calling this function
    >>> from Bio import Entrez
    >>> Entrez.email = 'a_real_mail_address'
    
    CREATED
    171001'''
    
    import re
    from Bio import Entrez
    Entrez.email = 'eugenio.fonzi2@unibo.it'
    
    # build query string
    query = "(" + geneSymbol + "[gene]) AND (Homo sapiens[orgn])"
    
    # send request to Entrez
    handle = Entrez.esearch(db="gene", term=query, tettype="gene", retmode="xml")
    
    # de-serialize output XML
    record = Entrez.read(handle)
    
    # close handle
    handle.close
    
    # extract entrez ID of gene
    if len(record['IdList']) > 0:
        
        ID = record["IdList"][0]
        
        # send request to Entrez
        handle = Entrez.efetch(db="gene" ,id=ID, tettype="gene", retmode="xml")
        
        # de-serialize output XML
        record2 = Entrez.read(handle)
        
        # close handle
        handle.close
    
        # extract cytoband info
        if 'Entrezgene_location' in record2[0].keys() and 'Maps_display-str' in record2[0]['Entrezgene_location'][0].keys():
            cytoband = record2[0]['Entrezgene_location'][0]['Maps_display-str']
    
        else:
            cytoband = ''

    else:
        cytoband = ''
        
    
    # get chromosome and chromosome arm
    if cytoband != '':
        
        # extract chromosome arm info from cytoband
        chr_arm = re.findall('^[^pq]+[pqcen]', cytoband)[0]
        
        # extract chromosome info from chromosome arm        
        chromosome = re.findall('^[^pqc]+', chr_arm)[0]
        
    else:
        chr_arm = ''
        chromosome = ''

    return (cytoband, chr_arm, chromosome)




def overRepresentationFromGeneSets_01(genes, geneSets, symbols=None):
    '''Perform over-representation test on a set of 'genes', given a 'geneSets' background.
    The test is based on the hypergeometric distribution and requires 4 parameters:
     - 'M' is the genetic background for the test (i.e. the total number of genes among all pathways)
     - 'N' is the number of sampled genes (i.e. the genes found altered in a sample) --> condition A
     - 'n' is the number of genes present in a pathway --> condition B
     - 'x' is the number of altered genes present in a pathway --> condition A & B
    
    INPUT
    genes: list of gene symbols on which we want to perform the test (i.e. all altered genes in a sample)
    geneSets: list of dictionaries. Each dictionary contains info for a different pathway and must have at least 3 keys:
        - 'pathID' --> pathway ID (string)
        - 'pathName' --> pathway name (string)
        - 'symbols' --> gene symbols contained in pathway (list of strings)

    OUTPUT
    a pandas dataframe with a row for each pathway and the columns: 
        - 'hits' --> symbol of altered genes that belong to the pathway (list of strings)
        - 'pathID' --> pathway ID (string)
        - 'pathName' --> pathway name (string)
        - 'x' --> number of altered genes that belong to the pathway (int)
        - 'pval' --> p-value of over-representation test (float)
        - 'fdr' --> adjusted p-value (Benjamini-Hochberg)

    The output dataframe is sorted by the column 'pval' in increasing order.
    
    CREATED
    170905
    
    MODIFIED
    170920
    '''
    
    import numpy as np
    import pandas as pd
    from scipy.stats import hypergeom as hypg
    from statsmodels.sandbox.stats.multicomp import multipletests
    
    
    # in case 'symbols' argument was not provided, derive it from 'geneSets'
    if not symbols:
        symbols = list(set([s for p in geneSets for s in p['symbols']]))
    
    
    # loop over pathways in 'geneSets' and check if they contain any of 'genes'
    # if so, save hit genes in a nested list (pathways with no hits will have an empty list)
    hits = [[g for g in genes if g in p['symbols']] for p in geneSets]
    
    # create another list with the counts of the hits per pathway
    x = [len(h) for h in hits]
    
    
    # create DF with all pathways from 'geneSets' and the columns: 'hits', pathID', 'pathName', 'x'
    df = pd.DataFrame({'pathID': [p['pathID'] for p in geneSets], 
             'pathName': [p['pathName'] for p in geneSets], 
             'hits': hits, 
             'x': x})
    
    # subset DF for pathways that have at least one hit (x>0)
    df = df[[h > 0 for h in df.x]]
    # reset DF indeces
    df = df.reset_index()
    # remove 'index' column, containing the old indeces
    df = df.drop('index', 1)
    
    # add empty column were 'pval' will be stored
    df['pval'] = ['.'] * len(df)
    
    # get 'M'
    M = len(symbols)
    
    # get 'N'
    N = len(genes)
    
    # loop over pathways and perform OVER-REPRESENTATION test
    for i in range(0, len(df)):
        
        # get current pathID
        pat = df.loc[i, 'pathID']
        
        # get 'n' for current pathID
        n = len([g for p in geneSets if p['pathID'] == pat for g in p['symbols']])
        
        # compute hypergeometric pval and save to 'pval' column
        df.loc[i, 'pval'] = hypg.sf(df.loc[i, 'x']-1, M, n, N)

    # adjust pvals (Benjamini-Hochberg)
    df['fdr'] = multipletests(df['pval'], alpha=0.05, method='fdr_bh')[1]

    # sort by adjusted pvals
    df = df.sort_values('pval')
    
    return(df)




def reactomeWebServiceEnrichment_01(genes, sampleName):
    '''Submit a request to REACTOME servers to perform an OVER-REPRESENTATION analysis, handle and return the output.
    
    INPUT
    genes: list of gene symbols.
    sampleName: string.
    
    OUTPUT
    a pandas dataframe with the columns 'fdr', 'pathID', 'pathName', 'pval' 
    
    CREATED
    170905.'''
    
    # dependencies
    import subprocess as sp
    import json
    import pandas as pd
    
    # define all pieces of the command string
    genes_arg = '-d "$(printf \'#' + sampleName + '\n' + '\n'.join(genes) + '\')"'
    header = '-H "Content-Type: text/plain"'
    post_url = 'X POST --url http://www.reactome.org/AnalysisService/identifiers/projection/'
    
    # paste them together
    cmd = ' '.join(['curl', header, genes_arg, post_url])
    
    # run command through SHELL and capture its STDOUT (a JSON)
    data = sp.Popen(cmd, shell=True, stdout=sp.PIPE).stdout.read().decode()
    
    # translate JSON to python data structures
    data = json.loads(data)
    
    # extract desired result info as a dict and wrap it in a pandas DF 
    data = pd.DataFrame({'fdr' : [p['entities']['fdr'] for p in data['pathways']], 
          'pval' : [p['entities']['pValue'] for p in data['pathways']], 
          'pathID' : [p['stId'] for p in data['pathways']], 
          'pathName' : [p['name'] for p in data['pathways']]})

    return(data)



def query_wes_anna_results_01(position1, position2):
    '''created on 170822 to query the variants called by Marco's pipeline on ANNA's WES data.
    Given a range of genomic coordinates, retrieve all variant calls lying in that range.
    Files from all filtering steps (from the raw calls to the filtered somatic variants) and both tools (Mutect and Varscan) 
    are queried.
    
    INPUT
    position1 = START chromosomic position (as integer)
    position2 = END chromosomic position (as integer)
    
    OUTPUT
    a list of pandas dataframes, each containing the hits from a different file
    
    QUERIED FILES
    1a) Mutect raw output
    1b) Varscan SNP raw output
    1c) Varscan INDEL raw output
    2a) Mutect ANNOSET*
    2b) Varscan SNP ANNOSET
    2c) Varscan INDEL ANNOSET
    3a) Mutect annotated variants
    3b) Varscan SNP annotated variants
    3c) Varscan INDEL annotated variants
    4a) all somatic variants (with polymorphisms)
    4b) all somatic variants (without polymorphisms)
    
    WARNINGS
    the path is specific for ANNA WES results
    the batch "nextera_merged" is excluded from the analysis'''

    import pandas as pd
    import os
    import fnmatch as fm
    
    # save pre and sub path
    pre_path = '/home/PERSONALE/eugenio.fonzi2/wes_anna/batches/'
    sub_path = '/Analysis/Results/'
    
    # get list of batch names
    batches = os.listdir(pre_path)

    # remove batch 'nextera_merged' because it has not yet been run
    batches = [b for b in batches if b!='nextera_merged']
    
    
    ### Create a series with batches names as indices and lists of files contained in each "Results" sub_folder as values
    
    # create an empty series with batches names as indices
    results = pd.Series(['NaN']*len(batches), index=batches)
    
    # loop over batches
    for batch in batches:
        
        # save as a list all file names contained in "Results" folder of "batch"
        results[batch] = os.listdir(pre_path + batch + sub_path)
    
    
    ### Create empty DFs to store query results     
    
    # create empty DF for 1a) MUTECT
    mutect = pd.DataFrame()
    # create empty DF for 1b) VARSCAN SNPS
    varsc_snp = pd.DataFrame()
    # create empty DF for 1c) VARSCAN INDELS
    varsc_ind = pd.DataFrame()
    
    # create empty DF for 2a) ANNOSET MUTECT
    aset_mut = pd.DataFrame()
    # create empty DF for 2b) ANNOSET VARSCAN SNPS
    aset_varsc = pd.DataFrame()
    # create empty DF for 2c) ANNOSET VARSCAN INDELS
    aset_varsc_ind = pd.DataFrame()
    
    # create empty DF for 3a) SOMATIC ANNOTATED MUTECT
    annotated_mut = pd.DataFrame()
    # create empty DF for 3b) SOMATIC ANNOTATED VARSCAN SNPS
    annotated_varsc = pd.DataFrame()
    # create empty DF for 3c) SOMATIC ANNOTATED VARSCAN INDELS
    annotated_varsc_ind = pd.DataFrame()
    
    # create empty DF for 4a) SOMATIC VARIANTS
    som_var = pd.DataFrame()
    # create empty DF for 4b) FILTERED SOMATIC VARIANTS
    filt_som_var = pd.DataFrame()
    

    ### Query all file types in all batches and concatenate rows that match position
    
    # loop over batches
    for batch in batches:
    
        # select and load files of 1a) MUTECT for batch
        files = [file for file in results[batch] if fm.fnmatch(file, 'call_stats*.out')]
        
        # loop over selected file names
        for file in files:
            
            # load file (first row is empty, so I specify "header=1")
            df = pd.read_table(pre_path + batch + sub_path + file, sep='\t', low_memory=False, header=1)
            
            # append to DF the variants (rows) that contain desired chromosomal position
            mutect = mutect.append(df[[position1 <= x <= position2 for x in list(df['position'])]])
    
    
    
        # select and load files of 1b) VARSCAN SNPS for batch
        files = [file for file in results[batch] if fm.fnmatch(file, '*snp.pairs*')]
        
        # loop over selected file names
        for file in files:
            
            # load file
            df = pd.read_table(pre_path + batch + sub_path + file, sep='\t', low_memory=False)
            
            # add file name at the end of each row
            df['file_name'] = [file]*len(df)
            
            # append to DF the variants (rows) that contain desired chromosomal position
            varsc_snp = varsc_snp.append(df[[position1 <= x <= position2 for x in list(df['position'])]])
    
     
    
        # select and load files of 1c) VARSCAN INDELS for batch
        files = [file for file in results[batch] if fm.fnmatch(file, '*indels.pairs*')]
        
        # loop over selected file names
        for file in files:
            
            # load file
            df = pd.read_table(pre_path + batch + sub_path + file, sep='\t', low_memory=False)
            
            # add file name at the end of each row
            df['file_name'] = [file]*len(df)
            
            # append to DF the variants (rows) that contain desired chromosomal position
            varsc_ind = varsc_ind.append(df[[position1 <= x <= position2 for x in list(df['position'])]])
    
    
        # (ANNOSET TABLES DON'T HAVE HEADER, SO I ADD IT WITH "names=" ARGUMENT)
    
        # load 2a) ANNOSET MUTECT table for batch
        df = pd.read_table(pre_path + batch + sub_path + 'annoSet.tsv', sep='\t', low_memory=False, names=range(0,37))
        # subset for chromosomal position
        aset_mut = aset_mut.append(df[[position1 <= x <= position2 for x in list(df[1])]])
    
        # load 2b) ANNOSET VARSCAN SNPS table for batch
        df = pd.read_table(pre_path + batch + sub_path + 'annoSetVarscan.tsv', names=range(0,25), sep='\t', low_memory=False)
        # subset for chromosomal position
        aset_varsc = aset_varsc.append(df[[position1 <= x <= position2 for x in list(df[1])]])
        
        # load 2c) ANNOSET VARSCAN INDELS table for batch
        df = pd.read_table(pre_path + batch + sub_path + 'annoSetVarscan_indels.tsv', sep='\t', low_memory=False, names=range(0,25))
        # subset for chromosomal position
        aset_varsc_ind = aset_varsc_ind.append(df[[position1 <= x <= position2 for x in list(df[1])]])
        
    
        # load 3a) SOMATIC ANNOTATED MUTECT table for batch
        df = pd.read_table(pre_path + batch + sub_path + 'somatic_mutation_hg19_annotated.tsv', sep='\t', low_memory=False)
        # subset for chromosomal position
        annotated_mut = annotated_mut.append(df[[position1 <= x <= position2 for x in list(df['Start'])]])
        
        # load 3b) SOMATIC ANNOTATED VARSCAN SNPS table for batch
        df = pd.read_table(pre_path + batch + sub_path + 'somatic_mutation_varscan_hg19_annotated.tsv', sep='\t', low_memory=False)
        # subset for chromosomal position
        annotated_varsc = annotated_varsc.append(df[[position1 <= x <= position2 for x in list(df['Start'])]])
       
        # load 3c) SOMATIC ANNOTATED VARSCAN INDELS table for batch
        df = pd.read_table(pre_path + batch + sub_path + 'somatic_indels_varscan_hg19_annotated.tsv', sep='\t', low_memory=False)
        # subset for chromosomal position
        annotated_varsc_ind = annotated_varsc_ind.append(df[[position1 <= x <= position2 for x in list(df['Start'])]])
    
        
        # load 4a) SOMATIC VARIANTS table for batch
        df = pd.read_table(pre_path + batch + sub_path + 'somatic_variants.tsv', sep='\t', low_memory=False)
        # subset for chromosomal position
        som_var = som_var.append(df[[position1 <= x <= position2 for x in list(df['START'])]])
        
        # load 4b) FILTERED SOMATIC VARIANTS table for batch
        df = pd.read_table(pre_path + batch + sub_path + 'somatic_variants_filtered.tsv', sep='\t', low_memory=False)
        # subset for chromosomal position
        filt_som_var = filt_som_var.append(df[[position1 <= x <= position2 for x in list(df['START'])]])
    
    return [mutect, varsc_snp, varsc_ind, aset_mut, aset_varsc, aset_varsc_ind, annotated_mut, annotated_varsc, annotated_varsc_ind, som_var, filt_som_var]
    



