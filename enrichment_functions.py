
def fisherTestOnGenes(ctrl, case, sort_by='pval'):
    
    '''
    For N gene symbols, perform Fisher's exact test N times. For each gene, the test compares the presence/absence of a certain event in that gene between two groups of samples.
    
    INPUT:
    'ctrl' and 'case' are two matrices with sample names on the rows and gene symbols on the columns. The cells are filled with 1 (if the gene has the event in that patient) or 0 (in case of no events). The gene symbols (i.e. the colum names) must be the same between 'ctrl' and 'case'. 
    '''
    
    from statsmodels.sandbox.stats.multicomp import multipletests
    import scipy.stats as ss
    import pandas as pd
    import numpy as np
    
    # set list of column names for result tables of FISHER
    columns = ['gene', 'ctrl_prop', 'case_prop', 'odds_ratio', 'pval']

    # create empty DF to store results
    DF = pd.DataFrame(columns=columns)
    
    # get list of gene symbols
    symbols = list(ctrl.columns)
    
    # loop over symbols
    for symb in symbols:
        
        # compute fields of 2x2 contingency table
        # t=top, b=bottom, l=left, r=right
        tl = sum(ctrl[symb])
        bl = len(ctrl[symb]) - tl
        tr = sum(case[symb])
        br = len(case[symb]) - tr

        # perform test
        fish = ss.fisher_exact([[tl, tr], [bl, br]])
        
        # create one-row DF with test results
        df = pd.DataFrame([[symb, tl/(tl+bl), tr/(tr+br), fish[0], fish[1]]], index=[symb], columns=columns)
        
        # append DF as new row of result table
        DF = DF.append(df)

        # apply FDR < 0.05
        DF['fdr'] = multipletests(DF['pval'], alpha=0.05, method='fdr_bh')[1]
        
        # compute log10 of FDR
        DF['log10_fdr'] = np.log10(DF['fdr'])
        
        # sort by un-adjusted p-values
        DF = DF.sort_values('pval')
        
        # sort by FDR in case specified in command options
        if sort_by == 'fdr':
            DF = DF.sort_values('fdr')
        
    return(DF)




def overRepresentationFromGeneSets_03(genes, geneSets, remove_no_hits=True):
    '''Perform over-representation test on a set of 'genes', given a 'geneSets' background.
    
    The test is based on the hypergeometric distribution and requires 4 parameters:
     - 'M' is the genetic background for the test (i.e. the total number of genes among all pathways)
     - 'N' is the number of sampled genes (i.e. the genes found altered in a sample) --> condition A
     - 'n' is the number of genes present in a pathway --> condition B
     - 'x' is the number of altered genes present in a pathway --> condition A & B
    
    INPUT
    genes: list of gene symbols on which we want to perform the test (i.e. all altered genes in a sample)
    geneSets: DATAFRAME with 3 columns:
        - 'pathID' --> pathway ID (string)
        - 'pathName' --> pathway name (string)
        - 'symbol' --> gene symbols contained in pathway, one row per symbol (string)
    
    OPTIONAL ARGUMENTS
    remove_no_hits: remove pathways with no hits among 'genes' (boolean, default True)

    OUTPUT
    a pandas dataframe with a row for each pathway and the columns: 
        - 'hits' --> symbols of altered genes that belong to the pathway (list of strings)
        - 'pathID' --> pathway ID (string)
        - 'pathName' --> pathway name (string)
        - 'hits_count' --> number of altered genes that belong to the pathway (int)
        - 'pval' --> p-value of over-representation test (float)   
        - 'pathSize' --> total number of genes contained in pathway (n)
        - 'total_genes' --> total number of altered genes
        - 'tested_genes' --> total number of tested genes (N)
        - 'background' --> total number of background genes (M)

    The output dataframe is sorted by the column 'pval' in increasing order.

    CHANGES TO PREVIOUS VERSION
    - changed structure of 'geneSets' input file, from list of dictionaries to DF. 
    - removed optional argument 'symbols'
    - added optional argument 'no_hits'
    - removed FDR computation
    
    CREATED
    171025
    '''

    import numpy as np
    import pandas as pd
    from scipy.stats import hypergeom as hypg
    from statsmodels.sandbox.stats.multicomp import multipletests
    
    
    # derive list of all 'symbols' from 'geneSets'
    symbols = list(set(geneSets.symbol))
    
    # get DF of unique pathID and pathName
    pathways = geneSets[['pathID', 'pathName']].drop_duplicates().copy()
    # reset DF indeces
    pathways = pathways.reset_index()

    # get 'M'
    M = len(symbols)
    
    # get 'N'
    N = len([g for g in genes if g in symbols])
    
    # create empty lists
    hits = []
    pathSize = []

    # loop over pathways
    for p in pathways.pathID:

        # get all symbols in pathway
        sym = geneSets[geneSets.pathID == p]['symbol'].tolist()

        # append number of 'sym' to 'pathSize'
        pathSize.append(len(sym))
        
        # get list of 'genes' that are contained in 'geneSets' and append it to 'hits'
        hits.append([g for g in genes if g in sym])
       
    # create another list with the counts of the hits per pathway
    hits_count = [len(h) for h in hits]
    
    # create DF with all pathways from 'geneSets' and following columns
    df = pd.DataFrame({'pathID': pathways.pathID, 
             'pathName': pathways.pathName, 
             'hits': hits, 
             'hits_count': hits_count,
             'pathSize': pathSize,
             'total_genes': [str(len(genes))] * len(pathways),
             'tested_genes': [str(N)] * len(pathways),
             'background': [str(M)] * len(pathways),
             'pval': [float('nan')] * len(pathways)})
    
    # separately store pathways with no hits 
    no_hits = df[df.hits_count == 0].copy()
    
    # subset DF for pathways that have at least one hit (x>0)
    df = df[df.hits_count > 0]
    # reset DF indeces
    df = df.reset_index()
    # remove 'index' column, containing the old indeces
    df = df.drop('index', 1)

    # loop over pathways and perform OVER-REPRESENTATION test
    for i in range(len(df)):
        
        # get 'n' for current pathID
        n = df.loc[i, 'pathSize']
        
        # compute hypergeometric pval and save to 'pval' column
        df.loc[i, 'pval'] = hypg.sf(df.loc[i, 'hits_count']-1, M, n, N)

    
    if remove_no_hits == False:
        # re-attach pathways with no hits
        df = pd.concat([df, no_hits])
        
    # sort by adjusted pvals
    df = df.sort_values('pval')
    
    # reset indeces
    df = df.reset_index()
    
    return(df)
