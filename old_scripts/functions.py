'''FIRST VERSION OF CUSTOM FUNCTIONS FILE. ONLY USED IN '170823_queryWesAnna.py'.
SUBSEQUENT FUNCTIONS ARE IN 'myFunctions.py'.'''


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