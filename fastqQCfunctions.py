
def parseFastqcZip(ziplist):
    
    '''
    Take a list of ZIP files that are FASTQC output (as absolute paths) and return....
    '''
    
    import pandas as pd
    import re
    from zipfile import ZipFile
    
    # loop over zipfiles
    for z in ziplist:

        # extract fastq name
        name = z.split('/')[-1].replace('_fastqc.zip', '')
        
        # open ZIP archive and open file with FASTQC result data
        # split data between all the different analyses of FASTQC
        # (each part ends with the string '>>END_MODULE')
        with ZipFile(z) as myzip:
            with myzip.open(name + '_fastqc/fastqc_data.txt') as myfile:
                fastqc = myfile.read().decode().split('>>END_MODULE')





def runFastqc(fastqlist, fastqcpath):
    
    '''
    Run FASTQC on an input list of fastq files (compressed or not) and save results to current directory

    'fastqlist' must comprise the absolute paths to the files
    'fastqcpath' must comprise the command name at the end
    '''
    
    import subprocess as sp
    
    l = len(fastqlist)
    
    # loop over files
    for i in range(l):
        
        if i%10 == 0:
            print('fastqc on file %d/%d' %(i, l))

        # get fastq path
        f = fastqlist[i]
        
        # get file name
        name = f.split('/')[-1]
        
        # build and run command (stderr is not needed)
        
        cmd = ' '.join([fastqcpath,
                        '-o .',
                        f])
        
        try:
            _ = sp.run(cmd, shell=True, stderr=sp.PIPE).stderr
            
        except:
            with open('failed_fastqc.log', 'a') as outfile:
                outfile.write(name + '\n')
            continue



    
def extractAdapterInfoFromFastqcZip(ziplist):
    
    '''
    Take a list of ZIP files that are FASTQC output (as absolute paths) and extract which adapters were found in 
    more than 1% of the reads
    Returns a DICT
    '''
    
    import pandas as pd
    import re
    from zipfile import ZipFile
    
    # create empty dictionary
    d = {}
    
    # loop over zipfiles
    for z in ziplist:

        # extract fastq name
        name = z.split('/')[-1].replace('_fastqc.zip', '')
        
        # open ZIP archive and open file with FASTQC result data
        # split data between all the different analyses of FASTQC
        # (each part ends with the string '>>END_MODULE')
        with ZipFile(z) as myzip:
            with myzip.open(name + '_fastqc/fastqc_data.txt') as myfile:
                fastqc = myfile.read().decode().split('>>END_MODULE')
        
        # extract results for ADAPTERS
        adapter = [e for e in fastqc if re.search('>>Adapter', e)][0].splitlines()[2:]
        
        # split each row by TAB
        adapter = [e.split('\t') for e in adapter]
        
        if not adapter:
            d[name] = ['impossible_to_determine']
            continue
            
        # save to DF, with adapters as column names
        df = pd.DataFrame(adapter[1:], columns=adapter[0], dtype=float)
        
        # create a SERIES containing whether any column contains values > 1
        # the values represent the percent of adapters among all the reads
        s = (df > 1).any()[1:]
        
        # add info to growing dictionary (key=sample name, value=list of adapters)
        d[name] = sorted(list(s[s==True].index))

    return(d)

    