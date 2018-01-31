

def flagstat(bamlist, samtoolspath):

    ### run FLAGSTAT for each BAM in BAMLIST
    ### bamlist is a list of bam file name with their absolute paths
    
    import subprocess as sp
    
    for bam in bamlist:
    
        # get sample name
        sample = bam.split('/')[-1].replace('.bam', '')
    
        # build command and run
        cmd = ' '.join([samtoolspath + 'samtools', 'flagstat',
                        bam,
                        '>', sample + '_flagstat.txt'])
    
        sp.run(cmd, shell=True)




def samtoolsDepth(bamlist, samtoolspath, bedfile, Q):
    
    ### run SAMTOOLS DEPTH for all BAM in a bam list
    ### bamlist is a list of bam file name with their absolute paths
    ### Q is the mapping quality threshold for read alignments to be considered
    
    import pandas as pd
    import subprocess as sp

    # join all BAM file names into one string
    bamlist_string = ' '.join(bamlist)
    
    # build command and run (capturing STDOUT)
    cmd = ' '.join([samtoolspath + 'samtools', 'depth',
                    '-b', bedfile,
                    '-Q', str(Q),
                    bamlist_string,
                    '> samtools_depth.tsv'])
    
    sp.run(cmd, shell=True)
        
    
    ### extract and summarize coverage data
    
    # get sample names
    sample_names = [b.split('/')[-1].replace('.bam', '') for b in bamlist]

    # transform to column names for coverage table
    colnames = ['chr', 'position'] + sample_names
    
    # load coverage table and
    coverage = pd.read_table('samtools_depth.tsv', sep='\t', header=None, names=colnames)
    
    summary = pd.DataFrame()
    for s in sample_names:
    
        # summarize data for sample and add to DF as new column
        summary = pd.concat([summary, pd.DataFrame(coverage.loc[:,s].describe())], axis=1)
    
    # add sample names as new columns
    summary.columns = sample_names
    
    # delete intermediate file
    sp.run('rm samtools_depth.tsv', shell=True)
    
    return(summary)




def offTargets(bamlist, bedtoolspath, bedfile, plots=False):
    
    '''Transform the BAM to BED, then intersect it to the exome target BED file to identify the read alignments that mapped OFF-TARGET and then compute the percentage of OFF-TARGET.
    OPTION 'plots' --> save plots of the distributions of the mapping qualities of the alignments for:
    - all reads alignment as absolute frequency
    - off-targets as absolute frequency
    - all reads alignment as relative frequency 
    - off-targets as relative frequency
    '''
    
    import pandas as pd
    import subprocess as sp
    
    if plots == True:
        import matplotlib.pyplot as plt

    # open connection to file and write column titles
    with open('read_alignment.log', 'w') as outfile:
        outfile.write('\t'.join(['sample', 'all_read_alignments', '%_off_target']) + '\n')
    
        
    for bam in bamlist:
        
        # get sample name
        sample = bam.split('/')[-1].replace('.bam', '')
    
        print('analyzing off-target for %s' %sample)    
    
        # transform BAM to BED
        cmd = ' '.join([bedtoolspath + 'bamToBed',
                        '-i', bam,
                        '>', sample + '.bed'])
        sp.run(cmd, shell=True)
    
    
        # get read alignments that don't overlap with EXOME BED
        cmd = ' '.join([bedtoolspath + 'intersectBed',
                        '-a', sample + '.bed',
                        '-b', bedfile,
                        '-v',
                        '>', sample + '_off.bed'])
        sp.run(cmd, shell=True)
    
    
        ### COUNT HOW MANY TOTAL AND OFF TARGET READS ARE THERE
        
        # count number of lines in file
        all_reads = sp.run(' '.join(['wc', '-l', sample + '.bed']), shell=True, stdout=sp.PIPE).stdout.decode()
        off_reads = sp.run(' '.join(['wc', '-l', sample + '_off.bed']), shell=True, stdout=sp.PIPE).stdout.decode()
        
        # split and transform string to float
        all_reads = float(all_reads.split(' ')[0])
        off_reads = float(off_reads.split(' ')[0])
    
        
        with open('read_alignment.log', 'a') as outfile:
            outfile.write('\t'.join([sample, str(all_reads), str(off_reads/all_reads*100)]) + '\n')
        

        
        if plots == True:
            
            ### COUNT READ ALIGNMENT BY THEIR MAPPING QUALITY FOR ALL READS AND OFF TARGET

            cmd = 'cut -f5 %s.bed | sort | uniq -c > %s.count' % (sample, sample)
            sp.run(cmd, shell=True)
        
            cmd = 'cut -f5 %s_off.bed | sort | uniq -c > %s_off.count' % (sample, sample)
            sp.run(cmd, shell=True)
            
            
            ### IMPORT COUNT TABLES AND NORMALIZE VALUES FOR ALL READS AND OFF TARGET
            
            count = pd.read_table(sample + '.count', sep='\s', header=None, engine='python')
            count = count.sort_values(1)
            count[2] = [e/sum(count[0]) for e in count[0]]
            
            count_off = pd.read_table(sample + '_off.count', sep='\s', header=None, engine='python')
            count_off = count_off.sort_values(1)
            count_off[2] = [e/sum(count_off[0]) for e in count_off[0]]

            ### SAVE PLOTS DISTRIBUTION OF READS MAPPING QUALITIES
            
            _ = plt.plot(count[1], count[0], '-', linewidth=0.7, label='all reads')
            _ = plt.plot(count_off[1], count_off[0], 'r-', linewidth=0.7, label='off target')
            _ = plt.title('Distribution of read alignments mapping qualities - absolute frequency')
            _ = plt.xlabel('Mapping quality')
            _ = plt.ylabel('Reads count')
            _ = plt.legend()
            plt.savefig(sample + '_mapq_abs.png')
            
            
            _ = plt.plot(count[1], count[2], '-', linewidth=0.7, label='all reads')
            _ = plt.plot(count_off[1], count_off[2], 'r-', linewidth=0.7, label='off target')
            _ = plt.title('Distribution of read alignments mapping qualities - relative frequency')
            _ = plt.xlabel('Mapping quality')
            _ = plt.ylabel('Reads proportion')
            plt.legend()
            plt.savefig(sample + '_mapq_rel.png')
    
    
        # remove intermediate files
        sp.run('rm %s*bed' %sample, shell=True)
        sp.run('rm *.count', shell=True)
        
        
        

def countAlignments(bamlist, samtoolspath, names=None):
    
    '''
    Take a list of BAM/SAM directories and count the number of alignments inside each.
    The results are stored into a table.
    
    Names to be used for each BAM/SAM can be specified by providing a list through 
    the option 'names'. The order must be the same as 'bamlist'.
    '''
    
    import subprocess as sp

    # open connection to file and write column names
    with open('alignment_counts.tsv', 'w') as outfile:
        outfile.write('file_name\tnumber_of_alignments\n')


    # loop over files
    for i in range(len(bamlist)):
        
        # get BAM path   
        bam = bamlist[i]
        
        # get name of BAM/SAM
        if names:
            name = names[i]
        else:
            name = bam.split('/')[-1][:-4]
    
        # count number of alignments (i.e. number of lines excluding the header) 
        cmd = '%s/samtools view %s | wc -l' %(samtoolspath, bam)
        lines = sp.run(cmd, shell=True, stdout=sp.PIPE).stdout.decode()
        
        # store file name and #alignments to file
        with open('alignment_counts.tsv', 'a') as outfile:
            outfile.write('%s\t%s' %(name, lines))


    


def compareBamPairwise(bamlist, samtoolspath, names=None, tmp=''):

    '''
    Take a list of BAM/SAM directories and compare all of them pairwise.
    For each pair, a .diff file is produced. If the file is empty, it means that
    all alignments were equal. Otherwise, the different alignments are stored
    into .diff.
    Only the first 9 fields of the BAM/SAM are compared.
    
    A list of 'names' to be used to identify each BAM/SAM can be provided and its
    order must match the order of paths in 'bamlist'.
    
    The 'sort' command requires a lot of disk space and default '/tmp/' folder could fill up;
    in this case is possible to specify a custom temporary directory (which must be previously created).
    '''    
    
    import subprocess as sp
    import itertools as it
    
    # in case a custom temporary directory was specified
    if tmp != '':
        tmp = '-T %s ' %tmp
    
    # loop over indeces
    for i in range(len(bamlist)):

        # get path and name of BAM/SAM    
        bam = bamlist[i]

        # get name of BAM/SAM        
        if names:
            name = names[i]
        else:
            name = bam.split('/')[-1][:-4]
        
        # keep only fields from 1 to 9 and save to file
        cmd = '%s/samtools view %s | cut -f1,2,3,4,5,6,7,8,9 > %s.cut1-9' %(samtoolspath, bam, name)
        sp.run(cmd, shell=True)


    # loop over file combinations (2 elements, no repetition)
    for c in it.combinations(range(len(bamlist)), 2):

        if names:
            # get names of the two samples in the combination
            name1 = names[c[0]]
            name2 = names[c[1]]
        else:
            name1 = bamlist[c[0]].split('/')[-1][:-4]            
            name2 = bamlist[c[1]].split('/')[-1][:-4]
    

        # merge the two files, sort and save the unique lines (alignments that differ)
        cmd = 'cat %s.cut1-9 %s.cut1-9 | sort %s| uniq -u > %s_%s.diff' %(name1, name2, tmp, name1, name2)
        sp.run(cmd, shell=True)  

