
# coding: utf-8

# In[1]:


import pandas as pd
import sys
import os
import yaml
import subprocess as sp
import argparse
import IPython

sys.path.append("/home/PERSONALE/eugenio.fonzi2")
import script.wes_pipeline_marco_functions as wpm

#### create dir tree
wpm.createDirTree()


##---------------------------------------------------------------------------------##



#### Parse command line arguments

parser = argparse.ArgumentParser()
parser.add_argument('-d', required=True, dest='dataset', action='store', help='whole pathway for dataset file')
parser.add_argument('-c', required=True, dest='config', action='store', help='whole pathway for config file')
parser.add_argument('-r', dest='run_cadaver', action='store_true', help='write and run a cadaver script to import fastq from NAS')
parser.add_argument('-w', dest='write_cadaver', action='store_true', help='write but not run a cadaver script to import fastq from NAS')
parser.add_argument('-a', dest='align', action='store_true', help='perform alignment to genomic and mitochondrial reference')
parser.add_argument('-v', dest='variant_call', action='store_true', help='perform variant calling step and downstream filtering')
parser.add_argument('-m', dest='custom_matching', action='store_true', help='define custom tumor-normal matching')
#parser.add_argument('-o', dest='outdir', action='store', help='specify desired output directory')

# get arguments
args = parser.parse_args()

# get dataset file path
dataset = args.dataset

# get config file path
config = args.config

# get True/False value of 'custom_matching'
run = args.run_cadaver

# get True/False value of 'custom_matching'
write = args.write_cadaver

# get True/False value of 'custom_matching'
align = args.align

# get True/False value of 'custom_matching'
variant = args.variant_call

# get True/False value of 'custom_matching'
match = args.custom_matching

# extract common main path
main = os.path.commonprefix([os.getcwd(), config])



##---------------------------------------------------------------------------------##



#### Load data

## load config file for WES
with open(config) as infile:
    cfg = yaml.load(infile)


## Load list of IDs of desired samples
with open('sample_list.tsv') as infile:
    sample_names = infile.read().splitlines()

## Load AML WES DATASET and set 'sample ids' as indeces
dataset = pd.read_table(dataset, sep = '\t', header=0, index_col='sample_id', dtype='str')

# set pcr primers
pcr1 = cfg['pcr1']
pcr2 = cfg['pcr2']


print('Following samples will be processed:')
print(sample_names)
print('          ')


##---------------------------------------------------------------------------------##


### If option '-w' was selected

if write:
    
    print('------------------------------------------------')
    print('WRITE cadaver script to download all fastq files')
    print('------------------------------------------------')
    
    # write a cadaver script to download all necessary fastq from NAS
    wpm.writeCadaverScript(sample_names, dataset, cfg, main)


    
##---------------------------------------------------------------------------------##



### If option '-r' was selected

if run:
    
    print('--------------------------------------------------------------------')
    print('WRITE and RUN cadaver script to download fastq files for each sample')
    print('--------------------------------------------------------------------')

    # loop over samples
    for name in sample_names:
        
        print('downloading all fastq for sample %s' %name)
        
        # write AND RUN a cadaver script to download from NAS the fastq for a sample
        wpm.writeRunCadaverScript(name, dataset, cfg, main)

        
##---------------------------------------------------------------------------------##


### If option '-a' was selected

if align:

    print('-----------------------------------------------------------')
    print('Align each sample to genomic and mitochondrial reference')
    print('-----------------------------------------------------------')

    for name in sample_names:
    
        print('> ALIGNMENT OF SAMPLE %s' %name)
        
        # get kit for library WES
        kit = dataset.loc[name, 'kit_wes']

        ### decompress fastq.gz
        print('decompressing fastq')
        fastq1, fastq2 = wpm.decompressFastq(name)
        
        ### fastqc
        print('fastqc')
        wpm.fastqc(fastq1)
        wpm.fastqc(fastq2)
        
        ### trim
        print('trimming fastq')
        wpm.trim(fastq1, fastq2, name, pcr1, pcr2)
    
        ### fastqc trimmed
        print('fastqc on trimmed fastq')
        wpm.fastqc(fastq1, on_trimmed=True)
        wpm.fastqc(fastq2, on_trimmed=True)
        
        
        #### GENOME ALIGNMENT
        
        ### alignment
        print('aligning to reference NUCLEAR GENOME')
        wpm.align(fastq1, fastq2, name, cfg['genome']['fastqpath'], cfg['genome']['alignpath'], cfg['genome']['reference'], kit)
        
        ### sorting
        print('sorting')
        wpm.sort(name, cfg['genome']['alignpath'])
        
        ### marking
        print('marking')
        wpm.mark(name, cfg['genome']['alignpath'])
        
        ### indexing
        print('indexing')
        wpm.index(name, cfg['genome']['alignpath'])
        
        ### indel realigning
        print('indel realignment')
        wpm.indelRealign(name, cfg['genome']['alignpath'], cfg['genome']['reference'], cfg['genome']['bed'][kit])
        
        ### bqsr
        print('bqsr')
        wpm.bqsr(name, cfg['genome']['alignpath'], cfg['genome']['reference'], cfg['genome']['bed'][kit])
        
        
        #### EXOME ALIGNMENT
        
        ### alignment
        print('aligning to reference EXOME')
        wpm.align(fastq1, fastq2, name, cfg['genome']['fastqpath'], cfg['exome']['alignpath'], cfg['exome']['reference'][kit], kit)
        
        ### sorting
        print('sorting')
        wpm.sort(name, cfg['exome']['alignpath'])
        
        
        #### EXTRACT UNMAPPED
        print('extracting unmapped reads')
        fastq1, fastq2 = wpm.extractUnmapped(name, cfg['exome']['alignpath'])
        
        
        #### MT ALIGNMENT
        
        ### alignment
        print('aligning unmapped reads to MITOCHONDRIAL GENOME')
        wpm.align(fastq1, fastq2, name, cfg['MT']['fastqpath'], cfg['MT']['alignpath'], cfg['MT']['reference'], kit)
        
        ### sorting
        print('sorting')
        wpm.sort(name, cfg['MT']['alignpath'])
        
        ### marking
        print('marking')
        wpm.mark(name, cfg['MT']['alignpath'])
        
        ### indexing
        print('indexing')
        wpm.index(name, cfg['MT']['alignpath'])
        
        ### indel realigning
        print('indel realignment')
        wpm.indelRealign(name, cfg['MT']['alignpath'], cfg['MT']['reference'], cfg['MT']['bed'])
        
        ### bqsr
        print('bqsr')
        wpm.bqsr(name, cfg['MT']['alignpath'], cfg['MT']['reference'], cfg['MT']['bed'])
        
        ### delete intermediate files
        print('deleting intermediate files')
        wpm.deleteIntermediate(name, step='alignment')
        
        print('                     ')
        print('END OF ALIGNMENT STEP')
        print('                     ')


##---------------------------------------------------------------------------------##



### if option '-v' was selected

if variant:
    
    print('------------------------------------------------')
    print('Call and filter somatic variants for each sample')
    print('------------------------------------------------')

    # create empty list to store tumor_normal name matches 
    name_list = []
    
    # if a custom 'tumor-normal' matching is needed
    if match:
        
        ## load CUSTOM MATCHING table
        matching = pd.read_table('matching.tsv', sep='\t', header=0, dtype='str')
    
        # get list of tumors that are both in sample list and custom list
        tumor = [s for s in sample_names if s in matching.index]
    
        # loop over them
        for t in tumor:
            
            # get normal samples that are matched to t in custom list
            normal = matching.loc[t, 'normal']
            
            # if there is only one normal, concatenate to tumor and add to tumor_normal list
            if type(normal) == str:
                name_list.append('_'.join([t, normal]))
    
            # if there are many normal, concatenate each of them to tumor and add to tumor_normal list
            else:
                for n in normal:
                    name_list.append('_'.join([t, n]))
            
            
    # if no custom matching is needed
    else:
    
        # get list of tumor contained in sample list
        tumor = [s for s in sample_names if dataset.loc[s, 'matching'] == 'tumor']
    
        # loop over them
        for t in tumor:
            
            # get patient ID
            patient = dataset.loc[t, 'patient_id']
        
            # get matched normal
            normal = dataset[(dataset.patient_id == patient) & (dataset.matching == 'normal')].index[0]
        
            # concatenate it to tumor and add to tumor_normal list
            name_list.append('_'.join([t, normal]))
    
    
    print('Following tumor-normal matches will be processed:')
    print(name_list)

    for name in name_list:    
    
        print('> VARIANTS CALL FOR %s MATCH' %name)

        # get kit for library WES
        kit = dataset.loc[name.split('_')[0], 'kit_wes']

        #### MUTECT
        # on genomic BAM
        print('running MuTect on genomic alignment')
        wpm.mutect(name, cfg['genome']['alignpath'], cfg['genome']['mutectpath'], 
                   cfg['genome']['reference'], cfg['genome']['bed'][kit])
        
        # on MT BAM
        print('running MuTect on mitochondrial alignment')
        wpm.mutect(name, cfg['MT']['alignpath'], cfg['MT']['mutectpath'], 
                   cfg['MT']['reference'], cfg['MT']['bed'])
        
        #### VARSCAN
        # on genomic BAM
        print('running VarScan on genomic alignment')
        wpm.varscan(name, cfg['genome']['alignpath'], cfg['genome']['varscanpath'], 
                    cfg['genome']['reference'], cfg['genome']['bed'][kit])
        
        # on MT BAM
        print('running VarScan on mitochondrial alignment')
        wpm.varscan(name, cfg['MT']['alignpath'], cfg['MT']['varscanpath'], 
               cfg['MT']['reference'], cfg['MT']['bed'])
    
    
        # filter somatic variants
        print('filtering somatic variants')
        wpm.filterSomatic(name, cfg['genome']['mutectpath'], cfg['MT']['mutectpath'],
                      cfg['genome']['varscanpath'], cfg['MT']['varscanpath'])
            
        # call annovar and merge [mutect + varscan]
        print('annotating with ANNOVAR')
        wpm.annovar(name)
    
        # apply filter on exonic, non-synonymous, polymorphisms
        print('filtering exonic/non-synonymous/non-polymorphic')
        wpm.filterExonicPolymorphic(name)
        
        #### delete intermediate files
        print('deleting intermediate files')
        wpm.deleteIntermediate(name, step='variant_call')

        print('                     ')        
        print('END OF VARIANT CALL STEP')
        print('                     ')

        
##---------------------------------------------------------------------------------##

print('------------')
print('SESSION INFO')
print('------------')

# print 'session info'
print(IPython.sys_info())

print('            ')
# print all installed modules with their version
print(sp.Popen("pip3 freeze", shell=True, stdout=sp.PIPE).stdout.read().decode().splitlines())

