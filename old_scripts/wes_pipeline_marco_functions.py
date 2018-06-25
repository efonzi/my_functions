# coding: utf-8


    
def createDirTree():
    
    import os
    
    # save tree structure as dictionary
    dirs = {'01_fastq': ['fastqc'],
            '02_fastq_trimmed': ['fastqc', 'logs'],
            '03_alignment_genome': ['01_intermediate', '02_bqsr', 'logs'],
            '04_alignment_exome': ['01_intermediate', '02_unmapped_fastq', 'logs'], 
            '05_alignment_MT': ['01_intermediate', '02_bqsr', 'logs'], 
            '06_mutect_genome': ['logs'],
            '07_mutect_MT': ['logs'],
            '08_varscan_genome': ['01_mpileup', 'logs'],
            '09_varscan_MT': ['01_mpileup', 'logs'], 
            '10_variants_filter': ['logs'], 
           }
    
    # loop over dictionary items and make each folder and subfolder
    [[os.makedirs('/'.join([d, d2]), exist_ok=True) for d2 in dirs[d]] for d in dirs]


    
    
def writeCadaverScript(sample_names, dataset, cfg, main):

########## write a script for cadaver to download the the fastq from NAS to desired local folder

    import os
    import re
    
    pwd = os.getcwd()
    outpath = pwd + '/01_fastq/'
    
    # open file connection
    with open('cadaver_get.sh', 'w') as outfile:

        for name in sample_names:

            # extract the path of fastq files from 'dataset' and 'cfg'
            fastq_path_nas = cfg[dataset.loc[name, 'fastq_position_wes']]
            
            # get list of all files present in the path
            all_fastq = os.listdir(main + 'NAS/' + fastq_path_nas)
            
            # select only the fastq files whose name contain 'fastq_id_wes' of current sample
            fastq_list = [f for f in all_fastq if re.search(dataset.loc[name, 'fastq_id_wes'], f)]    


            # loop over fastq files of current sample
            for fastq_name in fastq_list:
            
                # add a 'get' command for each fastq file, each in a different line
                outfile.write(' '.join(['get', 
                                        fastq_path_nas + fastq_name,
                                        outpath + fastq_name + '\n']))
            
        # add last line to close connection
        outfile.write('quit')



def writeRunCadaverScript(name, dataset, cfg, main):
    
    import os
    import re
    import subprocess as sp
    
    pwd = os.getcwd()
    outpath = pwd + '/01_fastq/'
    
    # extract the path of fastq files from 'dataset' and 'cfg'
    fastq_path_nas = cfg[dataset.loc[name, 'fastq_position_wes']]
    
    # get list of all files present in the path
    all_fastq = os.listdir(main + 'NAS/' + fastq_path_nas)
    
    # select only the fastq files whose name contain 'fastq_id_wes' of current sample
    fastq_list = [f for f in all_fastq if re.search(dataset.loc[name, 'fastq_id_wes'], f)] 
        
    # open file connection
    with open('cadaver_get.sh', 'w') as outfile:

        # loop over fastq files of current sample
        for fastq_name in fastq_list:
        
            # add a 'get' command for each fastq file, each in a different line
            outfile.write(' '.join(['get', 
                                        fastq_path_nas + fastq_name,
                                        outpath + fastq_name + '\n']))
            
        # add last line to close connection
        outfile.write('quit')
        
    # run cadaver script
    sp.call(' '.join(['cadaver',
                          'https://ngs-ptl.unibo.it:5006',
                          '-r cadaver_get.sh']), shell=True)

    # delete script
    sp.call('rm cadaver_get.sh', shell=True)



def decompressFastq(name):

    import os
    import re
    import subprocess as sp
    
    outpath = '01_fastq/'
        
    # list all downloaded fastq.gz files
    fastq_gz = [f for f in os.listdir(outpath) if f[-3:] == '.gz']
    
    # loop over them
    for gz in fastq_gz:
        
        # decompress it
        sp.call(' '.join(['gzip -d', outpath + gz]), shell=True)
    
    # list and join names of all decompressed R1 and R2 fastq files
    R1 = ' '.join([outpath + f for f in os.listdir(outpath) if re.search('R1', f)])
    R2 = ' '.join([outpath + f for f in os.listdir(outpath) if re.search('R2', f)])
    
    # create fastq names
    fastq1 = name + '_R1.fastq'
    fastq2 = name + '_R2.fastq'
    
    # concatenate multiple fastq (if there are)
    sp.call(' '.join(['cat', R1, '>', outpath + fastq1]), shell=True)
    sp.call(' '.join(['cat', R2, '>', outpath + fastq2]), shell=True)  

    return fastq1, fastq2



def fastqc(fastq, on_trimmed=False):
    
    import subprocess as sp

    if on_trimmed == False:
        inpath = '01_fastq/'
        outpath = '01_fastq/fastqc/'
    
    else:
        inpath = '02_fastq_trimmed/'
        outpath = '02_fastq_trimmed/fastqc/'

    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/FastQC/fastqc'

    
    
    cmd = ' '.join([toolpath,
                    '-o', outpath,
                    inpath + fastq])

    sp.call(cmd, shell=True)





def trim(fastq1, fastq2, name, pcr1, pcr2):
    
    import subprocess as sp
    
    inpath = '01_fastq/'
    outpath = '02_fastq_trimmed/'
    logpath = '02_fastq_trimmed/logs/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/AdapterRemoval'
    
    cmd = ' '.join([toolpath,
                  '--file1', inpath + fastq1, 
                  '--file2', inpath + fastq2,
                  '--pcr1', pcr1,
                  '--pcr2', pcr2, 
                  '--stats',
                  '--trimns',
                  '--trimqualities', 
                  '--minquality 20', 
                  '--minlength 80', 
                  '--output1', outpath + fastq1, 
                  '--output2', outpath + fastq2, 
                  '--discarded', logpath + name + '_discarded.log', 
                  '--outputstats', logpath + name + '_stats.log',
                  '--singleton', logpath + name + '_singleton.log'])
    
    sp.call(cmd, shell=True)


    
def align(fastq1, fastq2, name, fastqpath, alignpath, reference, kit):
    
    import subprocess as sp
        
    inpath = fastqpath
    outpath = alignpath + '01_intermediate/'
    logpath = alignpath + 'logs/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/bwa-0.7.12/'
    refpath = '/home/PERSONALE/eugenio.fonzi2/reference/'

    RG = '\'@RG\\tID:%s\\tSM:%s\\tPL:illumina\\tLB:%s\\tPU:PE\''%(name, name, kit)
    
    cmd = ' '.join([toolpath + 'bwa mem -t 20 -M',
                    '-R', RG,
                    refpath + reference,
                    inpath + fastq1,
                    inpath + fastq2,
                    '>', outpath + name + '.sam', 
                    '2>', logpath + name + '_alignment.log'])
    
    sp.call(cmd, shell=True)




    
def sort(name, alignpath):
    
    import subprocess as sp
    
    inoutpath = alignpath + '01_intermediate/'
    logpath = alignpath + 'logs/'
    tmppath = '/home/PERSONALE/eugenio.fonzi2/tmp/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/picard-tools-1.119/'
    
    cmd = ' '.join(['java -Djava.io.tmpdir=' + tmppath,
                    '-jar', toolpath + 'SortSam.jar',
                    'INPUT=' + inoutpath + name + '.sam',
                    'OUTPUT=' + inoutpath + name + '.sorted.bam',
                    'SORT_ORDER=coordinate',
                    'TMP_DIR=' + tmppath,
                    '2>', logpath + name + '_sorting.log'])
    
    sp.call(cmd, shell=True)
    
    

    
def mark(name, alignpath):
    
    import subprocess as sp
    
    inoutpath = alignpath + '01_intermediate/'
    logpath = alignpath + 'logs/'
    tmppath = '/home/PERSONALE/eugenio.fonzi2/tmp/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/picard-tools-1.119/'
    
    cmd = ' '.join(['java -Djava.io.tmpdir=' + tmppath, 
                    '-jar', toolpath + 'MarkDuplicates.jar',
                    'INPUT=' + inoutpath + name + '.sorted.bam',
                    'OUTPUT=' + inoutpath + name + '.marked.bam',
                    'METRICS_FILE=' + logpath + name + '_metrix.log',
                    'TMP_DIR=' + tmppath,
                    'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000',
                    '2>', logpath + name + '_marking.log'])

    sp.call(cmd, shell=True)
    




def index(name, alignpath):
    
    import subprocess as sp
    
    inoutpath = alignpath + '01_intermediate/'
    logpath = alignpath + 'logs/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/picard-tools-1.119/'
    
    cmd = ' '.join(['java -jar', toolpath + 'BuildBamIndex.jar',
                    'INPUT=' + inoutpath + name + '.marked.bam',
                    '2>', logpath + name + '_indexing.log'])

    sp.call(cmd, shell=True)
    
    
    
    
    
    
def indelRealign(name, alignpath, reference, bed):
    
    import subprocess as sp
    
    inoutpath = alignpath + '01_intermediate/'
    logpath = alignpath + 'logs/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/GenomeAnalysisTK-3.4-0/'
    refpath = '/home/PERSONALE/eugenio.fonzi2/reference/'
    
    cmd = ' '.join(['java -jar', toolpath + 'GenomeAnalysisTK.jar',
                    '-T RealignerTargetCreator',
                    '-nt 2',
                    '-R', refpath + reference,
                    '-I', inoutpath + name + '.marked.bam',
                    '-L', refpath + bed,
                    '-ip 50',
                    '-known', refpath + 'Mills_and_1000G_gold_standard.indels.b37.vcf',
                    '-o', logpath + name + '_target_intervals.list',
                    '2>', logpath + name + '_indelRealigning_01.log'])
    
    sp.call(cmd, shell=True)
    
    
    cmd = ' '.join(['java -jar', toolpath + 'GenomeAnalysisTK.jar',
                    '-T IndelRealigner',
                    '-R', refpath + reference,
                    '-I', inoutpath + name + '.marked.bam',
                    '-targetIntervals', logpath + name + '_target_intervals.list',
                    '-ip 50',
                    '-known', refpath + 'Mills_and_1000G_gold_standard.indels.b37.vcf',
                    '-o', inoutpath + name + '.indelRealigned.bam',
                    '2>', logpath + name + '_indelRealigning_02.log'])
    
    sp.call(cmd, shell=True)

    
    

    
def bqsr(name, alignpath, reference, bed):
    
    import subprocess as sp

    inpath = alignpath + '01_intermediate/'
    outpath = alignpath + '02_bqsr/'
    logpath = alignpath + 'logs/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/GenomeAnalysisTK-3.4-0/'
    refpath = '/home/PERSONALE/eugenio.fonzi2/reference/'

    cmd = ' '.join(['java -jar', toolpath + 'GenomeAnalysisTK.jar',
                    '-T BaseRecalibrator',
                    '-nct 2',
                    '-R', refpath + reference,
                    '-I', inpath + name + '.indelRealigned.bam',
                    '-L', refpath + bed,
                    '-ip 50',
                    '-knownSites', refpath + 'dbsnp_138.b37.vcf',
                    '-knownSites', refpath + 'Mills_and_1000G_gold_standard.indels.b37.vcf',
                    '-o', logpath + name + '_recal_data.table',
                    '2>', logpath + name + '_recalibrating_01.log'])
    
    sp.call(cmd, shell=True)

    
    cmd = ' '.join(['java -jar', toolpath + 'GenomeAnalysisTK.jar',
                    '-T PrintReads',
                    '-nct 2',
                    '-R', refpath + reference,
                    '-I', inpath + name + '.indelRealigned.bam',
                    '-BQSR', logpath + name + '_recal_data.table',
                    '-o', outpath + name + '.bam',
                    '2>', logpath + name + '_recalibrating_02.log'])
    
    sp.call(cmd, shell=True)
    
    

def extractUnmapped(name, alignpath):
    
    import subprocess as sp

    # extract unmapped reads
    inoutpath = alignpath + '01_intermediate/'
    logpath = alignpath + 'logs/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/samtools-1.5/bin/'
    refpath = '/home/PERSONALE/eugenio.fonzi2/reference/'
    
    cmd = ' '.join([toolpath + 'samtools view -b -f 4',
                    inoutpath + name + '.sorted.bam',
                    '>', inoutpath + name + '.unmapped.bam',
                    '2>', logpath + name + '_unmapped.log'])

    sp.call(cmd, shell=True)

                    
    # convert unmapped reads BAM to FASTQ
    inpath = alignpath + '01_intermediate/'
    outpath = alignpath + '02_unmapped_fastq/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/picard-tools-1.119/'

    cmd = ' '.join(['java -jar', toolpath + 'SamToFastq.jar',
                    'I=', inpath + name + '.unmapped.bam',
                    'F=', outpath + name + '_R1.unmapped.fastq',
                    'F2=', outpath + name + '_R2.unmapped.fastq',
                    'FU=', outpath + name + '_unpaired.unmapped.fastq',
                    '2>', logpath + name + '_unmapped_fastq.log'])    
    
    sp.call(cmd, shell=True)
    
    
    fastq1 = name + '_R1.unmapped.fastq'
    fastq2 = name + '_R2.unmapped.fastq'
    
    return fastq1, fastq2
    
    
    
def mutect(name, alignpath, mutectpath, reference, bed):

    import subprocess as sp

    tumor, normal = name.split('_')

    inpath = alignpath + '02_bqsr/'
    outpath = mutectpath
    logpath = mutectpath + 'logs/'
    javapath = '/home/PERSONALE/eugenio.fonzi2/tools/jdk1.6.0_35/jre/bin/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/mutect/'
    refpath = '/home/PERSONALE/eugenio.fonzi2/reference/'
    
    cmd = ' '.join([javapath + 'java -Xmx2g -jar',
                    toolpath + 'muTect-1.1.4.jar',
                    '--analysis_type MuTect',
                    '--reference_sequence', refpath + reference,
                    '--cosmic', refpath + 'b37_cosmic_v54_120711.vcf',
                    '--dbsnp', refpath + 'dbsnp_138.b37.vcf',
                    '--intervals', refpath + bed,
                    '--input_file:normal', inpath + normal + '.bam',
                    '--input_file:tumor', inpath + tumor + '.bam',
                    '--out', outpath + name + '.tsv',
                    '--coverage_file', outpath + name + '_cov_wig.log',
                    '>', logpath + name + '_mutect.log'])
    
    sp.call(cmd, shell=True)
    


def varscan(name, alignpath, varscanpath, reference, bed):
    
    import subprocess as sp

    tumor, normal = name.split('_')
    
    # build mpileup from BQSR BAM
    inpath = alignpath + '02_bqsr/'
    outpath = varscanpath + '01_mpileup/'
    logpath = varscanpath + 'logs/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/samtools-1.5/bin/'
    refpath = '/home/PERSONALE/eugenio.fonzi2/reference/'

    # for TUMOR
    cmd = ' '.join([toolpath + 'samtools mpileup -B -q 1',
                    '-f', refpath + reference,
                    '--positions', refpath + bed,
                    inpath + tumor + '.bam',
                    '>', outpath + tumor + '.mpileup',
                    '2>', logpath + tumor + '_mpileup.log'])

    sp.call(cmd, shell=True)

    
    # for NORMAL
    cmd = ' '.join([toolpath + 'samtools mpileup -B -q 1',
                    '-f', refpath + reference,
                    '--positions', refpath + bed,
                    inpath + normal + '.bam',
                    '>', outpath + normal + '.mpileup',
                    '2>', logpath + normal + '_mpileup.log'])

    sp.call(cmd, shell=True)


    # call VARSCAN
    inpath = varscanpath + '01_mpileup/'
    outpath = varscanpath
    logpath = varscanpath + 'logs/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/'
    refpath = '/home/PERSONALE/eugenio.fonzi2/reference/'

    cmd = ' '.join(['java -jar', toolpath + 'VarScan.v2.3.9.jar somatic',
                    inpath + normal + '.mpileup',
                    inpath + tumor + '.mpileup',
                    '--output-snp', outpath + name + '_snv.tsv',
                    '--output-indel', outpath + name + '_indel.tsv',
                    '--min-avg-qual 15',
                    '--strand_filter 1',
                    '--min-var-freq 0.05',
                    '--somatic-p-value 0.05',
                    '2>', logpath + name + '_varscan.log'])
                    
    sp.call(cmd, shell=True)
    

    
def filterSomatic(name, mutectpath, mutectpath_MT, varscanpath, varscanpath_MT):
    
    import pandas as pd
    import subprocess as sp
    
    toolpath = '/home/PERSONALE/eugenio.fonzi2/tools/'
    outpath = '10_variants_filter/'
    logpath = '10_variants_filter/logs/'
    
    # load mutect and varscan variants, for both genome and MT
    mutect = pd.read_table(mutectpath + name + '.tsv', sep='\t', header=1, dtype='str')
    mutect_MT = pd.read_table(mutectpath_MT + name + '.tsv', sep='\t', header=1, dtype='str')
    
    varscan_snp = pd.read_table(varscanpath + name + '_snv.tsv', sep='\t', header=0, dtype='str')
    varscan_snp_MT = pd.read_table(varscanpath_MT + name + '_snv.tsv', sep='\t', header=0, dtype='str')

    varscan_indel = pd.read_table(varscanpath + name + '_indel.tsv', sep='\t', header=0, dtype='str')
    varscan_indel_MT = pd.read_table(varscanpath_MT + name + '_indel.tsv', sep='\t', header=0, 
                                     dtype='str')
    
    
    # merge genome and MT variants
    mutect = pd.concat([mutect, mutect_MT], ignore_index=True)
    
    varscan_snp = pd.concat([varscan_snp, varscan_snp_MT], ignore_index=True)
    
    varscan_indel = pd.concat([varscan_indel, varscan_indel_MT], ignore_index=True)
    
    
    # merge all varscan variants
    varscan_all = pd.concat([varscan_snp, varscan_indel], ignore_index=True)
    
    
    # write to file UNFILTERED variants for mutect and varscan(snv + indel)
    mutect.to_csv(outpath + name + '_mutect.tsv', sep='\t', header=True, index=False)

    varscan_snp.to_csv(outpath + name + '_varscan_snv_temp1.tsv', sep='\t', header=True, index=False)
    varscan_indel.to_csv(outpath + name + '_varscan_indel_temp1.tsv', sep='\t', header=True, index=False)
    varscan_all.to_csv(outpath + name + '_varscan.tsv', sep='\t', header=True, index=False)


    # remove variants that did not pass mutect requirements
    mutect = mutect[(mutect.judgement == 'KEEP') & (mutect.covered == 'COVERED')]


    # apply somaticFilter to varscan SNV
    cmd = ' '.join(['java -jar', toolpath + 'VarScan.v2.3.9.jar somaticFilter',
                    outpath + name + '_varscan_snv_temp1.tsv', 
                    '--indel-file', outpath + name + '_varscan_indel_temp1.tsv',
                    '--min-coverage 1', 
                    '--min-reads2 2', 
                    '--min-var-freq 0.1', 
                    '--output-file', outpath + name + '_varscan_snv_temp2.tsv',
                    '2>', logpath + name + '_somaticFilter_snv_err.log'])
    
    sp.call(cmd, shell=True)
    
    
    # apply somaticFilter to varscan INDEL
    cmd = ' '.join(['java -jar', toolpath + 'VarScan.v2.3.9.jar somaticFilter',
                    outpath + name + '_varscan_indel_temp1.tsv',
                    '--min-coverage 1', 
                    '--min-reads2 2', 
                    '--min-var-freq 0.1', 
                    '--output-file', outpath + name + '_varscan_indel_temp2.tsv',
                    '2>', logpath + name + '_somaticFilter_indel_err.log'])
    
    sp.call(cmd, shell=True)
    
    
    # laod and merge filtered varscan tables
    varscan_snp = pd.read_table(outpath + name + '_varscan_snv_temp2.tsv', sep='\t', header=0)
    varscan_indel = pd.read_table(outpath + name + '_varscan_indel_temp2.tsv', sep='\t', header=0)

    varscan_all = pd.concat([varscan_snp, varscan_indel], ignore_index=True)
    
    
    # keep only variants that varscan called 'Somatic' or 'LOH'
    varscan_all = varscan_all[(varscan_all.somatic_status == 'Somatic') | (varscan_all.somatic_status == 'LOH')]
    
    # keep only variants with variant reads on both strands (for tumor)
    varscan_all = varscan_all[(varscan_all.tumor_reads2_plus > 0) & (varscan_all.tumor_reads2_minus > 0)]
    
    # write FILTERED mutect and varscan variants
    mutect.to_csv(outpath + name + '_mutect_somatic.tsv', sep='\t', header=True, index=False)
    varscan_all.to_csv(outpath + name + '_varscan_somatic.tsv', sep='\t', header=True, index=False)

    # delete temporary varscan files
    sp.call(' '.join(['rm', outpath + name + '_varscan_snv_temp1.tsv']), shell=True)
    sp.call(' '.join(['rm', outpath + name + '_varscan_indel_temp1.tsv']), shell=True)
    sp.call(' '.join(['rm', outpath + name + '_varscan_snv_temp2.tsv']), shell=True)
    sp.call(' '.join(['rm', outpath + name + '_varscan_indel_temp2.tsv']), shell=True)




# for next function
#    varscan_all[varscan_all.chrom == 'X'] = '23'
#    varscan_all[varscan_all.chrom == 'Y'] = '24'


def annovar(name):
    
    import subprocess as sp
    import pandas as pd
    
    inoutpath = '10_variants_filter/'
    logpath = '10_variants_filter/logs/'
    toolpath = '/home/PERSONALE/eugenio.fonzi2/annovar/'
        

    ### MUTECT

    # load lists of somatic filtered variants for mutect and varscan
    mutect = pd.read_table(inoutpath + name + '_mutect_somatic.tsv', sep='\t', header=0)

    # extract chr, start, end, ref, alt columns
    mutect_anno = pd.concat([mutect[['contig', 'position']], 
                             mutect['position'], 
                             mutect[['ref_allele', 'alt_allele']]
                            ], axis=1)
    
    # split genome and mt variants
    mutect_anno_genome = mutect_anno[mutect_anno.contig != 'MT'].copy()
    mutect_anno_mt = mutect_anno[mutect_anno.contig == 'MT'].copy()
    
    # write to file
    mutect_anno_genome.to_csv(inoutpath + name + '_mutect_annovar_genome.tsv', sep='\t', header=False, index=False)
    mutect_anno_mt.to_csv(inoutpath + name + '_mutect_annovar_mt.tsv', sep='\t', header=False, index=False)
    
    # call annovar for genome variants
    cmd_genome = ' '.join(['perl', toolpath + 'table_annovar.pl',
                    inoutpath + name + '_mutect_annovar_genome.tsv',
                    toolpath + 'humandb/',
                    '-buildver hg19',
                    '-out',  logpath + name + '_mutect_genome',
                    '-remove',
                    '-protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,snp138,1000g2014oct_all,exac03nontcga,ljb26_all',
                    '-operation g,r,r,f,f,f,f,f',
                    '2>', logpath + name + '_mutect_genome_annovar_err.log'])

    # call annovar for MT variants (only gene annotation, on 'ensemble')
    cmd_mt_g = ' '.join(['perl', toolpath + 'table_annovar.pl',
                    inoutpath + name + '_mutect_annovar_mt.tsv',
                    toolpath + 'humandb/',
                    '-buildver GRCh37_MT',
                    '-out',  logpath + name + '_mutect_mt_g',
                    '-remove',
                    '-protocol ensGene',
                    '-operation g',
                    '2>', logpath + name + '_mutect_mt_g_annovar_err.log'])

    sp.call(cmd_genome, shell=True)
    sp.call(cmd_mt_g, shell=True)

    
    # load annotated file
    mutect_anno_genome = pd.read_table(logpath + name + '_mutect_genome.hg19_multianno.txt', sep='\t', header=0)
    mutect_anno_mt_g = pd.read_table(logpath + name + '_mutect_mt_g.GRCh37_MT_multianno.txt', sep='\t', header=0)
    
    # select common columns
    mutect_anno_genome = mutect_anno_genome[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
                                             'ExonicFunc.refGene', 'AAChange.refGene', 'cytoBand', 'genomicSuperDups',
                                             'esp6500siv2_all', 'snp138', '1000g2014oct_all', 'ExAC_nontcga_ALL', 
                                             'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 
                                             'MutationAssessor_pred']]

    mutect_anno_mt_g = mutect_anno_mt_g[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.ensGene', 'Gene.ensGene',
                                             'ExonicFunc.ensGene', 'AAChange.ensGene']]
    
    # compute the difference of columns
    cols_diff = len(mutect_anno_genome.columns) - len(mutect_anno_mt_g.columns)
    
    # add to MT annovar table the missing columns
    mutect_anno_mt_g = pd.concat([mutect_anno_mt_g, pd.DataFrame(columns=['x']*cols_diff)], axis=1)
    
    # assign same column name
    mutect_anno_mt_g.columns = mutect_anno_genome.columns
    
    # merge by row
    mutect_anno = pd.concat([mutect_anno_genome, mutect_anno_mt_g], ignore_index=True)
    
    # merge total annovar output with mutect output by column
    mutect = pd.concat([mutect_anno,
                        pd.DataFrame({'somatic_status': [float('nan') * len(mutect)]}),
                        pd.DataFrame({'n_vaf': mutect['n_alt_count'] / (mutect['n_ref_count'] + mutect['n_alt_count'])}),
                        pd.DataFrame({'n_depth': mutect['n_ref_count'] + mutect['n_alt_count']}),
                        pd.DataFrame({'t_vaf': mutect['t_alt_count'] / (mutect['t_ref_count'] + mutect['t_alt_count'])}),
                        pd.DataFrame({'t_depth': mutect['t_ref_count'] + mutect['t_alt_count']}),
                        pd.DataFrame({'detection_method': ['mutect'] * len(mutect)}),
                       ], axis=1)
    
    
    
    ### VARSCAN

    # load lists of somatic filtered variants for mutect and varscan
    varscan = pd.read_table(inoutpath + name + '_varscan_somatic.tsv', sep='\t', header=0)

    # convert 'tumor_var_freq' from percentage to proportion
    varscan['tumor_var_freq'] = [float(freq[:-1])/100 for freq in varscan['tumor_var_freq']]
    
    # extract chr, start, end, ref, alt columns (renaming 'end' one)
    varscan_anno = pd.concat([varscan[['chrom', 'position']], 
                     pd.DataFrame({'end': varscan['position']}), 
                     varscan[['ref', 'var']]], axis=1)

    # get indeces of rows with DELETIONS and INSERTIONS
    deletion = varscan_anno[['-' in v for v in varscan_anno['var']]].index
    insertion = varscan_anno[['+' in v for v in varscan_anno['var']]].index
    
    ## INSERTIONS
    # change 'ref' column to '-'
    varscan_anno.loc[insertion, 'ref'] = ['-'] * len(insertion)
    # remove '-' from the beginning of strings in 'var' column
    varscan_anno.loc[insertion, 'var'] = [v[1:] for v in varscan_anno.loc[insertion, 'var']]
    # add 1 to 'start' and 'end' chromosomal positions
    varscan_anno.loc[insertion, 'position'] = varscan_anno.loc[insertion, 'position'] + 1
    varscan_anno.loc[insertion, 'end'] = varscan_anno.loc[insertion, 'end'] + 1
    
    ## DELETIONS
    # substitute strings in 'ref' with strings in 'var', removing the '-' at the beginning
    varscan_anno.loc[deletion, 'ref'] = [v[1:] for v in varscan_anno.loc[deletion, 'var']]
    # change 'var' column to '-'
    varscan_anno.loc[deletion, 'var'] = ['-'] * len(deletion)
    # add 1 to 'start' chromosomal positions
    varscan_anno.loc[deletion, 'position'] = varscan_anno.loc[deletion, 'position'] + 1
    # add length of string in 'ref' to 'end' chromosomal positions
    varscan_anno.loc[deletion, 'end'] = [varscan_anno.loc[i, 'end'] + 
                                         len(varscan_anno.loc[i, 'ref']) for i in deletion]

    # split between genomic variants and MT
    varscan_anno_genome = varscan_anno[varscan_anno.chrom != 'MT'].copy()
    varscan_anno_mt = varscan_anno[varscan_anno.chrom == 'MT'].copy()
    
    
    # write to file
    varscan_anno_genome.to_csv(inoutpath + name + '_varscan_annovar_genome.tsv', sep='\t', header=False, index=False)
    varscan_anno_mt.to_csv(inoutpath + name + '_varscan_annovar_mt.tsv', sep='\t', header=False, index=False)

    
    # call annovar for genomic variants    
    cmd_genome = ' '.join(['perl', toolpath + 'table_annovar.pl',
                    inoutpath + name + '_varscan_annovar_genome.tsv',
                    toolpath + 'humandb/',
                    '-buildver hg19',
                    '-out',  logpath + name + '_varscan_genome',
                    '-remove',
                    '-protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,snp138,1000g2014oct_all,exac03nontcga,ljb26_all',
                    '-operation g,r,r,f,f,f,f,f',
                    '2>', logpath + name + '_varscan_genome_annovar_err.log'])
        
    # call annovar for MT variants (only gene annotation, on 'ensemble')
    cmd_mt_g = ' '.join(['perl', toolpath + 'table_annovar.pl',
                    inoutpath + name + '_varscan_annovar_mt.tsv',
                    toolpath + 'humandb/',
                    '-buildver GRCh37_MT',
                    '-out',  logpath + name + '_varscan_mt_g',
                    '-remove',
                    '-protocol ensGene',
                    '-operation g',
                    '2>', logpath + name + '_varscan_mt_g_annovar_err.log'])
 
    sp.call(cmd_genome, shell=True)
    sp.call(cmd_mt_g, shell=True)

    
    # load annotated files
    varscan_anno_genome = pd.read_table(logpath + name + '_varscan_genome.hg19_multianno.txt', sep='\t', header=0)
    varscan_anno_mt_g = pd.read_table(logpath + name + '_varscan_mt_g.GRCh37_MT_multianno.txt', sep='\t', header=0)
    
    # select common columns
    varscan_anno_genome = varscan_anno_genome[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
                                             'ExonicFunc.refGene', 'AAChange.refGene', 'cytoBand', 'genomicSuperDups',
                                             'esp6500siv2_all', 'snp138', '1000g2014oct_all', 'ExAC_nontcga_ALL', 
                                             'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 
                                             'MutationAssessor_pred']]

    varscan_anno_mt_g = varscan_anno_mt_g[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.ensGene', 'Gene.ensGene',
                                             'ExonicFunc.ensGene', 'AAChange.ensGene']]
    
    # compute the difference of columns
    cols_diff = len(varscan_anno_genome.columns) - len(varscan_anno_mt_g.columns)
    
    # add to MT annovar table the missing columns
    varscan_anno_mt_g = pd.concat([varscan_anno_mt_g, pd.DataFrame(columns=['x']*cols_diff)], axis=1)
    
    # assign same column name
    varscan_anno_mt_g.columns = varscan_anno_genome.columns
    
    # merge by row
    varscan_anno = pd.concat([varscan_anno_genome, varscan_anno_mt_g], ignore_index=True)
    
    # merge total annovar output with varscan output by column
    varscan = pd.concat([varscan_anno,
                         varscan['somatic_status'],
                         pd.DataFrame({'n_vaf': varscan['normal_reads2']/(varscan['normal_reads1'] + varscan['normal_reads2'])}),
                         pd.DataFrame({'n_depth': varscan['normal_reads1'] + varscan['normal_reads2']}),
                         pd.DataFrame({'t_vaf': varscan['tumor_reads2']/(varscan['tumor_reads1'] + varscan['tumor_reads2'])}),
                         pd.DataFrame({'t_depth': varscan['tumor_reads1'] + varscan['tumor_reads2']}),
                         pd.DataFrame({'detection_method': ['varscan'] * len(varscan)}),
                        ], axis=1)
    

    # set common column names
    colnames = ['chr', 'start', 'end', 'ref', 'alt', 'exon', 'gene', 'var_type', 'aa_change', 'cytoband', 'segdups', 'esp', 
                'dbsnp_id', '1000g', 'exac', 'sift_score', 'sift_pred', 'poly_pred', 'mutasses_pred', 'somatic_status', 'n_vaf', 
                'n_depth', 't_vaf', 't_depth', 'detection_method']
    
    mutect.columns = colnames
    varscan.columns = colnames
    
    
    
    #### MERGE MUTECT WITH VARSCAN
    
    # set columns to merge
    colnames = ['chr', 'start', 'end', 'ref', 'alt']
    # concatenate the strings in these columns and save into columns 'variant'
    mutect['variant'] = ['_'.join([str(e) for e in mutect.loc[i, colnames]]) for i in mutect.index]
    varscan['variant'] = ['_'.join([str(e) for e in varscan.loc[i, colnames]]) for i in varscan.index]

    # set column 'variant' as index
    mutect.index = mutect.variant
    varscan.index = varscan.variant
    
    # identify shared and unique variants between mutect and varscan lists
    unique_mutect = [v for v in mutect.variant if v not in varscan.variant]
    unique_varscan = [v for v in varscan.variant if v not in mutect.variant]
    common = [v for v in varscan.variant if v in mutect.variant]

    # if a variant in varscan is in 'common', change its 'detection_method' field to 'both'
    varscan.loc[common, 'detection_method'] = 'both'
    
    # remove common variants from mutect
    mutect = mutect.loc[unique_mutect, :]
    
    # merge
    all_var = pd.concat([mutect, varscan], ignore_index=True)
    
    # write to file
    all_var.to_csv(inoutpath + name + '_all_somatic_annotated.tsv', sep='\t', header=True, index=False)
  

    # delete temporary annovar files
    sp.call(' '.join(['rm', inoutpath + name + '*_annovar_*.tsv']), shell=True)


    
    
def filterExonicPolymorphic(name):
    
    import pandas as pd
    import subprocess as sp
    
    inoutpath = '10_variants_filter/'
    logpath = '10_variants_filter/logs/'
    mainpath = '/home/PERSONALE/eugenio.fonzi2/'

    
    # load tables
    variants = pd.read_table(inoutpath + name + '_all_somatic_annotated.tsv', sep='\t', header=0)

    # keep exonic
    variants = variants[variants.exon == 'exonic']
    
    # remove synonymous
    variants = variants[variants.var_type != 'synonymous SNV']

    # remove polymorphisms according to 1000g, ESP and EXAC
    remove = variants[(variants['1000g'] >= 0.01) | (variants.esp >= 0.01) | (variants.exac >= 0.01)].index
    variants = variants.drop(remove)

    # select rsIDs and write to file
    rsID = variants.dbsnp_id[pd.notnull(variants.dbsnp_id)]
    rsID.to_csv('rsID.tsv', sep='\t', header=False, index=False)
    
    # run R script that retrieves MAF values from rsIDs
    cmd = ' '.join(['Rscript', mainpath + 'script/rsIDfilter.R',
                    '2>', logpath + name + '_rsID_err.log'])
    
    sp.call(cmd, shell=True)
                
    # load output of R script
    rsID = pd.read_table('rsID_maf.tsv', sep='\t', header=0)
    
    # if for any rsID a MAF value was retrieved
    if len(rsID) > 0:
        
        # select rsIDs with MAF >= 0.01
        rsID_remove = rsID[rsID.MAF >= 0.01]['Query']
        
        # get indeces of those rsIDs in variants table
        remove = [i for i in variants.index if variants.loc[i, 'dbsnp_id'] not in rsID_remove]
        
        # remove them
        variants = variants.drop(remove)
    
    # write to file
    variants.to_csv(inoutpath + name + '_all_somatic_annotated_filtered.tsv', sep='\t', header=True, index=False)
    
    sp.call('rm rsID*', shell=True)
    
    
    
    
def deleteIntermediate(name, step=None):
    
    import subprocess as sp
    
    if step == 'alignment':
    
        # save tree structure as dictionary
        dirs = {'01_fastq/': '',
                '02_fastq_trimmed/': '',
                '03_alignment_genome/': '01_intermediate/',
                '04_alignment_exome/': '01_intermediate/', 
                '05_alignment_MT/': '01_intermediate/', 
               }
    
    
    if step == 'variant_calling':
        
        # save tree structure as dictionary
        dirs = {'08_varscan_genome/': '01_mpileup/',
                 '09_varscan_MT/': '01_mpileup/', 
                }
    
    # loop over dictionary keys (parent directories)
    for d in dirs:
        
        # delete all content in each directory
        sp.call(' '.join(['rm', d + dirs[d] + '*']), shell=True)
 

