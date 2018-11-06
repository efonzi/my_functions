#####################################################
######################## ATOM #######################
#####################################################
## apm is Atom Package Manager, the best way to install packages

# AOSP proxy stops the installation of packages...
# the command 'apm config' allows to configure apm to use AOSP proxy
# 'http://192.168.65.4:8080' is the address of AOSP proxy
apm config set https-proxy http://192.168.65.4:8080



###################################################################
######################## PEDAGOGIA - ped432 #######################
###################################################################
## username: martinelli
## password: antonella
## ip: 172.20.45.196

# access from local
lftp martinelli@ped432
lftp martinelli@172.20.45.196
ftp 172.20.45.196 # then type username and password

# download files to local
lftp ftp://martinelli@172.20.45.196 -u martinelli,antonella -e "mirror --verbose SOURCE DEST; bye"
# SOURCE must be a folder, doesn't work with a single file --> ex. '20180521_Martinelli/MartinelliMetaboloma/'
# DEST --> ex. '/media/DATA1/ped/'
## DOESN'T WORK FROM NAS!!!!!!





###############################################################
####################### CONDA/JUPYTER #########################
###############################################################
# username: efonzi
# password: Anaconda1a.
####### detailed description of PIP, CONDA, JUPYTER
# http://jakevdp.github.io/blog/2017/12/05/installing-python-packages-from-jupyter/


#################  USEFUL COMMANDS ############
# https://conda.io/docs/user-guide/tasks/manage-environments.html#
conda create -n my_env
conda install -n my_env package_name # ex 'conda install -n my_env r=3.4.3'

# list available version of a package
conda search package_name

# list environments (1)
conda info --envs
# list environments (2)
conda env list

# list packages in an environment
conda list -n myenv

# remove an environment
conda remove --name myenv --all

###### add/append every channel
# add sets highest priority, while append sets lowest priority
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels statsmodels
conda config --append channels dranew


###### TO CREATE CONDA ENVIRONMENT 'snake' ######
conda create -y -n snake python=3.6.5 snakemake=5.1.5 ## until 181029
conda create -y -n snake python=3.6.5 snakemake=5.1.5 perl=5.26.2 perl-lwp-protocol-https ## since 181029, for snakefile_call_02.py

###### TO CREATE CONDA ENVIRONMENT 'py365_euge' ######
conda create -y -n py365_euge python=3.6.5 ipykernel jupyter pandas numpy matplotlib seaborn statsmodels matplotlib-venn pyyaml biopython rpy2 xlrd snakemake lxml samtools

source activate py365_euge
python -m ipykernel install --user --name py365_euge --display-name 'py365_euge'
source deactivate


###### TO BUILD 'r-rawcopy' PACKAGE FROM SCRATCH #######
conda install anaconda-client -y
anaconda login # then type username and password of ANACONDA account
conda install conda-build -y
##
mkdir ~/my_functions/r-rawcopy
> ~/my_functions/r-rawcopy/meta.yaml
> ~/my_functions/r-rawcopy/build.sh
> ~/my_functions/r-rawcopy/bld.bat
# fill these 3 files as in the ones saved in '~/my_functions/r-rawcopy/'
##
conda-build ~/my_functions/r-rawcopy
conda convert --platform all ~/miniconda3/conda-bld/linux-64/r-rawcopy-1.1-r3.3.2_0.tar.bz2 -o ~/miniconda3/conda-bld/
anaconda upload ~/miniconda3/conda-bld/linux-64/r-rawcopy-1.1-r3.3.2_0.tar.bz2
anaconda upload ~/miniconda3/conda-bld/linux-32/r-rawcopy-1.1-r3.3.2_0.tar.bz2
anaconda upload ~/miniconda3/conda-bld/osx-64/r-rawcopy-1.1-r3.3.2_0.tar.bz2
anaconda upload ~/miniconda3/conda-bld/linux-aarch64/r-rawcopy-1.1-r3.3.2_0.tar.bz2
anaconda upload ~/miniconda3/conda-bld/linux-armv6l/r-rawcopy-1.1-r3.3.2_0.tar.bz2
anaconda upload ~/miniconda3/conda-bld/linux-armv7l/r-rawcopy-1.1-r3.3.2_0.tar.bz2
anaconda upload ~/miniconda3/conda-bld/linux-ppc64le/r-rawcopy-1.1-r3.3.2_0.tar.bz2
#anaconda upload ~/miniconda3/conda-bld/win-32/r-rawcopy-1.1-r3.3.2_0.tar.bz2
#anaconda upload ~/miniconda3/conda-bld/win-64/r-rawcopy-1.1-r3.3.2_0.tar.bz2
# windows conversion didn't work


###### TO CREATE CONDA ENVIRONMENT 'r332_rawcopy' ######
conda create -y -n r332_rawcopy
conda install -y -n r332_rawcopy -c efonzi r-rawcopy
##
python3 ~/my_functions/conda_edit_pin_file.py -n r332_rawcopy



###### TO CREATE CONDA ENVIRONMENT 'py365_rsID' ######
# to also run ~/WESPipeline/scriptsrsIDquery.01.R
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels statsmodels
conda config --add channels efonzi
##
conda create -y -n py365_rsID python=3.6.5 ipykernel
conda install -y -n py365_rsID r-rsnps
conda install -y -n py365_rsID pandas numpy matplotlib seaborn statsmodels matplotlib-venn
source activate py365_rsID
python -m ipykernel install --user --name py365_rsID --display-name 'py365_rsID'
##
python3 ~/my_functions/conda_edit_pin_file.py -n py365_rsID
##
## alternatively, but creates an awful output everytime the env is de/activated
# conda install -n rsID r-plyr
# conda install -n rsID r-stringr
# conda install -n rsID r-httr
# conda install -n rsID r-xml



# to check the directories of the kernels currently installed
jupyter kernelspec list

jupyter notebook --no-browser --port=8231





####### CADAVER SCRIPT!!!
# automatically create script for cadaver and run it
# ----> https://gist.github.com/fheinle/8b78b0e83f0c8f4bd9da


####### if __name__ == '__main__':


####### WRITE AND LOAD YAML
import yaml
dictionary = {'one': 'a', 'two': 'b'}
with open('file.yaml', 'w') as outfile:
    yaml.dump(dictionary, outfile)
with open('file.yaml') as infile:
    dictionary = yaml.load(infile)


####### PATHWAY NAME MANIPULATION
import os
# to get parent directory of pwd
os.path.abspath(os.path.join(os.getcwd(), os.pardir))
# get common prefix of 2 strings (paths)
os.path.commonprefix([string1, string2])
# exclude directories
files = [f for f in files if os.path.isfile(path + f)]


####### SESSION INFO
import IPython
import subprocess as sp
# print 'session info'
print(IPython.sys_info())
# print all installed modules with their version
sp.Popen("pip3 freeze", shell=True, stdout=sp.PIPE).stdout.read().decode().splitlines()
# OR in jupyter
get_ipython().system('pip3 freeze')
# OR in jupyter
!pip3 freeze


####### FOUR EXAMPLES OF COLOR PALETTES IN SEABORN
# set notebook to show multiple outputs
get_ipython().run_line_magic('matplotlib', 'inline')
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
# print palettes
import seaborn as sns
sns.palplot(sns.color_palette('hls', 20))
sns.palplot(sns.color_palette('husl', 20))
sns.palplot(sns.color_palette('Paired', 20))
sns.palplot(sns.color_palette('Set2', 20))


####### check if a table file exists before loading it
# if it doesn't exist, create an empty pandas DF
try:
    table = pd.read_table('file_name', sep='\t', header=0)
except:
    table = pd.DataFrame()


####### TO RUN A LOOP OVER BATCHES OF A DETERMINED SIZE AT A TIME, INSTEAD OF THE ENTIRE OBJECT
# n = total number of elements (ex. 2389)
# d = desired number of elements in a batch (ex. 100)
# compute quotient and remainder of n/d
q, r = n//d, n%d
# loop over batches
for batch in range(q + 1):
    # determine size of batch (all of the same size, aside from the last one)
    if batch < q:
        size = d
    else:
        size = r
    # loop over batch indeces
    for index in range(size):
        # transform batch index to total index
        INDEX = batch * size + index

        # run some code
        # .............


####### write a list-like object to file
# each elements will be a different line in the file
file = open('file_name', 'w')
file.writelines( "%s\n" % e for e in objectToBeWritten )
file.close()


####### capture STDOUT of a bash command
sp.Popen("bash command", shell=True, stdout=sp.PIPE).stdout.read().decode().splitlines()



####### QUERY NCBI FOR CYTOBAND OF GENE
from Bio import Entrez
Entrez.email = 'eugenio.fonzi2@unibo.it'
geneSymbol = 'MCL1'
# build query string
query = "(" + geneSymbol + "[gene]) AND (Homo sapiens[orgn])"
# send request to Entrez
handle = Entrez.esearch(db="gene" ,term=query, tettype="gene", retmode="xml")
# de-serialize output XML
record = Entrez.read(handle)
# close handle
handle.close
# extract entrez ID of gene (ONLY THE FIRST IN THE LIST!!!!!!)
ID = record["IdList"][0]
# send request to Entrez
handle = Entrez.efetch(db="gene" ,id=ID, tettype="gene", retmode="xml")
# de-serialize output XML
record2 = Entrez.read(handle)
# close handle
handle.close
# extract cytoband were gene is located
cytoband = record2[0]['Entrezgene_location'][0]['Maps_display-str']


############################################
############### PANDAS
# returns TRUE if the element of an array/DF is 'NaN'
pd.isnull(df.loc['x', 'y'])
# returns TRUE if the element of an array/DF is NOT 'NaN'
pd.notnull(df.loc['x', 'y'])
# to check if the elements of a DF or an array-like are NaN
df.isnull()
# to check if each column of DF has any NaN
df.isnull().any()
# to check if a DF has any NaN
df.isnull().any().any()
# when slicing a DF assigning it to the same name, do like
df = df['columnX']
# when slicing a DF assigning it to a different name, do like
df2 = df['columnX'].copy()
# otherwise it will throw a warning when trying to slice or modify the new DF
# (see https://stackoverflow.com/questions/38147027/action-with-pandas-settingwithcopywarning/38147527#38147527)
# to merge two list-like objects without replicates
obj1.union(obj2)
# to check if two pandas DF are the same or not
df1.equals(df2)

# to slice a MultiIndex ('https://www.somebits.com/~nelson/pandas-multiindex-slice-demo.html')
# this will slice for the 2nd and 4th LEVEL of 4-level MultiIndex
df.loc[pd.IndexSlice[:, 'labelX', :, 'labelY'], :]


##################################################################
####### configure jupyter to print multiple output after the cell
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
# TO DO IT PERMANENTLY
# (see https://www.dataquest.io/blog/jupyter-notebook-tips-tricks-shortcuts/)
# create the file '~/.ipython/profile_default/ipython_config.py'
# and write in it the following:
c = get_config()
# Run all nodes interactively
c.InteractiveShell.ast_node_interactivity = "all"


####### increase size of plots (default is [6.0, 4.0])
plt.rcParams['figure.figsize'] = [12.0, 8.0]


####### Return the subset of the list of names that match pattern. It is the same as [n for n in names if fnmatch(n, pattern)], but implemented more efficiently.
fnmatch.filter(names, pattern)


####### How to capture the output of Python's help() function
import io
import sys
# Temporarily redirect stdout to a StringIO.
stdout = sys.stdout
s = io.StringIO()
sys.stdout = s
help(sys.intern) # <--- here put what you are looking "help()" for, as a string (ex. help('pandas'))
# Don't forget to reset stdout!
sys.stdout = stdout
# Read the StringIO for the help message.
s.seek(0)
help_string = s.read()
print(help_string)


####### start jupyter notebook
jupyter notebook --no-browser --port=8231

####### to install package as user (from shell)
$ python3.4 -m pip install --user packageName
