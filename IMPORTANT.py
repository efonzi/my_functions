

####### detailed description of PIP, CONDA, JUPYTER
# http://jakevdp.github.io/blog/2017/12/05/installing-python-packages-from-jupyter/


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
# extract entrez ID of gene
ID = record["IdList"][0]
# send request to Entrez
handle = Entrez.efetch(db="gene" ,id=ID, tettype="gene", retmode="xml")
# de-serialize output XML
record2 = Entrez.read(handle)
# close handle
handle.close
# extract cytoband were gene is located
cytoband = record2[0]['Entrezgene_location'][0]['Maps_display-str']


####### PANDAS
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
jupyter notebook --no-browser --port=8761/2/3/4


####### to install package as user (from shell)
$ python3.4 -m pip install --user packageName
