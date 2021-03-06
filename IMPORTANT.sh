#############################################################
############################## R ############################
#############################################################

###### INSTALL CUSTOM R PACKAGES FROM COMMAND LINE ########

# Create a custom path for every R package (mine will be ~/pkgs_R_euge/lib)
mkdir ~/pkgs_R_euge
mkdir ~/pkgs_R_euge/lib
echo "R_LIBS=~/pkgs_R_euge/lib" > ~/.Renviron

# from inside R, I can check the list of paths where R will look for libraries with
> .libPaths()

# install the package in the desired path (need to download a tarball first)
$ R CMD INSTALL -l ~/pkgs_R_euge/lib ~/path/to/package.tar.gz

## CAVEAT

# Packages will be installed for the version of R of the current active environment, so it won't work in environments with different R versions!!

# installing a package will overwrite any previous installations of that package, regardless of what version of the package they were

# DO NOT INCLUDE A '/' AT THE END OF THE PATH!!!
# '~/custom/path' --> OK
# '~/custom/path/' --> NO!!




######## R su linux MARCO ########
###170523
###version 3.3.2 is installed in usr/local/bin/ and it's the one to which the terminal and Rstudio connect automatically
###then, we installed version 3.3.3 in home/eugenio/ with the following command
./configure --enable-R-shlib=yes --with-readline=no --with-x=no
###these options are needed because the proxy "blocks" our installation
###as there is a problem when Rstudio tries to contact the version 3.3.2,
###to make Rstudio run we must tell it to look for R3.3.3 in its directory instead of R3.3.2
###to do this we added the following chunk
export RSTUDIO_WHICH_R=~/R-3.3.3/bin/R
###to the end of the file /etc/bash.bashrc
###/etc/bash.bashrc is a script that is run by ubuntu every time it opens a terminal
###################################

anyDuplicated() --> returns the index of the duplicated elements

when subsetting like
> x = matrix(1:9, ncol=3, nrow=3)
> x[1,]
an atomic vector is returned, instead of a 2-dimensional object. To avoid this use "drop" argument
---> x[1,, drop=FALSE] 


library(methods)
as(x, "matrix") ---> convert object "x" to matrix

##IRanges functions
start()
end()
width()
length() #because they are vectors
c() # to concatenate
reduce() # to convert an IRange to its normal representation (merges overlapping ranges)
disjoin() #quite the opposite of reduce(), it creates non-overlapping intervals
resize() # makes the intervals of the same chosen size, starting from the respective start or ending in the respective end. Or starting from the "center"

union() #outputs a "normal IRanges"
# union(ir1, ir2) is equal to reduce(c(ir1,ir2))

intersect()

findOverlaps(query, subject) # also queryHits(), subjectHits()

countOverlaps(query, subject) # faster and more memory efficient than findOverlaps()

nearest(query, subject) # which element in subject is the closest to each element in query

## GRanges specific functions ---> because they have strand
promoters() # default --> 2000 bp upstream up to 200 downstream from start

seqinfo()
seqlengths()

seqlevels()
seqnames() # can't add new name if it's not in seqlevels (like with factors)

genome() # error if trying to make operations on ranges labelled with different genomes

gaps() # gives all parts of a chr that are not covered by ranges (appears whole "*" strand)

sort() # the results is dependent on the order of the levels given by seqlevels()

DataFrame() # to make a special dataframe that can handle IRanges data

values() # for metadata columns, which are a separated DataFrame inside the GRange

findOverlaps() # "*" strand overlaps with both "+" and "-"

subsetByOverlaps(query, subject)

makeGRangesFromDataFrame(df) # "keep.extra.columns" argument to keep metadata

dropSeqLevels()
keepSeqLevels()
keepStandardChromosomes() 

> newStyle = mapSeqlevels(seqlevels(gr), "NCBI") # creates vector of convertion between our seqlevels and NCBI chromosome names
> gr = renameSeqlevels(gr, newStyle)

 
## Rle functions
Rle() # converts a numeric vector to an Rle

coverage() # converts an IRanges to a Rle

############################## SERVER FISICA ###############################
############################################################################

#----------------- CONFIGURE ACCESS TO SERVER BIO4 ------------------------

### open config file for ssh
nano ~/.ssh/config

### type following info in file
Host    bio4
        HostName 137.204.48.206
        User eugenio.fonzi2
        LocalForward 4444 localhost:8666

### save and close file

#--------------------- ACCESS SERVER  BIO4 ------------------------

### type
ssh bio4

### insert password UNIBO

#-------------------------- SCREEN -------------------------
### create a "screen"
screen -S screenName

### to leave a screen without halting any running pr ocess
# push "CTRL+A", then push "D"

### to go back to the screen
screen -r screenName

### to kill a screen
screen -X -S screenName quit

### to navigate inside one's own "screen"
screen -r ALL
screen -r ipython
screen -ls





#---------------------- JUPYTER --------------------

### 1. from remote, type
ipython notebook --no-browser

### 2. from local, type
ssh -N -L localhost:8888:localhost:8889 eugenio.fonzi2@137.204.48.206

### 3. from browser, access "localhost:8888"


# usually, ports 8888 and 8889 are automatically decided, but you can choose it manually with:
ipython notebook --no-browser --port=xxxx
ssh -N -L localhost:xxxx:localhost:xxxx eugenio.fonzi2@137.204.48.206

# if xxxx is busy, ipython will choose a yyyy port, then:
ssh -N -L localhost:xxxx:localhost:yyyy eugenio.fonzi2@137.204.48.206

# finally, access "localhost:xxxx" from browser

# detailed info at this page ---> https://coderwall.com/p/ohk6cg/remote-access-to-ipython-notebooks-via-ssh





# -------------------- PYTHON ---------------------

#### to install a python package only for one user of the server
cd ~/
python3.4 -m pip install --user packageName

# (adapt command 'pythonX.Y' to current version of Python)





# -------------------- RSTUDIO ---------------------

#### RStudio from server
- access bio4 in terminal
- open browser
- navigate to "localhost:4444"
- sign in with unibo username and password





#------------------ RSYNC ------------------------

### to syncronize folders/files between local and our server ###
$ rsync -av --delete SOURCE DEST

other options
-n  dry run
--exclude 'NAS'  will exclude the folder 'NAS' from the sync

remote server is specified like:
username@ip_address:/directory/path/
example:
eugenio.fonzi2@137.204.48.206:~/wes_analyses/bam/

# TO SPECIFY THE PASSWORD ON COMMAND LINE
# Note the space at the start of the command; in the bash shell
# this will stop the command (and the password) from being stored in the history.
$  sshpass -p "password" rsync -av /SOURCE /DEST



#----------------------------------------------------------------------------------------

### in the following folder we have many software tools
ls /home/condivisi/softwares

abacas         CONTIGuator_v2.7                             htslib-1.3         ncbi-igblast-1.6.0   SOAPdenovo2                      velvet-master
annovar        cufflinks-2.2.1.Linux_x86_64                 MaSuRCA-3.1.3      NGSQCToolkit_v2.3.3  sratoolkit.2.7.0-centos_linux64  VelvetOptimiser-master
bamtools       download.pl?toolkit=NGSQCToolkit_v2.3.3.zip  mauve              ngsutils             STAR
bedtools2      FastQC                                       MUMmer3.23         picard-tools-2.3.0   tabix-0.2.6
bigWigToWig    fastx-toolkit                                muscle             prokka               telegraf_1.1.0_amd64.deb
bowtie2-2.3.0  fossil                                       mutect             RSEM-1.3.0           Trimmomatic-0.36
bwa-0.7.13     gatk                                         ncbi-blast-2.3.0+  samtools-1.3         vcftools_0.1.13


### to call them I have to do like
/home/condivisi/softwares/mutect

### but for some of them, just the name will do, like
bwa

### others are like
java -jar $GATK
java -jar $TRIMMOMATIC
java -jar $PICARD

#----------------------------------------------------------------------------------------

####################################################################
############################## NAS #################################
####################################################################

### my account
user: efonzi
pass: Seragn0la1a.

### account Anto
user: apadella
pass: a(B32$

### account Marco
user: mmanfrini
pass: r86HYU0DIZ (could be lowercase)

################## WebDav ACCESS ##################
http --> 5005
https --> 5006

### from outside hospital network #####
https://ngs-ptl.unibo.it:5006
# it works from MacOSX's FINDER
# it doesn't work from ubuntu's Nautilus (it doesn't accept "https", only "dav(s)")
# it works from Cyberduck with https://ngs-ptl.unibo.it:5006 (as WebDAV)
# it works from ubuntu with command "sudo mount -t davfs https://ngs-ptl.unibo.it:5006 ~/path/to/folder", after installing "ca-certificates" with "sudo apt install ca-certificates"
# it works from Mac OS X with command "mount_webdav -i https://ngs-ptl.unibo.it:5006 ~/path/to/folder"
# to unmount "umount ~/path/to/folder"

### other mixed webdav info
ngs-ptl.personale.dir.unibo.it:5005
nasematologia.aosp.bo.it:5005
cadaver nasematologia.aosp.bo.it:5005  #oppure :5006
personale\str00964-ngs-ptl
personale\antonella.padella2
personale\eugenio.fonzi2
personale\giovanni.martinelli2


############### CREARE CHIAVE DI ACCESSO RSA PER SERVER ################
# a public key (shorter) is saved on the server and a private key (longer) is saved on the local
> ssh-keygen
> ssh-copy-id serverName # like “bio4”

# the public key is stored in “~/.ssh/id_rsa.pub” and the private key in “~/.ssh/id_rsa”

#=====================================
#=============== ECDSA ===============
#=====================================

## on October 2018 a cryptocurrency miner hacked into BIO4 using my UNIBO account, so I had to change password and SSH KEYS

# delete old keys
rm .ssh/id_rsa*
rm .ssh/known_hosts

# generate new keys with ECDSA algorithm
ssh-keygen -t ecdsa -b 521

# copy public key into BIO1
ssh-copy-id -i ~/.ssh/id_ecdsa eugenio.fonzi2@137.204.48.143

# copy public key into BIO3
ssh-copy-id -i ~/.ssh/id_ecdsa eugenio.fonzi2@137.204.48.205




#######################################################################
############################## OTHERS #################################
#######################################################################



### important job skills
Proficiency with Python and/or R (NumPy, SciPy, Pandas, Bioconductor and other data analysis technologies would also be beneficial)
