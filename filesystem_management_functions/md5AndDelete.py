
'''
Walk down the subdirectories of a chosen directory 'path' and HASH every file
in each subfolder with MD5SUM. Then delete every file aside from the newly created
.md5 and a list of chosen files.
'''

import subprocess as sp
import os

# set path
path = '/home/PERSONALE/eugenio.fonzi2/trial/mitochondrial_genes/'

# list of file names that must not be deleted after MD5SUM
to_save = ['every_file.md5', 'README.sh', 'stdout_stderr.txt', 'stderr.txt']

# os path walk
for root, dirs, files in os.walk(p):

    # if there are files in current root
    if files:
        # run MD5SUM of the files in root
        cmd = 'md5sum %s/* > %s/every_file.md5' %(root, root)
        sp.run(cmd, shell=True)

        # loop over files in root
        for file in files:

            # avoid files in list 'to_save'
            if file not in to_save:

                # build path to file
                f = root + '/' + file

                # remove file if it is not a symbolik link
                if not os.path.islink(f):
                    os.remove(f)
