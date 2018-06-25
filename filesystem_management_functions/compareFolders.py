
'''
Given two directories as input (folder1 and folder2), the script will compare
the files they contain. It will:
(1) print to STDOUT the list of the files unique to folder1
(2) print to STDOUT the list of the files unique to folder2
(3) run MD5 hashing on the common files (meaning the files with the same name)
(4) delete from folder1 the common files that also had identical content (same md5)
(5) print to STDOUT the list of the commong files that have different content

RUN LIKE:
$ python3 compareFolders.py -1 path/to/folder1 -2 path/to/folder2

The paths can be absolute or relative to CWD
'''

import os
import subprocess as sp
import argparse


def editFolderPath(f):
    '''
    function to derive the absolute directory path
    '''
    # get homepath
    homepath = os.path.expanduser('~') + '/'
    # get current path
    currentpath = os.getcwd() + '/'
    # in case a relative path was given as input
    if not os.path.commonprefix([homepath, f]):
        # append to current path
        f = currentpath + f
    # add '/' as suffix, if missing
    if f[-1:] != '/':
        f = f + '/'
    return(f)



#### Parse and edit command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-1', required=True, dest='folder1', action='store', help='path to folder1')
parser.add_argument('-2', required=True, dest='folder2', action='store', help='path to folder2')

args = parser.parse_args()

# derive absolute directory paths
folder1 = editFolderPath(args.folder1)
folder2 = editFolderPath(args.folder2)


# get list of files in folders
files1 = os.listdir(folder1)
files2 = os.listdir(folder2)

# print files unique to folder1
print('\nunique to folder1'.upper())
print('\n'.join([f for f in files1 if f not in files2]))

# print files unique to folder2
print('\nunique to folder2'.upper())
print('\n'.join([f for f in files2 if f not in files1]))


different = []

# loop over files in folder1
for f in files1:

    # if the file is also contained in folder2
    if f in files2:

        # run MD5 hashing on the file in folder1
        md51 = sp.run('md5sum %s%s' %(folder1, f), shell=True, stdout=sp.PIPE).stdout.decode().split('  ')[0]

        # run MD5 hashing on the file in folder2
        md52 = sp.run('md5sum %s%s' %(folder2, f), shell=True, stdout=sp.PIPE).stdout.decode().split('  ')[0]

        # if the hash codes are different, append file name to list of different files
        if md51 != md52:
            different.append(f)
        # else, delete file in folder1
        else:
            sp.run('rm %s%s' %(folder1, f), shell=True)


# print list of common files with different hashed codes
print('\nsame name but different'.upper())
print('\n'.join(different))


##### END #####
