"""
This script produces a file called 'pinned' which purpose is to prevent CONDA from automatically update the packages installed
in an environment. The file 'pinned' specifies the version that each package must be kept at. To create the file 'pinned' it 
is necessary to edit the output of the 'conda list -n my_env'
"""

import os
import argparse
import re
import subprocess as sp

homepath = os.path.expanduser('~') + '/'

# parse 'my_env' name from command line
parser = argparse.ArgumentParser()
parser.add_argument('-n', dest='my_env', required=True)
args = parser.parse_args()
my_env = args.my_env

# list packages installed in 'my_env' and capture STDOUT
cmd = 'conda list -n %s' %my_env
pins = sp.run(cmd, shell=True, stdout=sp.PIPE).stdout.decode().splitlines()[2:]

# in each line substitute ' ' with ' ==' (ex. 'packageX 5.6' becomes 'package ==5.6')
pins = '\n'.join([' =='.join(re.split('\s+', p)[:2]) for p in pins])

# write to file
dest = '%sminiconda3/envs/%s/conda-meta/pinned' %(homepath, my_env)
with open(dest, 'w') as outfile:
    outfile.write(pins)