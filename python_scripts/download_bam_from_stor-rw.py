# 181009

'''
For each sample in a list of tumor_normal pairs, search which remote directory
contains its BAM and download it to current local directory. It uses 'sshpass'
to provide password in command line (security unwise)
'''

#  sshpass -p "password" rsync -av /SOURCE /DEST
# Note the space at the start of the command; in the bash shell
# this will stop the command (and the password) from being stored in the history.


import os
import subprocess as sp

homepath = os.path.expanduser('~') + '/'

bash_command = 'sshpass -p "Seragn0la8a."'
remote = 'eugenio.fonzi2@137.204.48.206:'
storepath = '/mnt/stor-rw/users/eugenio_wes/'

runs = [r + '_align/' for r in ['180707', '180801', '180810', '180910']]
subfolder = {'nuclear/': '03_alignment_genome/02_bqsr/',
             'MT/': '05_alignment_MT/02_bqsr/'}

for sf in subfolder:
    os.makedirs(sf, exist_ok=True)

path = 'pair_list.tsv'
with open(path) as infile:
    pair_list = infile.read().splitlines()


for p in pair_list:

    for s in p.split('_'):

        for r in runs:

            for sf in subfolder:

                source = '%s%s%s%s%s.ba*' %(remote, storepath, r, subfolder[sf], s)

                # Note the space at the start of the command; in the bash shell
                # this will stop the command (and the password) from being stored in the history.
                cmd = ' %s rsync -av %s %s' %(bash_command, source, sf)

                try:
                    sp.run(cmd, shell=True)
                except:
                    continue

        continue
