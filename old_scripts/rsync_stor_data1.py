import os
import subprocess as sp
import json

with open('/media/DATA1/dic.json') as infile:
    dic = json.load(infile)


path = '/media/DATA1/stor-rw_backup/'

mt_bam = os.listdir(path + 'mt_bam/')
nuclear_bam = os.listdir(path + 'nuclear_bam/')
trim = os.listdir(path + 'trim_logs/')

for b in trim:
    print(b, dic['nuclear'][b])
    directory = dic['nuclear'][b]
    os.makedirs(directory, exist_ok=True)
    sp.run('mv %strim_logs/%s %s' %(path, b, directory), shell=True)

for b in mt_bam:
    print(b, dic['mt'][b])
    directory = dic['mt'][b]
    os.makedirs(directory, exist_ok=True)
    sp.run('mv %smt_bam/%s %s' %(path, b, directory), shell=True)

for b in nuclear_bam:
    print(b, dic['nuclear'][b])
    directory = dic['nuclear'][b]
    os.makedirs(directory, exist_ok=True)
    sp.run('mv %snuclear_bam/%s %s' %(path, b, directory), shell=True)
