
import os
import pandas as pd

import matplotlib.pyplot as plt
from zipfile import ZipFile
%matplotlib inline

#path = '/Users/Eugenio/Documents/PRE/'
path = '/home/eugenio/Scrivania/PRE/'

#z = path + '122bis_R1_fastqc.zip'

zlist = [path + z for z in os.listdir(path) if z[-4:] == '.zip']

dic = {}

for z in zlist:

    name = z.split('/')[-1].replace('_fastqc.zip', '')

    # open ZIP archive and open file with FASTQC result data
    # split data between all the different analyses of FASTQC
    # (each part ends with the string '>>END_MODULE')
    with ZipFile(z) as myzip:
        with myzip.open(name + '_fastqc/fastqc_data.txt') as myfile:
            fastqc = myfile.read().decode().split('>>END_MODULE')


    pbsq = [e for e in fastqc if '>>Per base sequence quality' in e][0].splitlines()[2:]
    pbsq = [e.split('\t') for e in pbsq]
    pbsq = pd.DataFrame(pbsq[1:], columns=pbsq[0])
    # if name == '386_R1':
    #     print(pbsq)
    perc10th = []

    for i, row in pbsq.iterrows():
        q = float(row['10th Percentile'])
        b = row['#Base']
        b = [int(e) for e in b.split('-')]
        b = len(range(b[0], b[-1])) + 1
        perc10th = perc10th + ([q] * b)

    dic[name] = perc10th

    plt.plot(range(len(perc10th)), perc10th)
    plt.ylim([0,40])
