
import os
import pandas as pd
import re
from zipfile import ZipFile
import matplotlib.pyplot as plt

z = '/media/DATA1/122bis_R1_fastqc.zip'

name = z.split('/')[-1].replace('_fastqc.zip', '')

# open ZIP archive and open file with FASTQC result data
# split data between all the different analyses of FASTQC
# (each part ends with the string '>>END_MODULE')
with ZipFile(z) as myzip:
    with myzip.open(name + '_fastqc/fastqc_data.txt') as myfile:
        fastqc = myfile.read().decode().split('>>END_MODULE')

fastqc

pbsq = [e for e in fastqc if re.search('>>Per base sequence quality', e)][0].splitlines()[2:]

pbsq = [e.split('\t') for e in pbsq]

pbsq = pd.DataFrame(pbsq[1:], columns=pbsq[0])

pd.Series(data=pbsq['Lower Quartile'].tolist(), index=pbsq['#Base'])


values = pbsq['Lower Quartile'].astype('float').tolist()

%matplotlib inline

fig = plt.figure()
ax = plt.axes()
plt.plot(range(len(pbsq)), values)
plt.ylim([0,40]);


#plt.show()
