### save session info to file
### to be used into bash scripts


import IPython
import subprocess as sp

info = IPython.sys_info().splitlines()

for l in info:
    
    with open('session_info', 'a') as outfile:
        outfile.write(l + '\n')

     
info = sp.run("pip3 freeze", shell=True, stdout=sp.PIPE).stdout.decode().splitlines()

for l in info:
    
    with open('session_info', 'a') as outfile:
        outfile.write(l + '\n')
