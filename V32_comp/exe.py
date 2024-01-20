#!/usr/bin/python3
#L.H.Miranda-Filho., ITP-CAS, lucmiranda@gmail.com, 2023

import numpy as np
import time
import subprocess



size_bac = 6

#sample_index = 3
##################################################################

processes = []
ll = 0

#for ii in range(19,24):

    #cmd = 'mkdir output_'+ str(ii).zfill(2)
    #print (cmd)
    #subprocess.call(cmd, shell=True)

    #cmd = '/bin/sleep 1'
    #print (cmd)
    #subprocess.call(cmd, shell=True)

    #cmd = 'mkdir ./output_'+str(ii).zfill(2)+'/conf'
    #print (cmd)
    #subprocess.call(cmd, shell=True)

    #cmd = '/bin/sleep 1'
    #print (cmd)
    #subprocess.call(cmd, shell=True)

ii = 5

cmd = 'mkdir ./output/output_'+ str(ii).zfill(2)
print (cmd)
subprocess.call(cmd, shell=True)

cmd = 'mkdir ./output/output_'+str(ii).zfill(2)+'/conf'
print (cmd)
subprocess.call(cmd, shell=True)

for jj in range(0,size_bac):

    cmd = './3d_hs_seq_main '+' '+str(jj)+' '+str(ii)
    print (cmd)

    processes.append('p'+ str(ll))
    processes[ll] = subprocess.Popen(cmd, shell=True)

    cmd = '/bin/sleep 2'
    print (cmd)
    print ()
    subprocess.call(cmd, shell=True)

    ll += 1

    if( ll%24==0 ):
        for mm in processes:
            mm.wait()

        del processes
        processes = []
        ll = 0

##################################################################

 
