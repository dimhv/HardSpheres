#!/usr/bin/python3
#L.H.Miranda-Filho., ITP-CAS, lucmiranda@gmail.com, 2023

import numpy as np
import time
import subprocess




#sample_index = 3
##################################################################

processes = []
ll = 0


samples = 88

for jj in range(0, samples):

    cmd = './3d_SG_parallel '+' '+str(jj+40)
    print (cmd)

    processes.append('p'+ str(ll))
    processes[ll] = subprocess.Popen(cmd, shell=True)

    cmd = '/bin/sleep 2'
    print (cmd)
    print ()
    subprocess.call(cmd, shell=True)

    ll += 1

    if( ll%1==0 ):
        for mm in processes:
            mm.wait()

        del processes
        processes = []
        ll = 0

##################################################################

 
