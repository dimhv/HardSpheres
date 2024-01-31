#!/usr/bin/python3
#L.H.Miranda-Filho., ITP-CAS, lucmiranda@gmail.com, 2023

import numpy as np
import time
import subprocess



#vec = [2,3,4,5,6,7,8,9,10,11,12,13,14,17,18,19,21]
#sample_index = 3
##################################################################

processes = []
ll = 0


ii = 2

for jj in range(0, 6):

    cmd = './3d_hs_parallel_main '+' '+str(jj)+' '+str(ii)
    print (cmd)

    processes.append('p'+ str(ll))
    processes[ll] = subprocess.Popen(cmd, shell=True)

    cmd = '/bin/sleep 2'
    print (cmd)
    print ()
    subprocess.call(cmd, shell=True)

    ll += 1

    if( ll%3==0 ):
        for mm in processes:
            mm.wait()

        del processes
        processes = []
        ll = 0

##################################################################

 
