#!/usr/bin/python3

import numpy as np
import sys

import concurrent.futures




#samples = 64
#delta_t = 100
#tobs = 20000

samples = int(sys.argv[1])
delta_t = int(sys.argv[2]) #delta_t consists the frame number in the main code
tobs = int(sys.argv[3])


time = []
jj = 0
while(jj < tobs):
    time.append(jj)
    jj += delta_t

print()
print("time:")
print(time)
print()

#chi = np.zeros(len(time))

M = int(0.5*samples*(samples-1))

#####################################################
def process_item(kk):

    #print(kk)

    mm = 0
    EA = np.zeros(M)

    sys.stdout.write(f"Processing item {kk: <1}... ")
    sys.stdout.flush()

    if((kk+1)%24==0):
        print()

    for ii in range(0, samples):
        for jj in range(ii + 1, samples):
            x1 = np.loadtxt('./conf_%.3d_%.5d.txt' % (ii, time[kk]), unpack=True)
            x2 = np.loadtxt('./conf_%.3d_%.5d.txt' % (jj, time[kk]), unpack=True)

            EA[mm] = np.dot(x1, x2) / float(len(x1))
            mm += 1

    chi_bac = np.sum(EA**2) * float(len(x1)) / float(M)

    return chi_bac

#####################################################

data = list(range(len(time)))


num_workers = 24

with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
    chi = list(executor.map(process_item, data))

print()

for ii in range(0, len(time)):
    time[ii] *= 10

#data = list(zip(time[1:], chi[1:]))
#np.savetxt('ztest.dat'.format(jj), data)

data = list(zip(time, chi))
np.savetxt('ztest.dat'.format(jj), data)
