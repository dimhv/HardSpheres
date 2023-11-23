#!/usr/bin/python3

import numpy as np
import random

N = 1000
L = 10

x = np.random.uniform(0,L,N)
y = np.random.uniform(0,L,N)
z = np.random.uniform(0,L,N)


d_min = L*L

while(d_min < 0.5):

    x = np.random.uniform(0,L,N)
    y = np.random.uniform(0,L,N)
    z = np.random.uniform(0,L,N)

    d_min = L*L

    for ii in range(0,N):
        for jj in range(ii+1, N):

            dx = abs( x[ii] - x[jj] )
            dy = abs( y[ii] - y[jj] )
            dz = abs( z[ii] - z[jj] )

            if(dx > 0.5*L):
                dx = L - dx

            if(dy > 0.5*L):
                dy = L - dy

            if(dz > 0.5*L):
                dz = L - dz


            d = np.sqrt(dx*dx + dy*dy + dz*dz)

            if(d < d_min):
                d_min = d

print("d_min = %d" % d_min)

# Define the probability density function f(x)
def pdf(w):
    return w**(-3)

# Define the cumulative distribution function (CDF)
def cdf(w):
    return 1 / (2 * w**2)

# Define the range
w_min = 0.7
w_max = w_min / 0.45

# Generate a list of numbers within the specified range following the desired distribution
def generate_numbers_00(num_samples):
    samples = []
    while len(samples) < num_samples:
        # Generate a random number between 0 and 1 from a uniform distribution
        u = random.random()
        # Scale u to the desired range
        w = w_min + u * (w_max - w_min)
        # Calculate the corresponding probability density value
        probability = pdf(w)
        # Generate another random number to accept or reject the sample
        acceptance = random.random()
        if acceptance <= probability:
            samples.append(w)
    return samples

def generate_numbers_01(num_samples):
    samples = []
    while len(samples) < num_samples:

        u = random.random()
        if( u < 0.5 ):
            samples.append(1.4)

        else:
            samples.append(1.)

    return samples


#sigma  = generate_numbers_00(N)
sigma  = generate_numbers_01(N)

sigma = sigma/(np.sum(sigma)/N)

tt = 0
amp_fac = abs(d_min - 0.001)
with open(f'./c_00000_000_{tt:03d}.xyz', 'w') as fout_conf:
    fout_conf.write("ITEM: TIMESTEP\n")
    fout_conf.write("0\n")

    fout_conf.write("ITEM: NUMBER OF ATOMS\n")
    fout_conf.write(f"{N}\n")

    fout_conf.write("ITEM: BOX BOUNDS pp pp pp\n")
    fout_conf.write(f"0 {(L/amp_fac):.10f}\n")
    fout_conf.write(f"0 {(L/amp_fac):.10f}\n")
    fout_conf.write(f"0 {(L/amp_fac):.10f}\n")

    fout_conf.write("ITEM: ATOMS id type x y z radius Transparency\n")

    for ii in range(N):
        fout_conf.write(f"{ii} 1 ")
        fout_conf.write(f"{x[ii]/amp_fac:.10f} ")
        fout_conf.write(f"{y[ii]/amp_fac:.10f} ")
        fout_conf.write(f"{z[ii]/amp_fac:.10f} ")
        fout_conf.write(f"{0.5 * sigma[ii]:.10f} 0.2\n")

np.savetxt(f'./c_00000_000_{tt:03d}.dat', list(zip(x,y,z,sigma)))

check = 0
for ii in range(0,N):
    for jj in range(ii+1, N):

        dx = abs( x[ii] - x[jj] )
        dy = abs( y[ii] - y[jj] )
        dz = abs( z[ii] - z[jj] )

        if(dx > 0.5*L):
            dx = L - dx

        if(dy > 0.5*L):
            dy = L - dy

        if(dz > 0.5*L):
            dz = L - dz


        d = np.sqrt(dx*dx + dy*dy + dz*dz)

        if( d < 0.5*amp_fac*(sigma[ii]+sigma[jj]) ):
            print("OVERLAP DETECTED")
            check = 1
            break
    if(check == 1):
        break

if(check==0):
    print("no overlap detected, amp_fac = %g" % amp_fac )
