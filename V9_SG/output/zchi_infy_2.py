#!/usr/bin/python3

import numpy as np
import sys

#chi_infy = 291.0431834291715
#chi_infy = 156.11911854087066
#chi_infy = 97.6538216260655
#chi_infy = 51.495947198645425
chi_infy = 21.407177227243782
time, chi = np.loadtxt('./ztest.dat', unpack=True)

chi = np.ones(len(chi)) - chi/chi_infy

data = list(zip(time[1:], chi[1:]))
np.savetxt('ztest_rlx.dat', data)
