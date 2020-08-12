#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 22:09:46 2020

@author: dreardon
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

n_samples = 1000000

eps1 = 4.4289220051202546193e-08
eps1_err = 0.00000000595139388989

eps2 = -9.776878955568826327e-08
eps2_err = 0.00000000339510407484

i = 86.465574649843698139 * np.pi/180
i_err = 0.00183480226343178329 * np.pi/180

asini = 1.8979911727791108372*sc.c
asini_err = 0.00000001274621563121*sc.c

eps1 = np.random.normal(loc=eps1, scale=eps1_err, size=n_samples)
eps2 = np.random.normal(loc=eps2, scale=eps2_err, size=n_samples)
asini = np.random.normal(loc=asini, scale=asini_err, size=n_samples)
i = np.random.normal(loc=i, scale=i_err, size=n_samples)


ratio = np.sqrt(1 - (eps1**2 + eps2**2))

difference = asini/(np.sin(i)) - asini/(np.sin(i)) * ratio

plt.hist(difference, bins=15)

print(np.median(difference), np.std(difference))