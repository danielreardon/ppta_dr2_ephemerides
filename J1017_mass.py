#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 01:02:30 2020

@author: dreardon


J1017's mass limit
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from os import path
from decimal import Decimal, InvalidOperation
import libstempo as T
import GalDynPsr
import scipy.constants as sc
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
from astropy import units as u
from scipy.signal import savgol_filter
from astropy.coordinates import SkyCoord, ICRS, BarycentricTrueEcliptic
import sys
sys.path.insert(0, '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/parameterComparisonScripts/')


"""
Definitions
"""

from math import log10, floor
def round_sig(x, sig=2, small_value=1.0e-9):
    return round(x, sig - int(floor(log10(max(abs(x), abs(small_value))))) - 1)

def wrms(data, weights):
    """
    Given some data and weights, return the weighted rms
    """
    return np.sqrt(np.cov(np.array(data).squeeze(),
                          aweights=np.array(weights).squeeze()))

def read_general2(filename, header=False):
    """
    Reads a general2 output into a numpy array
    """
    with open(filename, "r") as file:
        data = []
        files = []
        for line in file:
            if not header and not ('Finish' in line):
                files.append([(line.split('\t')[0])])
                data.append([float(x) for x in
                             (line.replace('\n', '').split('\t')[1:])])
            elif 'Starting general2 plugin' in line:
                header = False
    return np.array(data), files

def average_subbands(times, residuals, errs, freqs, files):
    """
    Given an array of residuals and frequencies, average the subbands together
    """

    times_avg = []
    residuals_avg = []
    errs_avg = []
    freqs_avg = []
    files = np.array(files).squeeze()

    ref_freq = -np.inf
    for file in np.unique(files):
        indicies = np.argwhere(files == file)

        times_avg.append(np.average(times[indicies],
                                    weights=1/(np.array(errs[indicies])**2)))
        residuals_avg.append(np.average(residuals[indicies],
                                        weights=1/(np.array(errs[indicies])**2)))
        errs_avg.append(np.average(errs[indicies],
                                   weights=1/(np.array(errs[indicies])**2))/np.sqrt(len(errs[indicies])))
        freqs_avg.append(np.average(freqs[indicies],
                                    weights=1/(np.array(errs[indicies])**2)))

    return np.array(times_avg).squeeze(), np.array(residuals_avg).squeeze(), \
        np.array(errs_avg).squeeze(), np.array(freqs_avg).squeeze()


def read_par(parfile):
    """
    Reads a par file and return a dictionary of parameter names and values
    """

    par = {}
    ignore = ['DMMODEL', 'DMOFF', "DM_", "CM_", 'CONSTRAIN',
              'CORRECT_TROPOSPHERE', 'PLANET_SHAPIRO', 'DILATEFREQ',
              'TIMEEPH', 'MODE', 'TZRMJD', 'TZRSITE', 'TZRFRQ', 'EPHVER',
              'T2CMETHOD']

    file = open(parfile, 'r')
    for line in file.readlines():
        err = None
        p_type = None
        sline = line.split()
        if len(sline) == 0 or line[0] == "#" or line[0:1] == "C " \
           or sline[0] in ignore:
            continue

        param = sline[0]
        if param == "E":
            param = "ECC"

        if param[0:3] == 'TNE' or param[0:4] == 'JUMP':
            param = sline[0] + '_' + sline[1] + '_' + sline[2]
            val = sline[3]
        else:
            val = sline[1]

        if len(sline) == 3 and sline[2] not in ['0', '1']:
            err = sline[2].replace('D', 'E')
        elif len(sline) == 4:
            err = sline[3].replace('D', 'E')

        try:
            val = int(val)
            p_type = 'd'
        except ValueError:
            try:
                val = float(Decimal(val.replace('D', 'E')))
                if 'e' in sline[1] or 'E' in sline[1].replace('D', 'E'):
                    p_type = 'e'
                else:
                    p_type = 'f'
            except InvalidOperation:
                p_type = 's'

        par[param] = val
        if err:
            par[param+"_ERR"] = float(err)

        if p_type:
            par[param+"_TYPE"] = p_type

    file.close()

    return par


def make_ell1(par, output):
    params = read_par(par)
    om = params['OM']*np.pi/180
    ecc = params['ECC']
    t0 = Decimal(params['T0'])
    nb = 2*np.pi/params['PB']

    eps1 = ecc*np.sin(om)
    eps2 = ecc*np.cos(om)
    tasc = t0 - Decimal(om/nb)

    file1 = open(par, "r")
    file2 = open(output, "w")

    for line in file1.readlines():

        # Separate each item in the line
        items = line.strip().split()

        if items[0] == 'ECC':
            items[0] = 'EPS1'
            items[1] = str(eps1)
        elif items[0] == 'OM':
            items[0] = 'EPS2'
            items[1] = str(eps2)
        elif items[0] == 'T0':
            items[0] = 'TASC'
            items[1] = str(tasc)

        string = '\t'.join(items)

        # Write to the file if conditions are met
        file2.write(string + "\n")

    print("wrote: ", output)
    return


def is_valid(array):
    """
    Returns boolean array of values that are finite an not nan
    """
    return np.isfinite(array)*(~np.isnan(array))


"""
Start of code
"""

datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/final/tempo2/'
outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/final/tempo2/output/'

parfiles = sorted(glob.glob(datadir + 'J1017*.par'))

psrnames = np.array(['J1017-7156'])
params = read_par(parfiles[0])


outfile = outdir + 'derived_params.txt'
if path.exists(outfile):
    os.remove(outfile)

n_samples = 10000000
# Define other useful constants
M_sun = 1.98847542e+30  # kg
Tsun = 4.926790353700459e-06  #s
rad_to_mas = 180*3600*1000/np.pi
parsec_to_m = 3.08567758e+16
sec_per_year = 86400*365.2425


Omega = np.random.rand(n_samples)*2*np.pi
cosi = np.random.rand(n_samples)
i_prior = np.arccos(cosi)*180/np.pi
i[i<0] = -i[i<0]

# XDOT
xdot = np.random.normal(loc=params["XDOT"], scale=params["XDOT_ERR"], size=n_samples)
a1 = np.random.normal(loc=params["A1"], scale=params["A1_ERR"], size=n_samples)
pmelat = np.random.normal(loc=params["PMELAT"], scale=params["PMELAT_ERR"], size=n_samples)
pmelong = np.random.normal(loc=params["PMELONG"], scale=params["PMELONG_ERR"], size=n_samples)
pm_tot = np.sqrt(pmelat**2 + pmelong**2)
pm = pm_tot/(sec_per_year*rad_to_mas)
i_limit =  np.abs(180/np.pi * np.arctan(a1 * pm / xdot))
i_lim_95 = np.percentile(i_limit, q=95)
sini_lim_95 = np.sin(np.percentile(i_limit, q=84)*np.pi/180)

A = (-pmelong/(sec_per_year*rad_to_mas))*np.sin(Omega) + (-pmelat/(sec_per_year*rad_to_mas))*np.cos(Omega)
i = np.arctan2(a1 * A, xdot)*180/np.pi
i[i<0] = -i[i<0]
#plt.figure()
plt.hist(i_prior, bins=1000, alpha=0.5)
#plt.show()

#plt.figure()
plt.hist(i, bins=1000, alpha=0.5)
#plt.show()

# H3 STIG
mtot2 = 0
h3 = params["H3"]*10**6
h3_err = params["H3_ERR"]*10**6
stig = params["STIG"]
stig_err = params["STIG_ERR"]
h3 = np.random.normal(loc=h3, scale=h3_err, size=n_samples)
stig = np.random.normal(loc=stig, scale=stig_err, size=n_samples)
h4 = stig * h3
sini = 2 * h3 * h4 / ( h3**2 + h4**2 )
m2 = h3**4 / h4**3
m2 = m2/(Tsun*10**6)
cut = np.argwhere((m2 > mass_func) * (m2 < 1.4))
m2 = m2[cut]
sini = sini[cut]
i = i[cut]
i_prior = i_prior[cut]
inc = np.arcsin(sini)*180/np.pi
mtot2 = (m2 * sini)**3 / mass_func
mp = np.sqrt(mtot2[is_valid(mtot2)*is_valid(m2)]) - m2[is_valid(mtot2)*is_valid(m2)]
mp = mp[is_valid(mp)]
mtot2 = mtot2[is_valid(mtot2)]
mtot = np.sqrt(mtot2[(mtot2 > 0)])
cut = is_valid(inc).squeeze()
inc = inc[cut]
i = i[cut]
i_prior = i_prior[cut]

#plt.figure()
plt.hist(inc, bins=1000, alpha=0.5)
plt.xlim(0, 90)
plt.show()


# Derive Mtot from OMDOT
ecc = params['ECC']
ecc_err = params['ECC_ERR']
ecc_posterior = np.random.normal(loc=params["ECC"], scale=params["ECC_ERR"], size=n_samples)
omdot_posterior = np.random.normal(loc=params["OMDOT"], scale=params["OMDOT_ERR"], size=n_samples)
cut = np.argwhere(omdot_posterior > 0)
omdot_posterior = omdot_posterior[cut]
ecc_posterior = ecc_posterior[cut]
Msun = 1.989*10**30  # kg
Tsun = sc.G * Msun / (sc.c**3)
Mtot = [1.6 if (np.std(mtot2)>0.2 or mtot2==0) else np.median(mtot2)]
Mtot = Mtot[0]
n = 2*np.pi/(params['PB']*86400)  # s
omdot_posterior = omdot_posterior / 86400 / 365.2425 * np.pi/180
Mtot = (omdot_posterior*(1 - ecc_posterior**2) / (3 * n**(5/3)))**(3/2)/Tsun

i = np.sort(i)
i_prior = np.sort(i_prior)
inc = np.sort(inc)

bins = 1000

plt.figure()
n, b, p = plt.hist(i_prior, bins=bins, alpha=0.5, range=(0, 90))
n2, b2, p = plt.hist(i, bins=bins, alpha=0.5, range=(0, 90))
n3, b3, p = plt.hist(inc, bins=bins, alpha=0.5, range=(0, 90))
plt.show()

plt.figure()
plt.plot(np.linspace(0,90,bins), n/n_samples)
plt.plot(np.linspace(0,90,bins), n2/n_samples)
plt.plot(np.linspace(0,90,bins), n3/n_samples)
combined = n*n2*n3
plt.plot(np.linspace(0,90,bins), combined/n_samples)
plt.xlim((0, 90))
plt.show()


