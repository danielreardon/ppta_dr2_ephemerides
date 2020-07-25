#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 20:34:01 2020

@author: dreardon


parfile_optimiser.py

Poiunt to a directory, or a specific .par file, and this tool will convert
to ELL1 if needed, and search for additional changes like Kopeikin terms and
Shapiro delays. It will then compute derived parameters.

"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from decimal import Decimal, InvalidOperation
import libstempo as T
import GalDynPsr
import scipy.constants as sc
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
import scipy.stats as stats
from astropy import units as u
from scipy.signal import savgol_filter
from astropy.coordinates import SkyCoord, ICRS, BarycentricTrueEcliptic

import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin/'
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rcParams.update({'font.size': 20})


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
            elif 'Finish' in line:
                return np.array(data), files
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

datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/'
outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/output/'
parfiles = sorted(glob.glob(datadir + '*.par'))

outfile = outdir + 'derived_params.txt'
#os.remove(outfile)

n_samples = 10000
# Define other useful constants
M_sun = 1.98847542e+30  # kg
Tsun = 4.926790353700459e-06  #s
rad_to_mas = 180*3600*1000/np.pi
parsec_to_m = 3.08567758e+16
sec_per_year = 86400*365.2425

print('Making some plots')

"""
Make a plot of J2241's noise
"""
data = np.loadtxt('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/2241/J2241-5236.tasc.txt', skiprows=1)

from astropy.time import Time
t = Time(data[:, 0], format='mjd')
yrs = t.byear  # observation year

plt.figure(figsize=(10,6))
plt.errorbar(yrs, (data[:, 1] - np.mean(data[:, 1]))*86400, yerr=data[:, 2]*86400, fmt='o', alpha=0.8, color='mediumblue')
p = np.polyfit(yrs, (data[:, 1] - np.mean(data[:, 1]))*86400, 10)
plt.xlabel(r'Year')
plt.ylabel(r'$\Delta T_{\rm asc}$ (s)')
plt.ylim([-0.25, 0.35])
xl = plt.xlim()
plt.plot([xl[0], xl[1]], [0, 0], color='k')
plt.xlim(xl)
plt.grid()
plt.savefig('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/J2241_orbit.pdf')
plt.show()

"""
Make a plot of Shapiro delays
"""
#shap_psrs = ['J0613-0200', 'J1017-7156','J1022+1001','J1125-6014', 'J1545-4550', 'J1600-3053', 'J1713+0747', 'J1732-5049', 'J1857+0943', 'J1909-3744', 'J2145-0750']
shap_psrs = ['J1125-6014']
#shap_psrs = ['J0613-0200', 'J1017-7156', 'J1022+1001', 'J1125-6014', 'J1600-3053', 'J1713+0747', 'J1857+0943', 'J1909-3744']
alpha=0.8

for psr in shap_psrs:
    print(psr)
    if psr in ['J1017-7156', 'J1713+0747', 'J1909-3744']:
        data_orig = datadir + 'shapiro/' +  psr + '.kop.par.out'
        data_noshap = datadir + '/shapiro/' +  psr + '.kop.par.no_shapiro.out'
    else:
        data_orig = datadir + 'shapiro/' +  psr + '.par.out'
        data_noshap = datadir + '/shapiro/' +  psr + '.par.no_shapiro.out'

    # No Shapiro data for plotting
    data, files = read_general2(data_noshap, header=True)
    bat = data[:, 0].squeeze()
    freqs = data[:, 1].squeeze()
    pre = data[:, 2].squeeze()*10**6
    post = data[:, 3].squeeze()*10**6
    errs = data[:, 4].squeeze()
    posttn = data[:, 5].squeeze()*10**6
    tndm = data[:, 6].squeeze()*10**6
    tnred = data[:, 7].squeeze()*10**6
    binphase = data[:, -1].squeeze()

    # Original data with Shapiro signal
    data_orig, files_orig = read_general2(data_orig, header=True)
    bat_orig = data_orig[:, 0].squeeze()
    freqs_orig = data_orig[:, 1].squeeze()
    pre_orig = data_orig[:, 2].squeeze()*10**6
    post_orig = data_orig[:, 3].squeeze()*10**6
    errs_orig = data_orig[:, 4].squeeze()
    posttn_orig = data_orig[:, 5].squeeze()*10**6
    tndm_orig = data_orig[:, 6].squeeze()*10**6
    tnred_orig = data_orig[:, 7].squeeze()*10**6
    binphase_orig = data_orig[:, -1].squeeze()

    xdata, ydata, new_errs, new_freqs = average_subbands(binphase, pre-tndm-tnred, errs, freqs, files)
    plt.subplots(3, 1, sharex=True, figsize=(10, 15))
    plt.subplot(3, 1, 1)
    indicies = np.argwhere(new_freqs > 2000)
    plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=2, color='mediumblue')
    indicies = np.argwhere(new_freqs < 1000)
    plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=0, color='crimson')
    indicies = np.argwhere((new_freqs > 1000)*(new_freqs < 2000))
    plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=1, color='darkcyan')

    plt.ylim((-20, 15))
    yl = plt.ylim()
    plt.plot([0.25, 0.25], yl, 'k--')
    plt.ylim(yl)
    plt.xlim((0, 1))
    plt.ylabel(r'Pre-fit residual  ($\mu s$), $M_c=0$')
    #plt.xlabel('Orbital phase')
    plt.tick_params(axis='x', which='both', labelbottom=False)
    plt.grid()

    xdata, ydata, new_errs, new_freqs = average_subbands(binphase, posttn, errs, freqs, files)
    plt.subplot(3, 1, 2)
    indicies = np.argwhere(new_freqs > 2000)
    plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=2, color='mediumblue')
    indicies = np.argwhere(new_freqs < 1000)
    plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=0, color='crimson')
    indicies = np.argwhere((new_freqs > 1000)*(new_freqs < 2000))
    plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=1, color='darkcyan')

    plt.ylim((-8, 8))
    yl = plt.ylim()
    plt.plot([0.25, 0.25], yl, 'k--')
    plt.ylim(yl)
    plt.xlim((0, 1))
    plt.ylabel(r'Post-fit residual  ($\mu s$), $M_c=0$')
    #plt.xlabel('Orbital phase')
    plt.tick_params(axis='x', which='both', labelbottom=False)
    plt.grid()


    xdata, ydata, new_errs, new_freqs = average_subbands(binphase_orig, posttn_orig, errs_orig, freqs_orig, files_orig)
    plt.subplot(3, 1, 3)
    indicies = np.argwhere(new_freqs > 2000)
    plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=2, color='mediumblue')
    indicies = np.argwhere(new_freqs < 1000)
    plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=0, color='crimson')
    indicies = np.argwhere((new_freqs > 1000)*(new_freqs < 2000))
    plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=1, color='darkcyan')
    plt.ylim((-8, 8))
    yl = plt.ylim()
    plt.plot([0.25, 0.25], yl, 'k--')
    #plt.plot([0, 1], [0, 0], 'k')
    plt.ylim(yl)
    plt.xlim((0, 1))
    plt.ylabel(r'Post-fit residual ($\mu s$)')
    plt.xlabel(r'Orbital phase')
    plt.grid()


    plt.tight_layout()
    if '1125' in psr:
        plt.savefig('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/J1125_Shapiro.pdf')
    plt.show()

    #plt.scatter(data_noshap[:, -1], data_noshap[:, 4]*10**6 - data[:, 4]*10**6)
    #plt.show()

"""
Make residual plots for each pulsar
"""
output_files = sorted(glob.glob('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2e/ecliptic/output/J*.out'))

for outfile in output_files:
    # No Shapiro data for plotting
    print(outfile)
    try:
        data, files = read_general2(outfile, header=True)
    except:
        continue

    plt.subplots(3, 1, sharex=True, figsize=(12, 12))

    psrname = outfile.split('/')[-1].split('.')[0].split('_')[0]

    errs = data[:, 4].squeeze()
    index = np.argwhere(errs < 100)

    bat = data[index, 0].squeeze()
    freqs = data[index, 1].squeeze()
    pre = data[index, 2].squeeze()*10**6
    post = data[index, 3].squeeze()*10**6
    errs = data[index, 4].squeeze()
    posttn = data[index, 5].squeeze()*10**6
    tndm = data[index, 6].squeeze()*10**6
    tnred = data[index, 8].squeeze()*10**6
    tnchrom = data[index, 10].squeeze()*10**6
    files = np.array(files).squeeze()[index]
    #binphase = data[:, -1].squeeze()

    from astropy.time import Time
    t = Time(bat, format='mjd')
    yrs = t.byear  # observation year

    zorder = np.array([0, 1, 2])

    xdata, ydata, new_errs, new_freqs = average_subbands(yrs, posttn, errs, freqs, files)
    plt.subplot(3, 1, 3)
    indicies_10 = np.argwhere(new_freqs > 2000)
    wrms_10 = wrms(ydata[indicies_10], 1/new_errs[indicies_10]**2)
    indicies_20 = np.argwhere((new_freqs > 1000)*(new_freqs < 2000))
    wrms_20 = wrms(ydata[indicies_20], 1/new_errs[indicies_20]**2)
    indicies_40 = np.argwhere(new_freqs < 1000)
    wrms_40 = wrms(ydata[indicies_40], 1/new_errs[indicies_40]**2)

    zsort = np.argsort([-wrms_10, -wrms_20, -wrms_40])

    plt.errorbar(xdata[indicies_10], ydata[indicies_10], yerr=new_errs[indicies_10], fmt='.', alpha=alpha, zorder=np.argwhere(zsort==0), color='mediumblue')
    plt.errorbar(xdata[indicies_20], ydata[indicies_20], yerr=new_errs[indicies_20], fmt='.', alpha=alpha, zorder=np.argwhere(zsort==1), color='darkcyan')
    plt.errorbar(xdata[indicies_40], ydata[indicies_40], yerr=new_errs[indicies_40], fmt='.', alpha=alpha, zorder=np.argwhere(zsort==2), color='crimson')
    yl = plt.ylim()
    #plt.plot([0, 1], [0, 0], 'k')
    leg_string = np.array([str(round_sig(wrms_10)) + r'$\,\mu\,$s', str(round_sig(wrms_20)) + r'$\,\mu\,$s',str(round_sig(wrms_40)) + r'$\,\mu\,$s'])
    plt.legend(leg_string, loc='upper left', framealpha=0.4, prop={'size': 17})
    plt.ylim(yl)
    plt.ylabel(r'Whitened residual ($\mu$\,s)')
    plt.xlim((1993.8, 2018.7))
    plt.xticks(ticks=[1994, 1996, 1998, 2000, 2002, 2004, 2006,2008,2010,2012,2014,2016,2018])
    plt.xlabel(r'Year')
    plt.grid()

    xdata, ydata, new_errs, new_freqs = average_subbands(yrs, post, errs, freqs, files)
    plt.subplot(3, 1, 1)
    indicies_10 = np.argwhere(new_freqs > 2000)
    wrms_10 = wrms(ydata[indicies_10], 1/new_errs[indicies_10]**2)
    indicies_20 = np.argwhere((new_freqs > 1000)*(new_freqs < 2000))
    wrms_20 = wrms(ydata[indicies_20], 1/new_errs[indicies_20]**2)
    indicies_40 = np.argwhere(new_freqs < 1000)
    wrms_40 = wrms(ydata[indicies_40], 1/new_errs[indicies_40]**2)

    zsort = np.argsort([-wrms_10, -wrms_20, -wrms_40])

    plt.errorbar(xdata[indicies_10], ydata[indicies_10], yerr=new_errs[indicies_10], fmt='.', alpha=alpha, zorder=np.argwhere(zsort==0), color='mediumblue')
    plt.errorbar(xdata[indicies_20], ydata[indicies_20], yerr=new_errs[indicies_20], fmt='.', alpha=alpha, zorder=np.argwhere(zsort==1), color='darkcyan')
    plt.errorbar(xdata[indicies_40], ydata[indicies_40], yerr=new_errs[indicies_40], fmt='.', alpha=alpha, zorder=np.argwhere(zsort==2), color='crimson')
    plt.title('PSR ' + psrname)

    yl = plt.ylim()
    plt.ylim(yl)
    plt.xlim((1993.8, 2018.7))
    plt.xticks(ticks=[1994, 1996, 1998, 2000, 2002, 2004, 2006,2008,2010,2012,2014,2016,2018])
    plt.ylabel(r'Post-fit residual ($\mu$\,s)')
    plt.tick_params(axis='x', which='both', labelbottom=False)
    leg_string = np.array([str(round_sig(wrms_10)) + r'$\,\mu\,$s', str(round_sig(wrms_20)) + r'$\,\mu\,$s',str(round_sig(wrms_40)) + r'$\,\mu\,$s'])
    plt.legend(leg_string, loc='upper left', framealpha=0.4, prop={'size': 17})
    plt.grid()

    xdata, ydata, new_errs, new_freqs = average_subbands(yrs, post-tndm-tnchrom - np.average(post-tndm-tnchrom, weights=1/errs**2), errs, freqs, files)
    plt.subplot(3, 1, 2)
    indicies_10 = np.argwhere(new_freqs > 2000)
    wrms_10 = wrms(ydata[indicies_10], 1/new_errs[indicies_10]**2)
    indicies_20 = np.argwhere((new_freqs > 1000)*(new_freqs < 2000))
    wrms_20 = wrms(ydata[indicies_20], 1/new_errs[indicies_20]**2)
    indicies_40 = np.argwhere(new_freqs < 1000)
    wrms_40 = wrms(ydata[indicies_40], 1/new_errs[indicies_40]**2)

    zsort = np.argsort([-wrms_10, -wrms_20, -wrms_40])

    plt.errorbar(xdata[indicies_10], ydata[indicies_10], yerr=new_errs[indicies_10], fmt='.', alpha=alpha, zorder=np.argwhere(zsort==0), color='mediumblue')
    plt.errorbar(xdata[indicies_20], ydata[indicies_20], yerr=new_errs[indicies_20], fmt='.', alpha=alpha, zorder=np.argwhere(zsort==1), color='darkcyan')
    plt.errorbar(xdata[indicies_40], ydata[indicies_40], yerr=new_errs[indicies_40], fmt='.', alpha=alpha, zorder=np.argwhere(zsort==2), color='crimson')

    yl = plt.ylim()
    plt.ylim(yl)
    plt.xlim((1993.8, 2018.7))
    plt.xticks(ticks=[1994, 1996, 1998, 2000, 2002, 2004, 2006,2008,2010,2012,2014,2016,2018])
    plt.ylabel(r'Post-fit residual\\ \\ $-$DM$(t)$ $-$CN$(t)$  ($\mu$\,s)')
    #plt.xlabel('Orbital phase')
    plt.tick_params(axis='x', which='both', labelbottom=False)
    leg_string = np.array([str(round_sig(wrms_10)) + r'$\,\mu\,$s', str(round_sig(wrms_20)) + r'$\,\mu\,$s',str(round_sig(wrms_40)) + r'$\,\mu\,$s'])
    plt.legend(leg_string, loc='upper left', framealpha=0.4, prop={'size': 17})
    plt.grid()



    plt.tight_layout()
    plt.savefig('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/output/' + psrname + '_res.pdf')
    plt.show()
