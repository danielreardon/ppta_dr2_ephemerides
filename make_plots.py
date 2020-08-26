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
from astropy.coordinates import SkyCoord, ICRS, BarycentricTrueEcliptic, Galactic
#from parfile_optimiser import read_par

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

#datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/ecliptic/'
#outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/output/ecliptic/'
datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/final/tempo2/'
outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/final/tempo2/output/'

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
        if not 'ecliptic' in datadir:
            data_orig = datadir + 'shapiro/' +  psr + '.par.out'
            data_noshap = datadir + '/shapiro/' +  psr + '.par.no_shapiro.out'
        else:
            data_orig = datadir + 'shapiro/' +  psr + '_ecliptic.par.out'
            data_noshap = datadir + '/shapiro/' +  psr + '_ecliptic.par.no_shapiro.out'

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
output_files = sorted(glob.glob('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/final/tempo2/*.out'))


dot_size = []
dot_names = []
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

    dot_names.append(psrname)
    dot_size.append(wrms(ydata, 1/new_errs**2))

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
    plt.savefig('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/final/tempo2/output/' + psrname + '_res.pdf')
    plt.show()


sys.exit()

"""
Sky plot
"""



datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/'
outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/output/'
parfiles = sorted(glob.glob(datadir + '*.par'))

datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2e/'
outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2e/output/'
parfiles2 = np.array(sorted(glob.glob(datadir + '*.par')))

fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)
factor = 100000000000000000000
rad_to_mas = 206264806.24709636

scale=3500000000000

shift = 8
dot_names = np.array(dot_names).squeeze()
dot_size = np.array(dot_size).squeeze()

scale_dots = 50

done = []

ra_array = []
dec_array = []
ra2_array = []
dec2_array = []
plot_dots = []
for par in parfiles:
    # print(par)
    if 'J0437' in par:
        continue
    psrname = par.split('/')[-1].split('.')[0]
    if psrname in ''.join(parfiles2) or psrname in ''.join(done):
        continue
    else:
        done.append(psrname)
    params = read_par(par)


    if 'ecliptic' in datadir:
        params_position = read_par(par)
        c = SkyCoord(params['ELONG'], params['ELAT'], frame=BarycentricTrueEcliptic,
                 unit=(u.deg, u.deg))
    else:
        c = SkyCoord(params['RAJ'] + ' ' + params['DECJ'],
                 unit=(u.hourangle, u.deg))
        pm = ICRS(ra=c.ra.deg*u.degree, dec=c.dec.deg*u.deg,
                  pm_ra_cosdec = params['PMRA']*u.mas/u.yr,
                  pm_dec = params['PMDEC']*u.mas/u.yr)
    #Convert to Galactic

    ra = round(c.ra.value * np.pi/180, 5) + shift * np.pi / 12
    dec = round(c.dec.value * np.pi/180, 5)

    pmdec = factor * pm.pm_dec.value / rad_to_mas
    pmra = factor * pm.pm_ra_cosdec.value / rad_to_mas

    pmra = round(pmra* np.pi/180, 5)
    pmdec = round(pmdec * np.pi/180, 5)

    if ra > np.pi:
        ra = ra - 2*np.pi

    ra_array.append(ra)
    dec_array.append(dec)
    ra2_array.append(ra + pmra)
    dec2_array.append(dec + pmdec)

    plot_dots.append(dot_size[np.argwhere(dot_names == psrname)][0])

plot_dots = np.array(plot_dots).squeeze()
ax.scatter(-np.array(ra_array), dec_array, color='crimson', clip_on=False, s=scale_dots/plot_dots, alpha=0.8)
plt.quiver(-np.array(ra_array), dec_array, -np.array(ra2_array), dec2_array, color='crimson', alpha=0.8, headwidth=5, headlength=4, headaxislength=4, width=0.002, scale=scale, clip_on=True)

ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])

#ax.grid(True)
#plt.show()

datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2e/'
outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2e/output/'
parfiles = sorted(glob.glob(datadir + '*.par'))


#fig = plt.figure(figsize=(8,6))
#ax = fig.add_subplot(111, projection="mollweide")
done = []
ra_array = []
dec_array = []
ra2_array = []
dec2_array = []
plot_dots = []
for par in parfiles:
    # print(par)
    if 'J0437' in par:
        continue
    psrname = par.split('/')[-1].split('.')[0]
    params = read_par(par)

    if psrname in ''.join(done):
        continue
    else:
        done.append(psrname)


    if 'ecliptic' in datadir:
        params_position = read_par(par)
        c = SkyCoord(params['ELONG'], params['ELAT'], frame=BarycentricTrueEcliptic,
                 unit=(u.deg, u.deg))
    else:
        c = SkyCoord(params['RAJ'] + ' ' + params['DECJ'],
                 unit=(u.hourangle, u.deg))
        pm = ICRS(ra=c.ra.deg*u.degree, dec=c.dec.deg*u.deg,
                  pm_ra_cosdec = params['PMRA']*u.mas/u.yr,
                  pm_dec = params['PMDEC']*u.mas/u.yr)
    #Convert to Galactic

    ra = round(c.ra.value * np.pi/180, 5) + shift * np.pi / 12
    dec = round(c.dec.value * np.pi/180, 5)

    pmdec = factor * pm.pm_dec.value / rad_to_mas
    pmra = factor * pm.pm_ra_cosdec.value / rad_to_mas

    pmra = round(pmra* np.pi/180, 5)
    pmdec = round(pmdec * np.pi/180, 5)

    if ra > np.pi:
        ra = ra - 2*np.pi

    ra_array.append(ra)
    dec_array.append(dec)
    ra2_array.append(ra + pmra)
    dec2_array.append(dec + pmdec)

    plot_dots.append(dot_size[np.argwhere(dot_names == psrname)][0])


plot_dots = np.array(plot_dots).squeeze()
ax.scatter(-np.array(ra_array), dec_array, color='mediumblue', clip_on=False, s=scale_dots/plot_dots, alpha=0.8)
plt.quiver(-np.array(ra_array), dec_array, -np.array(ra2_array), dec2_array, color='mediumblue', alpha = 0.8, headwidth=5, headlength=4, headaxislength=4, width=0.002, scale=scale, clip_on=True)

#ax.set_xticklabels(['6h','8h','10h','12h','14h','16h','18h','20h','22h','0h','4h'])
ax.set_xticklabels(['2h','0h','22h','20h','18h','16h','14h','12h','10h','8h','6h'])

#plt.xlabel(r'Right Ascension, $\alpha$')
#plt.ylabel(r'Declination, $\delta$')

ra = 1.2098040598407 + shift * np.pi / 12
if ra > np.pi:
        ra = ra - 2*np.pi

ax.scatter([-ra], [-0.8158972587], color='mediumblue', marker='x', s=200)
#plt.plot([8*np.pi/12, 8*np.pi/12], [-np.pi/2, np.pi/2], 'k--')

nsamp = 1000
#glat = np.zeroes((1, nsamp))
glong = np.linspace(-180, 180, nsamp)

ra = []
dec = []
ra2 = []
dec2 = []
for i in range(0, nsamp):
    l = glong[i]
    c = SkyCoord(l, -5, frame=Galactic,
                 unit=(u.deg, u.deg))
    ra.append(c.icrs.ra.value)
    dec.append(c.icrs.dec.value)
    c2 = SkyCoord(l, 5, frame=Galactic,
                 unit=(u.deg, u.deg))
    ra2.append(c2.icrs.ra.value)
    dec2.append(c2.icrs.dec.value)


ra = np.array(ra)*np.pi/180  + shift * np.pi / 12
ra[(ra > np.pi)] = ra[(ra > np.pi)] - 2*np.pi
dec = np.array(dec)*np.pi/180

ra2 = np.array(ra2)*np.pi/180  + shift * np.pi / 12
ra2[(ra2 > np.pi)] = ra2[(ra2 > np.pi)] - 2*np.pi
dec2 = np.array(dec2)*np.pi/180


ra_resamp = np.linspace(min(ra), max(ra), nsamp)
index = np.argsort(ra)
dec_resamp = np.interp(ra_resamp, ra[index], dec[index])
index = np.argsort(ra2)
dec2_resamp = np.interp(ra_resamp, ra2[index], dec2[index])

plt.fill_between(-ra_resamp, dec_resamp, dec2_resamp, zorder=0, alpha=0.3, color='k')

ra = []
dec = []
ra2 = []
dec2 = []
for i in range(0, nsamp):
    l = glong[i]
    c = SkyCoord(l, -10, frame=Galactic,
                 unit=(u.deg, u.deg))
    ra.append(c.icrs.ra.value)
    dec.append(c.icrs.dec.value)
    c2 = SkyCoord(l, 10, frame=Galactic,
                 unit=(u.deg, u.deg))
    ra2.append(c2.icrs.ra.value)
    dec2.append(c2.icrs.dec.value)


ra = np.array(ra)*np.pi/180  + shift * np.pi / 12
ra[(ra > np.pi)] = ra[(ra > np.pi)] - 2*np.pi
dec = np.array(dec)*np.pi/180

ra2 = np.array(ra2)*np.pi/180  + shift * np.pi / 12
ra2[(ra2 > np.pi)] = ra2[(ra2 > np.pi)] - 2*np.pi
dec2 = np.array(dec2)*np.pi/180


ra_resamp = np.linspace(min(ra), max(ra), nsamp)
index = np.argsort(ra)
dec_resamp = np.interp(ra_resamp, ra[index], dec[index])
index = np.argsort(ra2)
dec2_resamp = np.interp(ra_resamp, ra2[index], dec2[index])

plt.fill_between(-ra_resamp, dec_resamp, dec2_resamp, zorder=0, alpha=0.1, color='k')



c = SkyCoord(0, 0, frame=Galactic,
             unit=(u.deg, u.deg))
ra = c.icrs.ra.value *np.pi/180   + shift * np.pi / 12
dec = c.icrs.dec.value  *np.pi/180
if ra > np.pi:
    ra = ra - 2*np.pi
plt.scatter(-ra, dec, color='yellow', s=200, marker='*')

#index = np.argsort(-ra2)
#plt.plot(-ra2[index], dec2[index], 'k')

plt.tight_layout()
plt.savefig('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/sky_map.pdf', bbox_inches = 'tight',
    pad_inches = 0)
plt.show()


