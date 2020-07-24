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
import scipy.stats as stats
from astropy import units as u
from scipy.signal import savgol_filter
from astropy.coordinates import SkyCoord, ICRS, BarycentricTrueEcliptic



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

datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/'
outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/output/'
parfiles = sorted(glob.glob(datadir + '*.par'))

outfile = outdir + 'derived_params.txt'
plots = True
#os.remove(outfile)

n_samples = 10000
# Define other useful constants
M_sun = 1.98847542e+30  # kg
Tsun = 4.926790353700459e-06  #s
rad_to_mas = 180*3600*1000/np.pi
parsec_to_m = 3.08567758e+16
sec_per_year = 86400*365.2425

if plots:
    print('Making some plots')

    """
    Make a plot of J2241's noise
    """
    data = np.loadtxt('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2/2241/J2241-5236.tasc.txt', skiprows=1)

    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 18}

    matplotlib.rc('font', **font)

    plt.figure(figsize=(10,6))
    plt.errorbar(data[:, 0], (data[:, 1] - np.mean(data[:, 1]))*86400, yerr=data[:, 2]*86400, fmt='o', alpha=0.8)
    p = np.polyfit(data[:, 0], (data[:, 1] - np.mean(data[:, 1]))*86400, 10)
    plt.xlabel('MJD')
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
        plt.ylabel(r'Residual ($\mu s$)')
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
        plt.ylabel(r'Residual ($\mu s$)')
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
        plt.ylabel(r'Residual ($\mu s$)')
        plt.xlabel('Orbital phase')
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

    output_files = sorted(glob.glob('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2e/output/*.out'))

    for outfile in output_files:
        # No Shapiro data for plotting
        print(outfile)
        try:
            data, files = read_general2(outfile, header=True)
        except:
            continue

        plt.subplots(3, 1, sharex=True, figsize=(12, 12))

        psrname = outfile.split('/')[-1].split('.')[0].split('_')[0]
        bat = data[:, 0].squeeze()
        freqs = data[:, 1].squeeze()
        pre = data[:, 2].squeeze()*10**6
        post = data[:, 3].squeeze()*10**6
        errs = data[:, 4].squeeze()
        posttn = data[:, 5].squeeze()*10**6
        tndm = data[:, 6].squeeze()*10**6
        tnred = data[:, 8].squeeze()*10**6
        tnchrom = data[:, 10].squeeze()*10**6
        #binphase = data[:, -1].squeeze()

        from astropy.time import Time
        t = Time(bat, format='mjd')
        yrs = t.byear  # observation year

        xdata, ydata, new_errs, new_freqs = average_subbands(yrs, posttn, errs, freqs, files)
        plt.subplot(3, 1, 3)
        indicies = np.argwhere(new_freqs > 2000)
        wrms_10 = wrms(ydata[indicies], 1/new_errs[indicies]**2)
        plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=2, color='mediumblue')
        indicies = np.argwhere((new_freqs > 1000)*(new_freqs < 2000))
        wrms_20 = wrms(ydata[indicies], 1/new_errs[indicies]**2)
        plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=1, color='darkcyan')
        indicies = np.argwhere(new_freqs < 1000)
        wrms_40 = wrms(ydata[indicies], 1/new_errs[indicies]**2)
        plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=0, color='crimson')
        yl = plt.ylim()
        #plt.plot([0, 1], [0, 0], 'k')
        plt.legend([str(round_sig(wrms_10)) + r'$\,\mu\,$s', str(round_sig(wrms_40)) + r'$\,\mu\,$s', str(round_sig(wrms_20)) + r'$\,\mu\,$s'], loc='upper left', framealpha=0.4)
        plt.ylim(yl)
        plt.ylabel(r'Residual ($\mu s$)')
        plt.xlim((1993.8, 2018.7))
        plt.xticks(ticks=[1994, 1996, 1998, 2000, 2002, 2004, 2006,2008,2010,2012,2014,2016,2018])
        plt.xlabel('MJD')
        plt.grid()

        xdata, ydata, new_errs, new_freqs = average_subbands(yrs, post, errs, freqs, files)
        plt.subplot(3, 1, 1)
        indicies = np.argwhere(new_freqs > 2000)
        wrms_10 = wrms(ydata[indicies], 1/new_errs[indicies]**2)
        plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=2, color='mediumblue')
        indicies = np.argwhere((new_freqs > 1000)*(new_freqs < 2000))
        wrms_20 = wrms(ydata[indicies], 1/new_errs[indicies]**2)
        plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=1, color='darkcyan')
        indicies = np.argwhere(new_freqs < 1000)
        wrms_40 = wrms(ydata[indicies], 1/new_errs[indicies]**2)
        plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=0, color='crimson')
        plt.title('PSR ' + psrname)

        yl = plt.ylim()
        plt.ylim(yl)
        plt.xlim((1993.8, 2018.7))
        plt.xticks(ticks=[1994, 1996, 1998, 2000, 2002, 2004, 2006,2008,2010,2012,2014,2016,2018])
        plt.ylabel(r'Residual ($\mu s$)')
        plt.tick_params(axis='x', which='both', labelbottom=False)
        plt.legend([str(round_sig(wrms_10)) + r'$\,\mu\,$s', str(round_sig(wrms_40)) + r'$\,\mu\,$s', str(round_sig(wrms_20)) + r'$\,\mu\,$s'], loc='upper left', framealpha=0.4)
        plt.grid()

        xdata, ydata, new_errs, new_freqs = average_subbands(yrs, post-tndm-tnchrom, errs, freqs, files)
        plt.subplot(3, 1, 2)
        indicies = np.argwhere(new_freqs > 2000)
        wrms_10 = wrms(ydata[indicies], 1/new_errs[indicies]**2)
        plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=2, color='mediumblue')
        indicies = np.argwhere((new_freqs > 1000)*(new_freqs < 2000))
        wrms_20 = wrms(ydata[indicies], 1/new_errs[indicies]**2)
        plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=1, color='darkcyan')
        indicies = np.argwhere(new_freqs < 1000)
        wrms_40 = wrms(ydata[indicies], 1/new_errs[indicies]**2)
        plt.errorbar(xdata[indicies], ydata[indicies], yerr=new_errs[indicies], fmt='.', alpha=alpha, zorder=0, color='crimson')

        yl = plt.ylim()
        plt.ylim(yl)
        plt.xlim((1993.8, 2018.7))
        plt.xticks(ticks=[1994, 1996, 1998, 2000, 2002, 2004, 2006,2008,2010,2012,2014,2016,2018])
        plt.ylabel(r'Residual ($\mu s$)')
        #plt.xlabel('Orbital phase')
        plt.tick_params(axis='x', which='both', labelbottom=False)
        plt.legend([str(round_sig(wrms_10)) + r'$\,\mu\,$s', str(round_sig(wrms_40)) + r'$\,\mu\,$s', str(round_sig(wrms_20)) + r'$\,\mu\,$s'], loc='upper left', framealpha=0.4)
        plt.grid()



        plt.tight_layout()
        plt.savefig('/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/publish_collection/dr2e/output/' + psrname + '_res.pdf')
        plt.show()


for par in parfiles:
    #print(par)
    if 'J0437' in par:
        continue
    psrname = par.split('/')[-1].split('.')[0]

    #print('=========================')
    print(par.split('/')[-1])
    #print('=========================')
    #print(par.split('/')[-1].split('.')[0])
    params = read_par(par)

    with open(outfile, 'a+') as f:
        f.write(par.split('/')[-1] + '\n')


    if 'ecliptic' in datadir:
        params_position = read_par(par)
        c = SkyCoord(params['ELONG'], params['ELAT'], frame=BarycentricTrueEcliptic,
                 unit=(u.deg, u.deg))
    else:
        c = SkyCoord(params['RAJ'] + ' ' + params['DECJ'],
                 unit=(u.hourangle, u.deg))
    #Convert to Galactic

    ldeg = c.galactic.l.value
    sigl = 0
    bdeg = c.galactic.b.value
    sigb = 0

    if not 'ecliptic' in datadir:
        #Convert to Ecliptic
        #ecliptic transformation:
        elat = c.barycentrictrueecliptic.lat.value
        elon = c.barycentrictrueecliptic.lon.value
        #proper motion transformation:
        pm = ICRS(ra=c.ra.deg*u.degree, dec=c.dec.deg*u.deg,
                  pm_ra_cosdec = params['PMRA']*u.mas/u.yr,
                  pm_dec = params['PMDEC']*u.mas/u.yr)
        pm_ecliptic = pm.transform_to(BarycentricTrueEcliptic)
        pmlat = pm_ecliptic.pm_lat.value
        pmlon = pm_ecliptic.pm_lon_coslat.value
        with open(outfile, 'a+') as f:
            f.write("ELAT" + '\t' + str(elat) + '\n')
            f.write("ELONG" + '\t' + str(elon) + '\n')
            f.write("PMELAT" + '\t' + str(pmlat) + '\n')
            f.write("PMELONG" + '\t' + str(pmlon) + '\n')


    if 'BINARY' in params.keys():
        mass_func = (4*np.pi**2/sc.G) * (params['A1']*sc.c)**3/(params['PB']*86400)**2
        mass_func = mass_func/M_sun
        #print('Mass function = ', mass_func)
        with open(outfile, 'a+') as f:
            f.write("MASS_FUNC" + '\t' + str(mass_func) + '\n')

    # Check if pulsar needs ELL1:
    if 'ECC' in params.keys():
        if params['A1']*params['ECC']**2 < \
           0.01*params['TRES']*10**-6/np.sqrt(params['NTOA']):
            print(params['PSRJ'], "WARNING! Needs ELL1 model")
            # make ELL1
            make_ell1(par, outdir + params['PSRJ'] + '_ell1.par')
    elif 'EPS1' in params.keys():
        ecc = np.sqrt(params['EPS1']**2 + params['EPS2']**2)
        if params['A1']*ecc**2 > \
           0.01*params['TRES']*10**-6/np.sqrt(params['NTOA']):
            print('WARNING! Pulsar should NOT have ELL1 model')

    # Check if pulsar requires Kopeikin terms
#    if 'XDOT' in params.keys() and 'OMDOT' in params.keys():
        #print('Check for Kopeikin terms')


    if 'PX' in params.keys():
        if 'PX_ERR' in params.keys():
            D_prior = 1/np.random.normal(loc=params["PX"],
                                           scale=params["PX_ERR"], size=n_samples) # observed
            #plt.hist(D_prior, bins=100)
            #plt.show()
            D_array = D_prior[(D_prior > 0)*(D_prior < 100)]
            dkpc = np.median(D_prior[(D_prior > 0)*(D_prior < 100)])
            sigd = np.std(D_prior[(D_prior > 0)*(D_prior < 100)])
            #sigd = dkpc*params["PX_ERR"]/params["PX"]
            #print("Parallax distance (kpc) = ", round(dkpc, 3), " +/- ", round(sigd, 3))
            with open(outfile, 'a+') as f:
                f.write("D_PX(med/16th/84th)" + '\t' + str(np.median(D_array)) + '\t' + str(np.percentile(D_array, q=16)) + '\t' + str(np.percentile(D_array, q=84)) + '\n')


    if 'PBDOT' in params.keys():
        #print(' ')
        #print('=== Computing Shklovskii distance ===')
        # Pulsar distance from Shklovskii effect
        pbdot_posterior = np.random.normal(loc=params["PBDOT"],
                                           scale=params["PBDOT_ERR"], size=n_samples)  # observed
        pbdot_grav = 0
        if 'PX' in params.keys():
            D_prior = 1/np.random.normal(loc=params["PX"],
                                           scale=params["PX_ERR"], size=n_samples) # observed
            dkpc = np.median(D_prior[(D_prior > 0)*(D_prior < 100)])
            #sigd = np.std(D_prior[(D_prior > 0)*(D_prior < 100)])
            sigd = dkpc*params["PX_ERR"]/params["PX"]
        else:
            dkpc = 1
            sigd = 0.2*dkpc
        pb = params['PB']
        pb_err = params['PB_ERR']
        if 'ecliptic' in datadir:
            pmelat = params['PMELAT']
            pmelong = params['PMELONG']
            pm_tot = np.sqrt(pmelat**2 + pmelong**2)
            pm = pm_tot/(sec_per_year*rad_to_mas)
        else:
            pmra = params['PMRA']
            pmdec = params['PMDEC']
            pm_tot = np.sqrt(pmelat**2 + pmelong**2)
            pm = pm_tot/(sec_per_year*rad_to_mas)


        # Expected Shklovskii and Galactic terms
        Ex_pl =  GalDynPsr.modelLb.Expl(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term parallel to the Galactic plane
        Ex_z =  GalDynPsr.modelLb.Exz(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term perpendicular to the Galactic plane
        errpl = np.abs(0.03*Ex_pl) # 3% for circular Bovy et al. (2012)
        errz = np.abs(0.1*Ex_z) # 10% for vertical Holmberg and Flynn (2004)
        pbdot_gal = GalDynPsr.pdotint.PdotGal(Ex_pl, Ex_z, pb*86400)
        pbdot_gal_err =GalDynPsr.pdotint.ErrPdotGal(Ex_pl, errpl, Ex_z, errz, pb*86400, pb_err*86400)

        #print("Observed Pbdot = ", np.mean(pbdot_posterior), " +/- ", np.std(pbdot_posterior))
        #print("Galactic Pbdot contribution = ", pbdot_gal, " +/- ", pbdot_gal_err)
        pbdot_gal_posterior = np.random.normal(loc=pbdot_gal, scale=pbdot_gal_err, size=n_samples)  # sample randomly from pbdot_gal distribution
        pbdot_shklovskii = np.subtract(pbdot_posterior, pbdot_gal_posterior) - pbdot_grav
        #print("Observed Shklovskii contribution = ", np.mean(pbdot_shklovskii), " +/- ", np.std(pbdot_shklovskii))

        D_posterior = sc.c*pbdot_shklovskii/(pm**2*pb*86400)/parsec_to_m/1000
        D = np.mean(D_posterior)
        D_err = np.std(D_posterior)
        #print("Shklovskii distance (kpc) = ", round(D, 3), " +/- ", round(D_err,3))
        with open(outfile, 'a+') as f:
            #f.write("D_SHK" + '\t' + str(D) + '\t' + str(D_err) + '\n')
            f.write("D_SHK(med/16th/84th)" + '\t' + str(np.median(D_posterior)) + '\t' + str(np.percentile(D_posterior, q=16)) + '\t' + str(np.percentile(D_posterior, q=84)) + '\n')
        if 'PX' in params.keys():
            Davg = np.average([dkpc, D], weights=[1/sigd, 1/D_err])
            Davg_err = 1
            #print("Combined distance (kpc) = ", round(Davg, 3), " +/- ", round(Davg_err,3))

    D = 1
    if 'PBDOT' in params.keys():
        if 'PX' in params.keys():
            D = Davg
        else:
            D = D
    elif 'PX' in params.keys():
        if 'PX_ERR' in params.keys():
            D_prior = 1/np.random.normal(loc=params["PX"],
                                               scale=params["PX_ERR"], size=n_samples) # observed
            D = np.median(D_prior[(D_prior > 0)*(D_prior < 100)])

#    if 'F2' in params.keys():
#        #print(' ')
#        #rint('=== Computing radial velocity from F2 ===')
#        f0 = params['F0']
#        try:
#            f0err = params['F0_ERR']
#        except KeyError:
#            continue
#        f1 = params['F1']
#        f1_err = params['F1_ERR']
#        f2 = params['F2']
#        f2_err = params['F2_ERR']
#        pmra = params['PMRA']
#        pmdec = params['PMDEC']
#        pm = np.sqrt(pmra**2 + pmdec**2)/(sec_per_year*rad_to_mas)
#        p0_obs = 1/f0
#        p1_obs = -f1/f0**2
#        p2_obs = f1*(2*f1/f0**3) - (1/f0**2)*(f2)
#        #print(f2)
#        Ex_pl =  GalDynPsr.modelLb.Expl(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term parallel to the Galactic plane
#        Ex_z =  GalDynPsr.modelLb.Exz(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term perpendicular to the Galactic plane
#
#
#        dv_r = np.Inf
#        v_r = 0
#        # iterate Doppler correction
#        while dv_r > np.abs(0.01*v_r):
#            p0_int = p0_obs/(1 - v_r/sc.c) # approximate Doppler correction
#            p1_int = p1_obs - GalDynPsr.pdotint.PdotGal(Ex_pl, Ex_z, p0_int) # Shklovskii correction
#            v_r_new =  (2*p1_int*D*parsec_to_m - p2_obs*(sc.c/pm**2))/(3*p0_int) # km/s ... assuming p2_int=0
#            dv_r = np.abs(v_r - v_r_new)
#            v_r = v_r_new
#        #print("Observed P2 = ", p2_obs)
#        #print("Radial velocity = ", round(v_r/1000,1), " +/- ", round(f2_err/f2 * v_r/1000, 1), " km/s")
#        p2 =  (2*p1_int*D*parsec_to_m - v_r*(3*p0_int))/(sc.c/pm**2)
#        with open(outfile, 'a+') as f:
#            f.write("V_R" + '\t' + str(round(v_r,1)) + '\t' + str(round(f2_err/f2 * v_r, 1)) + '\n')

    if 'XDOT' in params.keys():
        #print(' ')
        #print('=== Computing limit on i from XDOT ===')
        xdot = np.random.normal(loc=params["XDOT"], scale=params["XDOT_ERR"], size=n_samples)
        a1 = np.random.normal(loc=params["A1"], scale=params["A1_ERR"], size=n_samples)

        if 'ecliptic' in datadir:
            pmelat = np.random.normal(loc=params["PMELAT"], scale=params["PMELAT_ERR"], size=n_samples)
            pmelong = np.random.normal(loc=params["PMELONG"], scale=params["PMELONG_ERR"], size=n_samples)
            pm_tot = np.sqrt(pmelat**2 + pmelong**2)
            pm = pm_tot/(sec_per_year*rad_to_mas)
        else:
            pmra = np.random.normal(loc=params["PMRA"], scale=params["PMRA_ERR"], size=n_samples)
            pmdec = np.random.normal(loc=params["PMDEC"], scale=params["PMDEC_ERR"], size=n_samples)
            pm_tot = np.sqrt(pmelat**2 + pmelong**2)
            pm = pm_tot/(sec_per_year*rad_to_mas)


        i_limit =  np.abs(180/np.pi * np.arctan(a1 * pm / xdot))
        #print("i <= ", int(np.ceil(i_limit)))
        with open(outfile, 'a+') as f:
            f.write("INC_LIM(med/std)" + '\t' + "<" + str(np.median(i_limit))+ '\t'+ str(np.std(i_limit)) +'\n')

    mtot2 = 0
    if 'H3' in params.keys():
        if 'H4' in params.keys():
            h3 = params["H3"]*10**6
            h3_err = params["H3_ERR"]*10**6
            h4 = params["H4"]*10**6
            h4_err = params["H4_ERR"]*10**6
            h3 = np.random.normal(loc=h3, scale=h3_err, size=n_samples)
            h4 = np.random.normal(loc=h4, scale=h4_err, size=n_samples)
            sini = 2 * h3 * h4 / ( h3**2 + h4**2 )
            m2 = h3**4 / h4**3  # in us
            m2 = m2/(Tsun*10**6)
            cut = np.argwhere((m2 > mass_func) * (m2 < 10))
            m2 = m2[cut]
            sini = sini[cut]
            inc = np.arcsin(sini)*180/np.pi
            mtot2 = (m2 * sini)**3 / mass_func
            mp = np.sqrt(mtot2[is_valid(mtot2)*is_valid(m2)]) - m2[is_valid(mtot2)*is_valid(m2)]
            mp = mp[is_valid(mp)]
            mtot2 = mtot2[is_valid(mtot2)]
            inc = inc[is_valid(inc)]

            with open(outfile, 'a+') as f:
                f.write("INC(med/16th/84th)" + '\t' + str(np.median(inc)) + '\t' + str(np.percentile(inc, q=16)) + '\t' + str(np.percentile(inc, q=84)) + '\n')
            with open(outfile, 'a+') as f:
                f.write("M2(med/16th/84th)" + '\t' + str(np.median(m2)) + '\t' + str(np.percentile(m2, q=16)) + '\t' + str(np.percentile(m2, q=84)) + '\n')
            with open(outfile, 'a+') as f:
                f.write("MP(med/16th/84th)" + '\t' + str(np.median(mp)) + '\t' + str(np.percentile(mp, q=16)) + '\t' + str(np.percentile(mp, q=84)) + '\n')
            with open(outfile, 'a+') as f:
                f.write("MTOT(med/16th/84th)" + '\t' + str(np.median(mtot2)) + '\t' + str(np.percentile(mtot2, q=16)) + '\t' + str(np.percentile(mtot2, q=84)) + '\n')

            #plt.hist(mp, bins=100)
            #plt.xlabel('Pulsar mass')
            #plt.show()
            #print('Median m2: ', np.median(m2), ", inclination: ", np.median(inc), ", Pulsar mass: ", np.median(mp))
            #print('1-sigma range:', np.percentile(mp, q=16), np.percentile(mp, q=84))
            #plt.scatter(m2, inc, alpha=0.5)
            #plt.xlabel('companion mass')
            #plt.ylabel('inclination (deg)')
            #plt.show()
        elif 'STIG' in params.keys():
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
            cut = np.argwhere((m2 > mass_func) * (m2 < 10))
            m2 = m2[cut]
            sini = sini[cut]
            inc = np.arcsin(sini)*180/np.pi
            mtot2 = (m2 * sini)**3 / mass_func
            mp = np.sqrt(mtot2[is_valid(mtot2)*is_valid(m2)]) - m2[is_valid(mtot2)*is_valid(m2)]
            mp = mp[is_valid(mp)]
            mtot2 = mtot2[is_valid(mtot2)]
            inc = inc[is_valid(inc)]

            try:
                with open(outfile, 'a+') as f:
                    f.write("INC(med/16th/84th)" + '\t' + str(np.median(inc)) + '\t' + str(np.percentile(inc, q=16)) + '\t' + str(np.percentile(inc, q=84)) + '\n')
                with open(outfile, 'a+') as f:
                    f.write("M2(med/16th/84th)" + '\t' + str(np.median(m2)) + '\t' + str(np.percentile(m2, q=16)) + '\t' + str(np.percentile(m2, q=84)) + '\n')
                with open(outfile, 'a+') as f:
                    f.write("MP(med/16th/84th)" + '\t' + str(np.median(mp)) + '\t' + str(np.percentile(mp, q=16)) + '\t' + str(np.percentile(mp, q=84)) + '\n')
                with open(outfile, 'a+') as f:
                    f.write("MTOT(med/16th/84th)" + '\t' + str(np.median(mtot2)) + '\t' + str(np.percentile(mtot2, q=16)) + '\t' + str(np.percentile(mtot2, q=84)) + '\n')
            except Exception as e:
                print(e)

            #plt.hist(mp, bins=100)
            #plt.xlabel('Pulsar mass')
            #plt.show()
            #print('Median m2: ', np.median(m2), ", inclination: ", np.median(inc), ", Pulsar mass: ", np.median(mp))
            #print('1-sigma range:', np.percentile(mp, q=16), np.percentile(mp, q=84))
            #plt.scatter(m2, inc, alpha=0.5)
            #plt.xlabel('companion mass')
            #plt.ylabel('inclination (deg)')
            #plt.show()

    if 'M2' in params.keys():
        m2 = np.random.normal(loc=params["M2"], scale=params["M2_ERR"], size=n_samples)
        if 'KIN' in params.keys():
            kin = np.random.normal(loc=params["KIN"], scale=params["KIN_ERR"], size=n_samples)
            sini = np.sin(kin*np.pi/180)
            inc=kin
        else:
            sini = np.random.normal(loc=params["SINI"], scale=params["SINI_ERR"], size=n_samples)
            inc = np.arcsin(sini)*180/np.pi
        cut = np.argwhere((m2 > mass_func) * (m2 < 10))
        m2 = m2[cut]
        sini = sini[cut]
        inc = np.arcsin(sini)*180/np.pi
        mtot2 = (m2 * sini)**3 / mass_func
        mp = np.sqrt(mtot2[is_valid(mtot2)*is_valid(m2)]) - m2[is_valid(mtot2)*is_valid(m2)]
        mp = mp[is_valid(mp)]
        mtot2 = mtot2[is_valid(mtot2)]
        inc = inc[is_valid(inc)]

        if not 'KIN' in params.keys():
            with open(outfile, 'a+') as f:
                f.write("INC(med/16th/84th)" + '\t' + str(np.median(inc)) + '\t' + str(np.percentile(inc, q=16)) + '\t' + str(np.percentile(inc, q=84)) + '\n')
        with open(outfile, 'a+') as f:
            f.write("MP(med/16th/84th)" + '\t' + str(np.median(mp)) + '\t' + str(np.percentile(mp, q=16)) + '\t' + str(np.percentile(mp, q=84)) + '\n')
        with open(outfile, 'a+') as f:
            f.write("MTOT(med/16th/84th)" + '\t' + str(np.median(mtot2)) + '\t' + str(np.percentile(mtot2, q=16)) + '\t' + str(np.percentile(mtot2, q=84)) + '\n')

        #plt.hist(mp, bins=100)
        #plt.xlabel('Pulsar mass')
        #plt.show()
        #print('Median m2: ', np.median(m2), ", inclination: ", np.median(inc), ", Pulsar mass: ", np.median(mp))
        #print('1-sigma range:', np.percentile(mp, q=16), np.percentile(mp, q=84))
        #plt.scatter(m2, inc, alpha=0.5)
        #plt.xlabel('companion mass')
        #plt.ylabel('inclination (deg)')
        #plt.show()


    if 'OMDOT' in params.keys() or 'EPS1DOT' in params.keys() or 'EPS2DOT' in params.keys():
        #print(' ')
        #print('=== Computing GR contribution ===')

        if 'ECC' in params.keys():
            ecc = params['ECC']
            ecc_err = params['ECC_ERR']
        elif 'EPS1' in params.keys():
            eps1_posterior = np.random.normal(loc=params["EPS1"],
                                           scale=params["EPS1_ERR"], size=n_samples)  # observed
            eps2_posterior = np.random.normal(loc=params["EPS2"],
                                           scale=params["EPS2_ERR"], size=n_samples)  # observed
            ecc_posterior = np.sqrt(eps1_posterior**2 + eps2_posterior**2)


        if 'EPS1DOT' in params.keys():
            eps1dot_posterior = np.random.normal(loc=params["EPS1DOT"],
                                           scale=params["EPS1DOT_ERR"], size=n_samples)  # observed
            eps2dot_posterior = np.random.normal(loc=params["EPS2DOT"],
                                           scale=params["EPS2DOT_ERR"], size=n_samples)  # observed

            omdot_posterior = (86400*365.2425*180/np.pi)*np.sqrt(eps1dot_posterior**2 + eps2dot_posterior**2)/ecc
            omdot = np.mean(omdot_posterior)
            omdot_err = np.std(omdot_posterior)
            #print("Observed OMDOT = ", omdot, " +/- ", omdot_err)
            with open(outfile, 'a+') as f:
                f.write("OMDOT" + '\t' + str(omdot) + '\t' + str(omdot_err) + '\n')

        Msun = 1.989*10**30  # kg
        Tsun = sc.G * Msun / (sc.c**3)
        Mtot = [1.6 if (np.std(mtot2)>0.2 or mtot2==0) else np.median(mtot2)]
        Mtot = Mtot[0]
        n = 2*np.pi/(params['PB']*86400)  # s
        omdot_gr = 3 * ((Tsun*Mtot)**(2/3)) * (n**(5/3)) / (1 - ecc**2)
        omdot_gr = round(omdot_gr * 86400 * 365.2425 * 180/np.pi, 5)
        #print("OMDOT_gr = {0}, for Mtot = {1}".format(omdot_gr, Mtot))
        with open(outfile, 'a+') as f:
            f.write("OMDOT_GR" + '\t' + str(omdot_gr) + ' for M_TOT = ' + str(Mtot) + '\n')

    with open(outfile, 'a+') as f:
        f.write('\n')

    print(" ")



