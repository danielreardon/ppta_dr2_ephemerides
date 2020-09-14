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

"""
Definitions
"""

import os
from matplotlib import rc
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin/'
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rcParams.update({'font.size': 20})

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

def derive_combined_mass(params, plot=False, n_samples=10000000):

    if 'PSRJ' in params:
        psrname = params['PSRJ']
    elif 'PSR' in params:
        psrname = params['PSR']

    # Define other useful constants
    M_sun = 1.98847542e+30  # kg
    Tsun = 4.926790353700459e-06  #s
    rad_to_mas = 180*3600*1000/np.pi
    parsec_to_m = 3.08567758e+16
    sec_per_year = 86400*365.2425

    mass_func = (4*np.pi**2/sc.G) * (params['A1']*sc.c)**3/(params['PB']*86400)**2
    mass_func = mass_func/M_sun

    Omega = np.random.rand(n_samples)*2*np.pi
    cosi = np.random.rand(n_samples)
    i_prior = np.arccos(cosi)*180/np.pi
    i_prior[i_prior<0] = -i_prior[i_prior<0]

    #i_prior = np.random.rand(n_samples)*90

    # XDOT
    if 'XDOT' in params.keys():
        xdot = np.random.normal(loc=params["XDOT"], scale=params["XDOT_ERR"], size=n_samples)
        if np.abs(params['XDOT']) > 1e-10:
            xdot *= 1e-12
        a1 = np.random.normal(loc=params["A1"], scale=params["A1_ERR"], size=n_samples)
        pmelat = np.random.normal(loc=params["PMELAT"], scale=params["PMELAT_ERR"], size=n_samples)
        pmelong = np.random.normal(loc=params["PMELONG"], scale=params["PMELONG_ERR"], size=n_samples)
        pm_tot = np.sqrt(pmelat**2 + pmelong**2)
        pm = pm_tot/(sec_per_year*rad_to_mas)

        A = (-pmelong/(sec_per_year*rad_to_mas))*np.sin(Omega) + (-pmelat/(sec_per_year*rad_to_mas))*np.cos(Omega)
        i = np.arctan2(a1 * A, xdot)*180/np.pi
        i[i>90] = i[i>90] - 180
        i[i<0] = -i[i<0]
        i[i>90] = i[i>90] - 180
        i[i<0] = -i[i<0]

    elif 'KIN' in params.keys():
        i = np.random.normal(loc=params["KIN"], scale=params["KIN_ERR"], size=n_samples)
        i[i>90] = i[i>90] - 180
        i[i<0] = -i[i<0]
        i[i>90] = i[i>90] - 180
        i[i<0] = -i[i<0]
    #plt.figure()
    #plt.hist(i_prior, bins=1000, alpha=0.5)
    #plt.show()

    #plt.figure()
    #plt.hist(i, bins=1000, alpha=0.5)
    #plt.show()

    if 'H3' in params.keys():
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
        if '1017' in psrname:
            #sini_lim = 0.7880737788017274
            sini_lim = 1
        elif '1022' in psrname:
            #sini_lim = 0.96928029962
            sini_lim = 1
        else:
            sini_lim = 1
        cut = np.argwhere((m2 > mass_func) * (m2 < 1.4) * (sini < sini_lim))
        m2 = m2[cut]
        sini = sini[cut]
        i = i[cut]
        i_prior = i_prior[cut]
        inc = np.arcsin(sini)*180/np.pi
        mtot2 = (m2 * sini)**3 / mass_func
        mp = np.sqrt(mtot2[is_valid(mtot2)*is_valid(m2)]) - m2[is_valid(mtot2)*is_valid(m2)]
        mp = mp[is_valid(mp)]
        mtot2 = mtot2[is_valid(mtot2)]
        mtot = np.sqrt(mtot2[(mtot2 > mass_func)])
        cut = is_valid(inc).squeeze()
        inc = inc[cut]
        i = i[cut]
        i_prior = i_prior[cut]
        shap_label = r'Shapiro delay, $h_3$ and $\zeta$'
    else:
        sini_lim = 1
        m2 = np.random.normal(loc=params['M2'], scale=params['M2_ERR'], size=n_samples)
        if 'SINI_ERR' in params.keys():
            sini = np.random.normal(loc=params['SINI'], scale=params['SINI_ERR'], size=n_samples)
        else:
            sini = np.sin(i*np.pi/180)
        cut = np.argwhere((m2 > mass_func) * (m2 < 1.4) * (sini < sini_lim))
        m2 = m2[cut]
        sini = sini[cut]
        i = i[cut]
        i_prior = i_prior[cut]
        inc = np.arcsin(sini)*180/np.pi
        mtot2 = (m2 * sini)**3 / mass_func
        mp = np.sqrt(mtot2[is_valid(mtot2)*is_valid(m2)]) - m2[is_valid(mtot2)*is_valid(m2)]
        mp = mp[is_valid(mp)]
        mtot2 = mtot2[is_valid(mtot2)]
        mtot = np.sqrt(mtot2[(mtot2 > mass_func)])
        cut = is_valid(inc).squeeze()
        inc = inc[cut]
        i = i[cut]
        i_prior = i_prior[cut]
        shap_label = r'Shapiro delay, $M_c$ and $\sin{i}$'


    #plt.figure()
    #plt.hist(inc, bins=1000, alpha=0.5)
    #plt.xlim(0, 90)
    #plt.show()


    # Derive Mtot from OMDOT
    if 'OMDOT' in params.keys():
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
    else:
        Mtot = 5*np.random.rand(n_samples)  # uniform

    i = np.sort(i)
    i_prior = np.sort(i_prior)
    inc = np.sort(inc)

    bins = 1000

    plt.figure()
    n, b, p = plt.hist(i_prior, bins=bins, alpha=0.5, range=(0, 90))
    n2, b2, p = plt.hist(i, bins=bins, alpha=0.5, range=(0, 90))
    n3, b3, p = plt.hist(inc, bins=bins, alpha=0.5, range=(0, 90))
    plt.close()

    combined = n*n2*n3
    mx = np.max(combined)

    i_posterior = np.random.choice(np.linspace(0,90,bins), size=(1, n_samples), p=(combined/np.sum(combined))).squeeze()

    sini = np.sin(i_posterior*np.pi/180).squeeze()
    m2 = m2.squeeze()
    mtot2 = np.power(np.multiply(m2.squeeze(), sini[0:np.size(m2)]), 3) / mass_func
    mtot2 = mtot2.squeeze()
    mtot = np.sqrt(mtot2)
    mtot = mtot[is_valid(mtot)]

    plt.figure()
    n4, b, p = plt.hist(Mtot, bins=bins, alpha=0.5, range=(0, 5))
    n5, b2, p = plt.hist(mtot, bins=bins, alpha=0.5, range=(0, 5))
    plt.close()




    plt.figure(figsize=(18, 6))
    plt.subplot(1, 2, 1)
    x = np.linspace(0,90, bins)
    dx = x[1] - x[0]
    scalefactor = 1 / dx
    plt.plot(x, scalefactor*savgol_filter(n2/np.sum(n2), 51, 3), color='mediumblue', linewidth=2)
    #plt.fill_between(np.linspace(0,90,bins), savgol_filter(n2/np.sum(n2), 51, 3), color='crimson', alpha=0.3)
    scalefactor = 1 / dx
    plt.plot(x, scalefactor*savgol_filter(n3/np.sum(n3), 51, 3), color='crimson', linewidth=2, alpha=0.8)
    #plt.fill_between(np.linspace(0,90,bins), savgol_filter(n3/np.sum(n3), 51, 3), color='darkorange', alpha=0.3)
    scalefactor = 1 / dx
    plt.plot(x, scalefactor*savgol_filter(n/np.sum(n), 51, 3), color='dimgrey', linestyle='--', linewidth=2)
    #plt.fill_between(np.linspace(0,90,bins), savgol_filter(n/np.sum(n), 51, 3), color='mediumblue', alpha=0.3)


    #plt.plot(np.linspace(0,90,bins), savgol_filter(combined/np.sum(combined), 51, 3), color='k', linewidth=2)
    scalefactor = 1 / dx
    plt.fill_between(x, scalefactor*savgol_filter(combined/np.sum(combined), 51, 3), color='k', alpha=0.15)
    #plt.plot(np.linspace(0,90,bins), n/mx*np.max(n))
    #plt.plot(np.linspace(0,90,bins), n2/mx*np.max(n2))
    #plt.plot(np.linspace(0,90,bins), n3/mx*np.max(n3))
    if 'XDOT' in params.keys():
        plt.legend([r'Kopeikin, $\dot{x}$', shap_label,r'Uniform $\cos{i}$ prior',  r'Total'], loc='upper left')
    else:
        plt.legend([r'Kopeikin, $i$', shap_label,r'Uniform $\cos(i)$ prior', r'Total'] )#, loc='upper left')
    plt.xlim((0, 90))
    yl = plt.ylim()
    plt.ylim((0, yl[1]))
    plt.xlabel(r'Inclination, $i$ ($^\circ$)')
    plt.ylabel(r'Probability')
    #locs, labels = plt.yticks()            # Get locations and labels
    #labels[labels == '0.000'] = '0'
    #plt.yticks(locs, labels)  # Set locations and labels


    combined = n4*n5
    plt.subplot(1, 2, 2)
    x = np.linspace(0,5, bins)
    dx = x[1] - x[0]
    scalefactor = 1 / dx
    plt.plot(x, scalefactor*savgol_filter(n4/np.sum(n4), 51, 3), color='midnightblue', linewidth=2)
    #plt.fill_between(np.linspace(0,5, bins), savgol_filter(n4/np.sum(n4), 51, 3), color='crimson', alpha=0.3)
    scalefactor = 1 / dx
    plt.plot(x, scalefactor*savgol_filter(n5/np.sum(n5), 51, 3), color='darkmagenta', linewidth=2, alpha=0.9)
    #plt.fill_between(np.linspace(0,5, bins), savgol_filter(n5/np.sum(n5), 51, 3), color='darkorange', alpha=0.3)
    #plt.plot(np.linspace(0,5, bins), savgol_filter(combined/np.sum(combined), 51, 3), color='k', linewidth=2)
    scalefactor = 1 / dx
    plt.fill_between(x, scalefactor*savgol_filter(combined/np.sum(combined), 51, 3), color='k', alpha=0.15)
    plt.xlim((0, 5))
    yl = plt.ylim()
    if '1600' in psrname:
        plt.ylim((0, 1))
    else:
        plt.ylim((0, yl[1]))
    plt.xlabel(r'Total system mass, $M_{\rm tot}$ ($M_\odot$)')
    plt.ylabel(r'Probability')
    #locs, labels = plt.yticks()            # Get locations and labels
    #plt.yticks(locs, [])  # Set locations and labels
    plt.legend(['Periastron advance $\dot{\omega}$', 'Inclination, $i$', 'Total'] )#, loc='upper left')
    plt.tight_layout()
    if plot:
        plt.savefig('{0}_mass.pdf'.format(psrname), bbox_inches='tight', pad_inches=0.1,)
        plt.show()
    else:
        plt.close()

    mtot_posterior = np.random.choice(np.linspace(0,5,bins), size=(1, n_samples), p=(combined/np.sum(combined))).squeeze()
    mtot_posterior = mtot_posterior[(mtot_posterior > mass_func)]

    mp_new = mtot_posterior[0:np.size(m2)] - m2
    mp_new = mp_new[is_valid(mp_new)]


#    plt.figure()
#    plt.hist(mp_new, bins=100)
#    plt.xlim((0, 4))
#    plt.xlabel('Pulsar mass')
#    if plot:
#        plt.show()
#    else:
#        plt.close()


    mp_med = np.percentile(mp_new, q=50)
    high_error = np.percentile(mp_new, q=84) - np.percentile(mp_new, q=50)
    lower_error = np.percentile(mp_new, q=50) - np.percentile(mp_new, q=16)

    #return mp_med, high_error,lower_error
    mp = mp_new
    mtot = mtot_posterior
    inc = i_posterior

    # print(np.percentile(mp, q=5.0))

    return {'M2': np.median(m2), 'M2_std': np.std(m2), 'M2_lolim': np.percentile(m2, q=16.0), 'M2_uplim': np.percentile(m2, q=84.0), 'Mpsr': np.median(mp), 'Mpsr_std': np.std(mp), 'Mpsr_lolim': np.percentile(mp, q=16.0), 'Mpsr_uplim': np.percentile(mp, q=84.0), 'Mtot': np.median(mtot), 'Mtot_std': np.std(mtot), 'Mtot_lolim': np.percentile(mtot, q=16.0), 'Mtot_uplim': np.percentile(mtot, q=84.0), 'inc': np.median(inc), 'inc_std': np.std(inc), 'inc_lolim': np.percentile(inc, q=16.0), 'inc_uplim': np.percentile(inc, q=84.0)}




