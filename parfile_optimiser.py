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
import sys
from decimal import Decimal, InvalidOperation
import libstempo as T
import GalDynPsr
import scipy.constants as sc
import matplotlib.pyplot as plt
import scipy.stats as stats
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS, BarycentricTrueEcliptic



"""
Definitions
"""


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



"""
Start of code
"""

datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/partim/dr2_boris/new_params_ver1/'
outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/partim/dr2_boris/new_params_ver1/output/'
parfiles = sorted(glob.glob(datadir + 'J*.par'))

outfile = outdir + 'derived_params.txt'

n_samples = 10000
# Define other useful constants
M_sun = 1.98847542e+30  # kg
rad_to_mas = 180*3600*1000/np.pi
parsec_to_m = 3.08567758e+16
sec_per_year = 86400*365.2425

for par in parfiles:
    if 'J0437' in par:
        continue
    print(par.split('/')[-1].split('.')[0])
    params = read_par(par)

    with open(outfile, 'a+') as f:
        f.write(par.split('/')[-1].split('.')[0] + '\n')


    #Convert to Galactic
    c = SkyCoord(params['RAJ'] + ' ' + params['DECJ'],
                 unit=(u.hourangle, u.deg))
    ldeg = c.galactic.l.value
    sigl = 0
    bdeg = c.galactic.b.value
    sigb = 0

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
    if 'XDOT' in params.keys() and 'OMDOT' in params.keys():
        print('Check for Kopeikin terms')


    if 'PX' in params.keys():
        D_prior = 1/np.random.normal(loc=params["PX"],
                                       scale=params["PX_ERR"], size=n_samples) # observed
        #plt.hist(D_prior, bins=100)
        #plt.show()
        dkpc = np.median(D_prior[(D_prior > 0)*(D_prior < 100)])
        sigd = np.std(D_prior[(D_prior > 0)*(D_prior < 100)])
        #sigd = dkpc*params["PX_ERR"]/params["PX"]
        #print("Parallax distance (kpc) = ", round(dkpc, 3), " +/- ", round(sigd, 3))
        with open(outfile, 'a+') as f:
            f.write("D_PX" + '\t' + str(dkpc) + '\t' + str(sigd) + '\n')


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
        pmra = params['PMRA']
        pmdec = params['PMDEC']

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
        pm = np.sqrt(pmra**2 + pmdec**2)/(sec_per_year*rad_to_mas)
        D_posterior = sc.c*pbdot_shklovskii/(pm**2*pb*86400)/parsec_to_m/1000
        D = np.mean(D_posterior)
        D_err = np.std(D_posterior)
        #print("Shklovskii distance (kpc) = ", round(D, 3), " +/- ", round(D_err,3))
        with open(outfile, 'a+') as f:
            f.write("D_SHK" + '\t' + str(D) + '\t' + str(D_err) + '\n')
        if 'PX' in params.keys():
            Davg = np.average([dkpc, D], weights=[1/sigd, 1/D_err])
            Davg_err = 1
            #print("Combined distance (kpc) = ", round(Davg, 3), " +/- ", round(Davg_err,3))

    if 'F2' in params.keys():
        #print(' ')
        #rint('=== Computing radial velocity from F2 ===')
        f0 = params['F0']
        try:
            f0err = params['F0_ERR']
        except KeyError:
            continue
        f1 = params['F1']
        f1_err = params['F1_ERR']
        f2 = params['F2']
        f2_err = params['F2_ERR']
        pmra = params['PMRA']
        pmdec = params['PMDEC']
        pm = np.sqrt(pmra**2 + pmdec**2)/(sec_per_year*rad_to_mas)
        p0_obs = 1/f0
        p1_obs = -f1/f0**2
        p2_obs = f1*(2*f1/f0**3) - (1/f0**2)*(f2)
        Ex_pl =  GalDynPsr.modelLb.Expl(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term parallel to the Galactic plane
        Ex_z =  GalDynPsr.modelLb.Exz(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term perpendicular to the Galactic plane

        dv_r = np.Inf
        v_r = 0
        # iterate Doppler correction
        while dv_r > np.abs(0.01*v_r):
            p0_int = p0_obs/(1 - v_r*1000/sc.c) # approximate Doppler correction
            p1_int = p1_obs - GalDynPsr.pdotint.PdotGal(Ex_pl, Ex_z, p0_int) # Shklovskii correction
            v_r_new =  (2*p1_int*D*parsec_to_m - p2_obs*(sc.c/pm**2))/(3*p0_int)/1000 # km/s ... assuming p2_int=0
            dv_r = np.abs(v_r - v_r_new)
            v_r = v_r_new
        #print("Observed P2 = ", p2_obs)
        #print("Radial velocity = ", round(v_r/1000,1), " +/- ", round(f2_err/f2 * v_r/1000, 1), " km/s")
        with open(outfile, 'a+') as f:
            f.write("V_R" + '\t' + str(round(v_r/1000,1)) + '\t' + str(round(f2_err/f2 * v_r/1000, 1)) + '\n')

    if 'XDOT' in params.keys():
        #print(' ')
        #print('=== Computing limit on i from XDOT ===')
        i_limit =  np.abs(180/np.pi * np.arctan(params['A1'] * np.sqrt(params['PMRA']**2 + params['PMDEC']**2)/(rad_to_mas*sec_per_year) / params['XDOT']))
        #print("i <= ", int(np.ceil(i_limit)))
        with open(outfile, 'a+') as f:
            f.write("INC" + '\t' + "<" + str(i_limit)+ '\n')

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
        Mtot = 1.6
        n = 2*np.pi/(params['PB']*86400)  # s
        omdot_gr = 3 * ((Tsun*Mtot)**(2/3)) * (n**(5/3)) / (1 - ecc**2)
        omdot_gr = round(omdot_gr * 86400 * 365.2425 * 180/np.pi, 5)
        #print("OMDOT_gr = {0}, for Mtot = {1}".format(omdot_gr, Mtot))
        with open(outfile, 'a+') as f:
            f.write("OMDOT_GR" + '\t' + str(omdot_gr) + ' for M_TOT = ' + str(Mtot) + '\n')

    with open(outfile, 'a+') as f:
        f.write('\n')

    print(" ")

