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

datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/partim/dr2_boris/new_params/'
outdir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/partim/dr2_boris/output/'
parfiles = sorted(glob.glob(datadir + 'J*.par'))

for par in parfiles:
    if 'J0437' in par:
        continue
    print("Working on: ", par.split('/')[-1].split('.')[0])
    params = read_par(par)

    # Check if pulsar needs ELL1:
    if 'ECC' in params.keys():
        if params['A1']*params['ECC']**2 < \
           0.01*params['TRES']*10**-6/np.sqrt(params['NTOA']):
            print(params['PSRJ'], " Needs ELL1 model")
            # make ELL1
            make_ell1(par, outdir + params['PSRJ'] + '_ell1.par')
    elif 'EPS1' in params.keys():
        ecc = params['EPS1']**2 + params['EPS2']**2
        if params['A1']*ecc**2 > \
           0.01*params['TRES']*10**-6/np.sqrt(params['NTOA']):
            print('WARNING! Pulsar should NOT have ELL1 model')

    if 'XDOT' in params.keys() and 'OMDOT' in params.keys():





