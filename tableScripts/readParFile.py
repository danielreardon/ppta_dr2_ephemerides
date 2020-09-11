
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Will read in a par file to a dictionary. Taken from Daniel Reardon's
script (parfile_optimiser.py in ppta_dr2_ephemeris git)
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import sys
from decimal import Decimal, InvalidOperation
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord


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






def get_derived_params(parfilename):

    """
    For a given par file, find the derived parameters in the "derived_params.txt" file.
    (searches for the par file)
    """
    derivedParamsFile = '../final/tempo2/output/derived_params.txt'

    file = open(derivedParamsFile, 'r')

    # this is only true when reading the lines we want
    take = False

    par={}

    for line in file.readlines():

        sline = line.split()

        if len(sline)>0 and take==False:
            if sline[0]==parfilename:
                take = True
        elif len(sline)==0 and take==True:
            take = False
        else: pass


        err = None
        p_type = None
        sixteenth=None
        eighty_fourth=None

        if len(sline)>1 and line[0] != "J":
            val = sline[1]
        if len(sline)>0:
            param = sline[0]


        try:
          if param != "ELAT" and param != "ELONG" and param != "PMELAT" and param != "PMELONG" and param != "MASS_FUNC" and param != "OMDOT_GR" and len(sline)>=3:
            sixteenth = sline[2]
            eighty_fourth = sline[3]

        except:
          if param == "INC_LIM(med/std)":
            err = sline[2]
          elif param == "ECC(med/std)":
            err = sline[2]
          elif param == "OM(med/std)":
            err = sline[2]
          elif param == "T0(med/std)":
            err = sline[2]



        if take==True and line[0]!="J":

            # fill par with the parameter values 
            par[param] = val
            if err:
                par[param+"_ERR"] = float(err)
            if p_type:
                par[param+"_TYPE"] = p_type
            if sixteenth:
                par[param+"_16th"] = float(sixteenth)
            if eighty_fourth:
                par[param+"_84th"] = float(eighty_fourth)


    return par




