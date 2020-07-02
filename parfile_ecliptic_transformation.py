#!/usr/bin/python

"""
This scripts takes a .par file and 
- replace the equatorial coordinates RA and DEC with ecliptic coordinates ELONG and ELAT
- replace PMRA  and PMDEC (the proper motion in the equatorial coordinate system) with 
  PMLONG and PMLAT (the proper motion in the ecliptic coordinate system)

usage: python *.py *par

author: Mohsen Shamohammadi
email: msh.ph.ir@gmail.com
"""

from __future__ import print_function
import numpy as np
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS, BarycentricTrueEcliptic
from astropy import _erfa as erfa
from astropy.coordinates.funcs import get_jd12
from astropy.coordinates.builtin_frames.ecliptic import EQUINOX_J2000
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product


def _true_ecliptic_rotation_matrix(equinox):
    """
    # This code calls pnm06a from ERFA, which retrieves the precession
    # matrix (including frame bias) according to the IAU 2006 model, and
    # including the nutation. This family of systems is less popular
    # (see https://github.com/astropy/astropy/pull/6508).

    # This function is originated from Astropy package:
    astropy.coordinates.builtin_frames.ecliptic_transforms._true_ecliptic_rotation_matrix(equinox)
    
    You can find the equinox for other epochs using Time method in astropy.time
    Example: the equinox of J1950 is:

        from astropy.time import Time
        equinox = Time('J1950')
    """
    jd1, jd2 = get_jd12(equinox, 'tt')
    rnpb = erfa.pnm06a(jd1, jd2)
    _, nut_obl = erfa.nut06a(jd1, jd2)*u.radian
    obl = erfa.obl06(jd1, jd2)*u.radian + nut_obl  # calculate the true obliquity of the ecliptic
    return matrix_product(rotation_matrix(obl, 'x'), rnpb)


def equatorial_to_ecliptic(ra, ra_err, dec, dec_err, pmra, pmdec):
    #convert Equatorial coordinates (RAJ, DECJ) to Barycentric true ecliptic coordinates (ELONG, ELAT):
    c = SkyCoord(str(ra) + ' ' + str(dec), unit=(u.hourangle, u.deg))
    ra = c.ra.value
    dec = c.dec.value
    elat = c.transform_to(BarycentricTrueEcliptic).lat.value
    elong = c.transform_to(BarycentricTrueEcliptic).lon.value

    #convert the units of equatorial coordinate uncertainties (RAJ_ERR, DECJ_ERR) to degree:
    c_err = SkyCoord(str(ra_err) + ' ' + str(dec_err), unit=(u.hourangle, u.deg))
    ra_err = c_err.ra.value
    dec_err = c_err.dec.value

    #proper motion transformation:
    pm = ICRS(ra=c.ra.deg*u.degree, dec=c.dec.deg*u.deg,
              pm_ra_cosdec = np.float64(pmra)*u.mas/u.yr,
              pm_dec = np.float64(pmdec)*u.mas/u.yr)
    pm_ecliptic = pm.transform_to(BarycentricTrueEcliptic)
    pmlat = pm_ecliptic.pm_lat.value
    pmlong = pm_ecliptic.pm_lon_coslat.value
    return ra, ra_err, dec, dec_err, elong, elat, pmlong, pmlat



def equatorial_to_ecliptic_uncertainty(ra, ra_err, dec, dec_err, elong, elat, A):
    """
    Convert Equatorial coordinate uncertainties to Ecliptic using Euler rotations
    Ecliptic Coords = A * Eqiatorial Coords
    Error[Eqs] = Error[([A00, A01, A02], [A10, A11, A12], [A20, A21, A22]) * ECs]
    """

    #convert the coordinates units from degree to radian:
    [ra, dec, elong, elat] = [item * np.pi/180. for item in [ra, dec, elong, elat]]

    #Error analysis: 
    #EqX, EqY, EqZ = Error(A * Eqiatorial Coords) 
    EqX = np.sqrt(A[0,0] ** 2 * np.cos(dec) ** 2 * ra_err ** 2 + (A[0,0] * np.cos(ra) + A[0,1] * np.sin(ra)) ** 2 * dec_err ** 2 + A[0,1] ** 2 * np.cos(dec) ** 2 * ra_err ** 2 + A[0,2] ** 2 * dec_err ** 2)
    EqY = np.sqrt(A[1,0] ** 2 * np.cos(dec) ** 2 * ra_err ** 2 + (A[1,0] * np.cos(ra) + A[1,1] * np.sin(ra)) ** 2 * dec_err ** 2 + A[1,1] ** 2 * np.cos(dec) ** 2 * ra_err ** 2 + A[1,2] ** 2 * dec_err ** 2)
    EqZ = np.sqrt(A[2,0] ** 2 * np.cos(dec) ** 2 * ra_err ** 2 + (A[2,0] * np.cos(ra) + A[2,1] * np.sin(ra)) ** 2 * dec_err ** 2 + A[2,1] ** 2 * np.cos(dec) ** 2 * ra_err ** 2 + A[2,2] ** 2 * dec_err ** 2)


    ### calculating uncertainties in Ecliptic transformation: 
    # ELAT_ERR:
    ELAT_ERR = EqZ
    
    #ELONG_ERR from EcY:
    q2 = EqX ** 2 - (ELAT_ERR * np.cos(elong)) ** 2
    q = q2 / np.cos(elat) ** 2
    ELONG_ERR_1 = np.sqrt(q)

    #ELONG_ERR from EcY:
    q2 = EqY ** 2 - (ELAT_ERR * np.sin(elong)) ** 2
    q = q2 / np.cos(elat) ** 2
    ELONG_ERR_2 = np.sqrt(q)

    #rms error of ELONG_ERR_{1,2}:
    ELONG_ERR = np.sqrt(np.mean(np.square([ELONG_ERR_1, ELONG_ERR_2])))
    return ELONG_ERR, ELAT_ERR



def readit(filename):
    #reading par file:
    with open(filename, 'r') as f:
        data = f.readlines()
        for line in data:
            if 'RAJ ' in line:
                sline = line.split()
                RAJ = sline[1]
                RAJ_ERR = "0:0:" + sline[3]
            if 'DECJ ' in line:
                sline = line.split()
                DECJ = sline[1]
                DECJ_ERR = "0:0:" + sline[3]
            if 'PMRA ' in line:
                sline = line.split()
                PMRA = sline[1]
            if 'PMDEC ' in line:
                sline = line.split()
                PMDEC = sline[1]
    return data, RAJ, RAJ_ERR, DECJ, DECJ_ERR, PMRA, PMDEC


def writeit(filename, mat):
    #reading a par file and get coordinates and proper motion:
    data, raj, raj_e, decj, decj_e, pmra, pmdec= readit(filename)


    #finding float format of equatorial parameters and calculating ecliptic parameters:
    RAJ, RAJ_ERR, DECJ, DECJ_ERR, ELONG, ELAT, PMLONG, PMLAT = equatorial_to_ecliptic(raj, raj_e, decj, decj_e, pmra, pmdec)

    #finding ecliptic uncertainties:
    ELONG_ERR, ELAT_ERR = equatorial_to_ecliptic_uncertainty(RAJ, RAJ_ERR, DECJ, DECJ_ERR, ELONG, ELAT, mat)

    #writing a new par file:
    newpar = filename.replace('.par','')
    with open('{}_ecliptic.par'.format(newpar), 'w') as g:
        k = ""
        for line in data:
            if 'RAJ ' in line:
                line = 'ELONG          {}           1  {}\n'.format(ELONG, ELONG_ERR)
            if 'DECJ ' in line:
                line = 'ELAT           {}           1  {}\n'.format(ELAT, ELAT_ERR)
            if 'PMRA ' in line:
                line = 'PMELONG        {}  1 \n'.format(PMLONG)
            if 'PMDEC ' in line:
                line = 'PMELAT         {}  1 \n'.format(PMLAT)
            k += "{}".format(line)
        g.write(k)
    print('{}_ecliptic.par is written'.format(newpar))



if __name__ == '__main__':

    #finding the true ecliptic rotation matrix using EQUINOX_J2000 (the instance of equinox):
    mat = _true_ecliptic_rotation_matrix(EQUINOX_J2000)

    for df in sys.argv[1:]:
        writeit(df, mat)
