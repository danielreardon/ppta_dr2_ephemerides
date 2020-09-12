#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 16:57:06 2020

@author: dreardon
"""

from astropy import units as u
import numpy as np
from numpy import cos, sin, tan, pi, arcsin, arccos, arctan2
from astropy.coordinates import SkyCoord, ICRS, BarycentricTrueEcliptic

rad_to_mas = 180*3600*1000/pi
parsec_to_m = 3.08567758e+16
sec_per_year = 86400*365.2425

# Equatorial
ra = '10:00:51.9032239'
dec = '-17:53:49.38722'
pmra = 1.0 / rad_to_mas / sec_per_year
pmdec = 1.0 / rad_to_mas / sec_per_year
c = SkyCoord(ra + ' ' + dec, unit=(u.hourangle, u.deg))
pm = ICRS(ra=c.ra.deg*u.degree, dec=c.dec.deg*u.deg,
          pm_ra_cosdec=pmra*u.rad/u.s,
          pm_dec=pmdec*u.rad/u.s)
pm_ecliptic = pm.transform_to(BarycentricTrueEcliptic)

# Ecliptic
elat = c.barycentrictrueecliptic.lat.value
elong = c.barycentrictrueecliptic.lon.value
pmelat = pm_ecliptic.pm_lat.value
pmelong = pm_ecliptic.pm_lon_coslat.value
dec = c.dec.value * pi /180
ra = c.ra.value * pi/180
eps = 0.409092804223

print(elat, elong, pmelat, pmelong)

# Now, convert back!
elat = elat*pi/180
elong = elong*pi/180
pmelat = pmelat / rad_to_mas / sec_per_year
pmelong = pmelong / rad_to_mas / sec_per_year

dec_new = arcsin(cos(eps)*sin(elat) + sin(eps)*cos(elat)*sin(elong))
pmdec_new = \
    (1/cos(dec))*((cos(eps)*cos(elat) - sin(eps)*sin(elat)*sin(elong))*pmelat + \
    sin(eps)*cos(elat)*cos(elong)*(pmelong / cos(elat)))

print(dec_new*180/pi, pmdec_new * rad_to_mas * sec_per_year)

ra_new = arctan2( cos(eps)*sin(elong) - sin(eps)*tan(elat), cos(elong) )
pmra_new = cos(ra)**2 * ( ( (sin(elong)/cos(elong)**2)*(cos(eps)*sin(elong) - sin(eps)*tan(elat)) + cos(eps)) * (pmelong / cos(elat)) - \
           sin(eps)*pmelat/(cos(elong)*cos(elat)**2)) * cos(dec)

print(ra_new*24/pi, pmra_new * rad_to_mas * sec_per_year)




c1=cos(elat)*cos(eps)-sin(elat)*sin(elong)*sin(eps)
c2=cos(elat)*cos(elong)*sin(eps)
#c2=sin(elong)*sin(eps)
cdelta = np.sqrt(c1*c1 + c2*c2)
c1/=cdelta
c2/=cdelta

#pmelong /= cos(elat)

pmra_new2 = c1*pmelong -c2*pmelat
pmdec_new2 = c2*pmelong + c1*pmelat

print(pmra_new2  * rad_to_mas * sec_per_year, pmdec_new2 * rad_to_mas * sec_per_year)


ra = np.random.rand(20)*pi
dec = np.random.rand(20)*pi

pos_array = []
# Set up unit vectors
for ira in ra:
    for idec in dec:
        pos_array.append([np.cos(ira) * np.cos(idec), np.sin(ira) * np.cos(idec), np.sin(idec)])

diff_array = []
hd_array = []
dipole_array = []
monopole_array = []
for pos1 in pos_array:
    for pos2 in pos_array:
        if np.all(pos1 == pos2):
            hd = 1
            dipole_array.append(1 + 1e-5)
        else:
            omc2 = (1 - np.dot(pos1, pos2)) / 2
            hd = 1.5 * omc2 * np.log(omc2) - 0.25 * omc2 + 0.5
            dipole_array.append(np.dot(pos1, pos2))
        hd_array.append(hd)
        diff_array.append( np.arccos(np.dot(pos1, pos2) / np.sqrt(np.dot(pos1, pos1)) / np.sqrt(np.dot(pos2, pos2))))

import matplotlib.pyplot as plt
plt.scatter(diff_array, hd_array)
plt.scatter(diff_array, dipole_array)
plt.show()


