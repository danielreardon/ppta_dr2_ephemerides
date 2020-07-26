import glob
import sys, os
from decimal import Decimal, InvalidOperation
import GalDynPsr
import matplotlib.pyplot as PLT
import numpy as NP
import scipy.constants as FCNST
import matplotlib
import scipy.stats as stats
from astropy import units as U
from scipy.signal import savgol_filter
from astropy.coordinates import SkyCoord, ICRS, BarycentricTrueEcliptic

n_samples = 10000
# Define other useful constants
M_sun = 1.98847542e+30  # kg
Tsun = 4.926790353700459e-06  #s
rad_to_mas = 180*3600*1000/NP.pi
parsec_to_m = 3.08567758e+16
sec_per_year = 86400*365.2425

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

def distance_from_parallax(psrparms):
    if 'PX' in psrparms:
        if 'PX_ERR' in psrparms:
            d_prior = 1/NP.random.normal(loc=psrparms["PX"], scale=psrparms["PX_ERR"], size=n_samples) # observed
            d_array = d_prior[(d_prior > 0)*(d_prior < 100)]
            dkpc = NP.median(d_array)
            sigd = NP.std(d_array)
    
            return {'dpsr': dkpc, 'dpsr_std': sigd, 'dpsr_lolim': NP.percentile(d_array, q=16.0), 'dpsr_uplim': NP.percentile(d_array, q=84.0)}
        else:
            return -1
    return -1

def distance_from_pbdot(psrparms):

    if ('ELAT' in psrparms) and ('ELONG' in psrparms):
        c = SkyCoord(psrparms['ELONG'], psrparms['ELAT'], frame=BarycentricTrueEcliptic, unit=(U.deg, U.deg))
    else:
        c = SkyCoord(psrparms['RAJ'] + ' ' + psrparms['DECJ'],
                 unit=(U.hourangle, U.deg))

    #Convert to Galactic

    ldeg = c.galactic.l.value
    sigl = 0.0
    bdeg = c.galactic.b.value
    sigb = 0.0

    if not 'ecliptic' in datadir:

        #Convert to Ecliptic
        #ecliptic transformation:

        elat = c.barycentrictrueecliptic.lat.value
        elon = c.barycentrictrueecliptic.lon.value

        #proper motion transformation:

        pm = ICRS(ra=c.ra.deg*U.degree, dec=c.dec.deg*U.deg, pm_ra_cosdec = psrparms['PMRA']*U.mas/U.yr, pm_dec = psrparms['PMDEC']*U.mas/U.yr)
        pm_ecliptic = pm.transform_to(BarycentricTrueEcliptic)
        pmlat = pm_ecliptic.pm_lat.value
        pmlon = pm_ecliptic.pm_lon_coslat.value
    
    if 'PBDOT' in psrparms:
        pbdot_posterior = NP.random.normal(loc=psrparms["PBDOT"], scale=psrparms["PBDOT_ERR"], size=n_samples)  # observed
        pbdot_grav = 0.0
        if 'PX' in psrparms:
            dkpc, sigd, dkpc16, dkpc84 = distance_from_parallax(psrparms)
        else:
            dkpc = 1.0
            sigd = 0.2*dkpc
        pb = psrparms['PB']
        pb_err = psrparms['PB_ERR']
        if ('PMELAT' in psrparms) and ('PMELONG' in psrparms):
        # if 'ecliptic' in datadir:
            pmelat = psrparms['PMELAT']
            pmelong = psrparms['PMELONG']
            pm_tot = NP.sqrt(pmelat**2 + pmelong**2)
            pm = pm_tot/(sec_per_year*rad_to_mas)
        else:
            pmra = psrparms['PMRA']
            pmdec = psrparms['PMDEC']
            pm_tot = NP.sqrt(pmelat**2 + pmelong**2)
            pm = pm_tot/(sec_per_year*rad_to_mas)
    
        # Expected Shklovskii and Galactic terms

        Ex_pl =  GalDyNP.r.modelLb.Expl(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term parallel to the Galactic plane
        Ex_z =  GalDynPsr.modelLb.Exz(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term perpendicular to the Galactic plane
        errpl = NP.abs(0.03*Ex_pl) # 3% for circular Bovy et al. (2012)
        errz = NP.abs(0.1*Ex_z) # 10% for vertical Holmberg and Flynn (2004)
        pbdot_gal = GalDynPsr.pdotint.PdotGal(Ex_pl, Ex_z, pb*86400)
        pbdot_gal_err = GalDynPsr.pdotint.ErrPdotGal(Ex_pl, errpl, Ex_z, errz, pb*86400, pb_err*86400)

        pbdot_gal_posterior = NP.random.normal(loc=pbdot_gal, scale=pbdot_gal_err, size=n_samples)  # sample randomly from pbdot_gal distribution
        pbdot_shklovskii = NP.subtract(pbdot_posterior, pbdot_gal_posterior) - pbdot_grav

        d_posterior = FCNST.c*pbdot_shklovskii/(pm**2*pb*86400)/parsec_to_m/1000
        d = NP.mean(d_posterior)
        d_err = NP.std(d_posterior)

        return {'dpsr': NP.median(d_posterior), 'dpsr_std': NP.std(d_posterior), 'dpsr_lolim': NP.percentile(d_posterior, q=16.0), 'dpsr_uplim': NP.percentile(d_posterior, q=84.0)}
    else:
        return -1
    return -1

def mass_from_psrparms(psrparms):
    mtot = 0
    if 'H3' in psrparms:
        h3 = psrparms["H3"]*10**6
        h3_err = psrparms["H3_ERR"]*10**6
        h3_arr = NP.random.normal(loc=h3, scale=h3_err, size=n_samples)
        if 'H4' in psrparms:
            h4 = psrparms["H4"]*10**6
            h4_err = psrparms["H4_ERR"]*10**6
            h4_arr = NP.random.normal(loc=h4, scale=h4_err, size=n_samples)
        elif 'STIG' in psrparms:
            stig = psrparms["STIG"]
            stig_err = psrparms["STIG_ERR"]
            stig_arr = NP.random.normal(loc=stig, scale=stig_err, size=n_samples)
            h4_arr = stig_arr * h3_arr
        sini = 2 * h3_arr * h4_arr / ( h3_arr**2 + h4_arr**2 )
        m2 = h3_arr**4 / h4_arr**3  # in us
        m2 = m2/(Tsun*10**6)
    elif 'M2' in psrparms:
        m2 = NP.random.normal(loc=psrparms["M2"], scale=psrparms["M2_ERR"], size=n_samples)
        if 'KIN' in psrparms:
            kin = NP.random.normal(loc=psrparms["KIN"], scale=psrparms["KIN_ERR"], size=n_samples)
            sini = NP.sin(kin*NP.pi/180)
            inc=kin
        else:
            sini = NP.random.normal(loc=psrparms["SINI"], scale=psrparms["SINI_ERR"], size=n_samples)
            inc = NP.arcsin(sini)*180/NP.pi

        cut = NP.argwhere((m2 > mass_func) * (m2 < 10))
        m2 = m2[cut]
        sini = sini[cut]
        inc = NP.arcsin(sini)*180/NP.pi
        mtot = (m2 * sini)**3 / mass_func
        mp = NP.sqrt(mtot[is_valid(mtot)*is_valid(m2)]) - m2[is_valid(mtot)*is_valid(m2)]
        mp = mp[is_valid(mp)]
        mtot = mtot[is_valid(mtot)]
        inc = inc[is_valid(inc)]
        
    return {'M2': NP.median(m2), 'M2_std': NP.std(m2), 'M2_lolim': NP.percentile(m2, q=16.0), 'M2_uplim': NP.percentile(m2, q=84.0), 'Mpsr': NP.median(mp), 'Mpsr_std': NP.std(mp), 'Mpsr_lolim': NP.percentile(mp, q=16.0), 'Mpsr_uplim': NP.percentile(mp, q=84.0), 'Mtot': NP.median(mtot), 'Mtot_std': NP.std(mtot), 'Mtot_lolim': NP.percentile(mtot, q=16.0), 'Mtot_uplim': NP.percentile(mtot, q=84.0), 'inc': NP.median(inc), 'inc_std': NP.std(inc), 'inc_lolim': NP.percentile(inc, q=16.0), 'inc_uplim': NP.percentile(inc, q=84.0)}
