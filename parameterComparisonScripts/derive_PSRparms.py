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

n_samples = 1000000
# Define other useful constants
M_sun = 1.98847542e+30  # kg
Tsun = 4.926790353700459e-06  #s
rad_to_mas = 180*3600*1000/NP.pi
parsec_to_m = 3.08567758e+16
sec_per_year = 86400*365.2425

font = {'family' : 'serif',
        'weight' : 'medium',
        'size'   : 12}

matplotlib.rc('font', **font)

from matplotlib import rc
import os
os.environ["PATH"] += os.pathsep + '/usr/bin/'
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

def is_valid(array):
    """
    Returns boolean array of values that are finite an not nan
    """
    return NP.isfinite(array)*(~NP.isnan(array))

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
    elif ('BETA' in psrparms) and ('LAMBDA' in psrparms):
        c = SkyCoord(psrparms['LAMBDA'], psrparms['BETA'], frame=BarycentricTrueEcliptic, unit=(U.deg, U.deg))
    elif ('RAJ' in psrparms) and ('DECJ' in psrparms):
        c = SkyCoord(psrparms['RAJ'] + ' ' + psrparms['DECJ'], unit=(U.hourangle, U.deg))
    else:
        return -1

    #Convert to Galactic

    ldeg = c.galactic.l.value
    sigl = 0.0
    bdeg = c.galactic.b.value
    sigb = 0.0

    if (('ELAT' not in psrparms) or ('ELONG' not in psrparms)) and (('BETA' not in psrparms) or ('LAMBDA' not in psrparms)):

        #Convert to Ecliptic
        #ecliptic transformation:

        elat = c.barycentrictrueecliptic.lat.value
        elon = c.barycentrictrueecliptic.lon.value

        #proper motion transformation:

        if ('PMRA' not in psrparms) or ('PMDEC' not in psrparms):
            return -1

        pm = ICRS(ra=c.ra.deg*U.degree, dec=c.dec.deg*U.deg, pm_ra_cosdec = psrparms['PMRA']*U.mas/U.yr, pm_dec = psrparms['PMDEC']*U.mas/U.yr)
        pm_ecliptic = pm.transform_to(BarycentricTrueEcliptic)
        pmlat = pm_ecliptic.pm_lat.value
        pmlon = pm_ecliptic.pm_lon_coslat.value

    if ('PBDOT' in psrparms) and ('PBDOT_ERR' in psrparms):
        pbdot_posterior = NP.random.normal(loc=psrparms["PBDOT"], scale=psrparms["PBDOT_ERR"], size=n_samples)  # observed
        if 'PSRJ' in psrparms:
            psrname = psrparms['PSRJ']
        elif 'PSR' in psrparms:
            psrname = psrparms['PSR']

        if psrname == 'J1909-3744':
            pbdot_grav = -2.763*10**(-15)
        else:
            pbdot_grav = 0

        if ('PX' in psrparms) and 'PX_ERR' in psrparms:
            d_prior = 1/NP.random.normal(loc=psrparms["PX"], scale=psrparms["PX_ERR"], size=n_samples) # observed
            dkpc = NP.median(d_prior[(d_prior > 0)*(d_prior < 100)])
            sigd = dkpc*psrparms["PX_ERR"]/psrparms["PX"]
        else:
            dkpc = 1.0
            sigd = 0.2*dkpc
        pb = psrparms['PB']
        pb_err = psrparms['PB_ERR']
        if ('PMELAT' in psrparms) and ('PMELONG' in psrparms):
            pmelat = psrparms['PMELAT']
            pmelong = psrparms['PMELONG']
            pm_tot = NP.sqrt(pmelat**2 + pmelong**2)
            pm = pm_tot/(sec_per_year*rad_to_mas)
        elif ('PMBETA' in psrparms) and ('PMLAMBDA' in psrparms):
            pmelat = psrparms['PMBETA']
            pmelong = psrparms['PMLAMBDA']
            pm_tot = NP.sqrt(pmelat**2 + pmelong**2)
            pm = pm_tot/(sec_per_year*rad_to_mas)
        elif ('PMRA' in psrparms) and ('PMDEC' in psrparms):
            pmra = psrparms['PMRA']
            pmdec = psrparms['PMDEC']
            pm_tot = NP.sqrt(pmra**2 + pmdec**2)
            pm = pm_tot/(sec_per_year*rad_to_mas)

        # Expected Shklovskii and Galactic terms

        Ex_pl =  GalDynPsr.modelLb.Expl(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term parallel to the Galactic plane
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
    mtot2 = 0
    if 'BINARY' in psrparms:
        if ('A1' not in psrparms) or ('PB' not in psrparms):
            return -1
        mass_func = (4*NP.pi**2/FCNST.G) * (psrparms['A1']*FCNST.c)**3/(psrparms['PB']*86400)**2
        mass_func = mass_func/M_sun

    if (('H3' not in psrparms) or ('H3_ERR' not in psrparms)) and (('M2' not in psrparms) or ('M2_ERR' not in psrparms)):
        return -1

    if ('H3' in psrparms) and ('H3_ERR' in psrparms):
        h3 = psrparms["H3"]*10**6
        h3_err = psrparms["H3_ERR"]*10**6
        h3_arr = NP.random.normal(loc=h3, scale=h3_err, size=n_samples)
        if ('H4' in psrparms) and ('H4_ERR' in psrparms):
            h4 = psrparms["H4"]*10**6
            h4_err = psrparms["H4_ERR"]*10**6
            h4_arr = NP.random.normal(loc=h4, scale=h4_err, size=n_samples)
        elif ('STIG' in psrparms) and ('STIG_ERR' in psrparms):
            stig = psrparms["STIG"]
            stig_err = psrparms["STIG_ERR"]
            stig_arr = NP.random.normal(loc=stig, scale=stig_err, size=n_samples)
            h4_arr = stig_arr * h3_arr
        else:
            return -1
        sini = 2 * h3_arr * h4_arr / ( h3_arr**2 + h4_arr**2 )
        m2 = h3_arr**4 / h4_arr**3  # in us
        m2 = m2/(Tsun*10**6)
    elif ('M2' in psrparms) and ('M2_ERR' in psrparms):
        m2 = NP.random.normal(loc=psrparms["M2"], scale=psrparms["M2_ERR"], size=n_samples)
        if ('KIN' in psrparms) and ('KIN_ERR' in psrparms):
            kin = NP.random.normal(loc=psrparms["KIN"], scale=psrparms["KIN_ERR"], size=n_samples)
            sini = NP.sin(kin*NP.pi/180)
            inc=kin
        else:
            sini = NP.random.normal(loc=psrparms["SINI"], scale=psrparms["SINI_ERR"], size=n_samples)
            inc = NP.arcsin(sini)*180/NP.pi

    if ('XDOT' in psrparms):
        xdot = NP.random.normal(loc=params["XDOT"], scale=params["XDOT_ERR"], size=n_samples)
        a1 = NP.random.normal(loc=params["A1"], scale=params["A1_ERR"], size=n_samples)
        # get proper motion
        if 'ELAT' in params.keys():
            pm1 = NP.random.normal(loc=params["PMELAT"], scale=params["PMELAT_ERR"], size=n_samples)
            pm2 = NP.random.normal(loc=params["PMELONG"], scale=params["PMELONG_ERR"], size=n_samples)
        else:
            pm1 = NP.random.normal(loc=params["PMRA"], scale=params["PMRA_ERR"], size=n_samples)
            pm2 = NP.random.normal(loc=params["PMDEC"], scale=params["PMDEC_ERR"], size=n_samples)
        pm_tot = NP.sqrt(pm1**2 + pm2**2)
        pm = pm_tot/(sec_per_year*rad_to_mas)
        i_limit =  np.abs(np.arctan(a1 * pm / xdot))
        sini_lim = np.sin(np.median(i_limit))

    cut = NP.argwhere((m2 > mass_func) * (m2 < 1.4) * (sini < sini_lim))
    m2 = m2[cut]
    sini = sini[cut]
    inc = NP.arcsin(sini)*180/NP.pi
    mtot2 = (m2 * sini)**3 / mass_func
    mp = NP.sqrt(mtot2[is_valid(mtot2)*is_valid(m2)]) - m2[is_valid(mtot2)*is_valid(m2)]
    mp = mp[is_valid(mp)]
    mtot2 = mtot2[is_valid(mtot2)]
    mtot = NP.sqrt(mtot2[(mtot2 > 0)])
    inc = inc[is_valid(inc)]

    return {'M2': NP.median(m2), 'M2_std': NP.std(m2), 'M2_lolim': NP.percentile(m2, q=16.0), 'M2_uplim': NP.percentile(m2, q=84.0), 'Mpsr': NP.median(mp), 'Mpsr_std': NP.std(mp), 'Mpsr_lolim': NP.percentile(mp, q=16.0), 'Mpsr_uplim': NP.percentile(mp, q=84.0), 'Mtot': NP.median(mtot), 'Mtot_std': NP.std(mtot), 'Mtot_lolim': NP.percentile(mtot, q=16.0), 'Mtot_uplim': NP.percentile(mtot, q=84.0), 'inc': NP.median(inc), 'inc_std': NP.std(inc), 'inc_lolim': NP.percentile(inc, q=16.0), 'inc_uplim': NP.percentile(inc, q=84.0)}
