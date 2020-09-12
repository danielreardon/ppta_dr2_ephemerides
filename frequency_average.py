#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 11:11:55 2020

@author: dreardon

Read white noise parameters from .par file and apply to TOAs .tim file

Average subbands together, and apply ECORR

Optionally, remove a linear trend as a function of frequency for each
    system, before averaging

The procedure is as follows:
    0) We assume the white noise model is correct, and the DM is fixed at
        a sensible value to allow accurate measurement of the FD parameters
    1) Read white noise parameters from .par file and apply EFACs, EQUADs
        to errors, and (optionally) FD parameters to TOAs.
    2) Optionally, for each system, fit and remove a straight line as a
        function of frequency, to correct for system-dependent profile
        evolution
    3) For each filename, collect all subbands and perform a weighted average
        of their centre frequencies, TOAs, and TOA errors
    4) Apply ECORRs to the averaged TOAs
    5) Save a new .tim file with averaged TOAs

"""
from decimal import Decimal, InvalidOperation
import numpy as np
from lmfit import Parameters, Minimizer
import matplotlib.pyplot as plt
import glob


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


def read_tim(timfile):
    """
    Reads a .tim file and saves non-commented lines to an array
    """
    f = open(timfile, 'r')
    lines = f.readlines()
    tim = []
    for line in lines:
        line = line.strip()
        # if not header or commented out
        if line[0] not in ['#', 'C', 'F', 'M']:
            tim.append(line)

    return tim


def fitter(model, params, args):
    # Do fit
    func = Minimizer(model, params, fcn_args=args)
    results = func.minimize()
    return results


def polynomial(params, x, y, w):
    model = params['p0'] + params['p1']*x + params['p2']*x**2 + \
         + params['p3']*x**3 + params['p4']*x**4 + params['p5']*x**5
    return (y - model) * w


def get_data(tim, flag="-group"):
    """
    Returns necessary data from tim array:
        filename, frequency, MJD, TOAerr, flag_value
    """
    files = []
    freqs = []
    toas = []
    errs = []
    flag_vals = []
    for line in tim:
        line = line.split()
        files.append(line[0])
        freqs.append(Decimal(line[1]))
        toas.append(Decimal(line[2]))
        errs.append(Decimal(line[3]))
        ind = np.argwhere([l.strip() == flag for l in line])
        flag_vals.append(line[int(ind + 1)])
    files = np.array(files).squeeze()
    freqs = np.array(freqs).squeeze()
    toas = np.array(toas).squeeze()
    errs = np.array(errs).squeeze()
    flag_vals = np.array(flag_vals).squeeze()

    return files, freqs, toas, errs, flag_vals


def dedisperse(toas, freqs, params, dm=None, ref_freq=1000, reverse=False):
    """
    Add or subtract DM delay from TOAs
    """
    if dm is None:
        dm = Decimal(params['DM'])
    dm_const = Decimal(2.41e-4)
    dm_delay = dm / dm_const / (freqs / Decimal(ref_freq))**2

    if reverse:
        return toas + dm_delay / Decimal(86400) / Decimal(1e6)
    else:
        return toas - dm_delay / Decimal(86400) / Decimal(1e6)


def apply_efac_equad(params, errs, flag_vals):
    """
    Applies efacs (TNEF) and equads (TNEQ) from .par file to toa errors

        errs_new = sqrt( (TNEF*errs)**2 + TNEQ )
    """
    # Apply EFACs first
    for key in params.keys():
        if "TNEF" in key and not ("ERR" in key or "TYPE" in key):
            efac = Decimal(params[key])  # efac value
            group = '_'.join(key.split('_')[2:])

            print('Applying EFAC of {0} to {1}'.format(round(efac, 3), group))

            inds = np.argwhere(flag_vals == group)
            errs[inds] *= efac
    print(" ")

    # Now apply EQUADs
    for key in params.keys():
        if "TNEQ" in key and not ("ERR" in key or "TYPE" in key):
            equad = Decimal(10**(params[key] + 6))  # equad value
            group = '_'.join(key.split('_')[2:])

            print('Applying EQUAD of {0}us to {1}'.format(round(equad, 3),
                  group))

            inds = np.argwhere(flag_vals == group)
            errs[inds] = np.sqrt(errs[inds]**2 + equad**2)
    print(" ")

    return errs


def apply_fd(params, toas, freqs, ref_freq=1000):
    """
    Read FD parameters from .par file and apply to TOAs
    Applies FD parameters with respect to reference frequencyref_freq (MHz)
    """
    polyorder = []
    fd_vals = []
    for key in params.keys():
        if key[0:2] == 'FD' and len(key) == 3:
            polyorder.append(int(key[2]))
            fd_vals.append(float(params[key]))  # FD value
    polyorder = np.array(polyorder)
    fd_vals = np.array(fd_vals)

    if polyorder.size == 0:
        return toas

    print("Applying FD parameters")
    print("FD polynomial orders:", polyorder)
    print("Corresponding FD values:", fd_vals)

    correction = np.zeros(np.shape(freqs))
    freqs = np.array([float(f) for f in freqs]).squeeze()
    for ii in range(0, polyorder.size):
        p = polyorder[ii]
        c = fd_vals[ii]
        correction += c * np.log(freqs/ref_freq)**p

    mx = round(max(correction*10**6), 3)
    mn = round(min(correction*10**6), 3)
    print("Correction ranges from {0}us to {1}us".format(mx, mn))
    print(" ")

    # convert correction into Decimal days to apply to TOAs in MJD
    correction = np.array([Decimal(c) for c in correction]).squeeze()
    correction /= 86400

    return toas + correction


def correct_systems_fd(tim, toas, files, freqs, errs, flag='-h', plot=True,
                       max_polyorder=5, downweight_outliers=True, nsig=5,
                       maxoutliers=0.05):
    """
    For each group of TOAs defined by a flag value, fit and remove a
        linear trend as a function of frequency

    Remove time-dependence by subtracting the weighted mean TOA from each
        group of subbands in a single epoch
    """
    # Get flag_values
    flag_vals = []
    for line in tim:
        line = line.split()
        ind = np.argwhere([l.strip() == flag for l in line])
        flag_vals.append(line[int(ind + 1)])
    flag_vals = np.array(flag_vals).squeeze()

    for val in np.unique(flag_vals):
        print("FD system correction for {0} {1}".format(flag, val))
        flaginds = np.argwhere(flag_vals == val).squeeze()
        ifiles = files[flaginds].squeeze()
        itoas = toas[flaginds].squeeze()
        ifreqs = freqs[flaginds].squeeze()
        ierrs = errs[flaginds].squeeze()

        # Now remove all time-structure by subtracting mean TOA
        for file in np.unique(ifiles):
            inds = np.argwhere(ifiles == file).squeeze()
            if np.size(np.array(inds)) > 1:
                avg = np.average(itoas[inds], weights=1/ierrs[inds]**2)
                itoas[inds] -= avg
            else:
                itoas[inds] = 0
        if plot:
            plt.errorbar(ifreqs, itoas * 86400 * 10**6,
                         yerr=ierrs, fmt='.', color='mediumblue', alpha=0.5)

        x = np.array([float(f) for f in ifreqs]).squeeze()
        y = np.array([float(t) for t in itoas * 86400 * 10**6]).squeeze()
        w = np.array([float(1/e) for e in ierrs]).squeeze()

        chisqr = np.inf
        for o in range(0, max_polyorder + 1):
            params = Parameters()
            params.add('p0', value=0, vary=True, min=-np.inf, max=np.inf)
            params.add('p1', value=0, vary=False, min=-np.inf, max=np.inf)
            if o >= 1:
                params['p1'].vary = True
            params.add('p2', value=0, vary=False, min=-np.inf, max=np.inf)
            if o >= 2:
                params['p2'].vary = True
            params.add('p3', value=0, vary=False, min=-np.inf, max=np.inf)
            if o >= 3:
                params['p3'].vary = True
            params.add('p4', value=0, vary=False, min=-np.inf, max=np.inf)
            if o >= 4:
                params['p4'].vary = True
            params.add('p5', value=0, vary=False, min=-np.inf, max=np.inf)
            if o >= 5:
                params['p5'].vary = True

            results_new = fitter(polynomial, params, (x, y, w))
            if results_new.chisqr < chisqr - 2:
                chisqr = results_new.chisqr
                results = results_new
                order = o
            else:
                break
        print("Removed order-{0} polynomial".format(order))

        # subtract model from toas
        model = -polynomial(results.params, x, np.zeros(np.shape(y)),
                            np.ones(np.shape(w))).squeeze()

        if downweight_outliers:
            print("Down-weighting {0}-sigma outlying subbands".format(nsig))
            outlier_inds = np.argwhere(np.abs(y - model) > nsig * ierrs).\
                squeeze()
            if np.size(outlier_inds) > maxoutliers*np.size(itoas):
                print(min(toas[flaginds][outlier_inds]))
                print(max(toas[flaginds][outlier_inds]))
                outlier_inds = []
                print("Too many outliers. Flag selection requires checking")
            print("{0} sub-bands ({1}%) downweighted".
                  format(np.size(outlier_inds),
                         round(100*np.size(outlier_inds)/np.size(itoas), 3)))
            errs = np.array([float(e) for e in errs]).squeeze()
            extra = np.zeros(np.shape(errs[flaginds]))
            extra[outlier_inds] = 1e6
            errs[flaginds] += extra
            errs = np.array([Decimal(e) for e in errs]).squeeze()

            if plot:
                plt.errorbar(ifreqs[outlier_inds], itoas[outlier_inds]
                             * 86400 * 10**6, yerr=ierrs[outlier_inds],
                             fmt='x', color='crimson')

        model = np.array([Decimal(m) / 86400 / Decimal(10**6) for m
                          in model]).squeeze()
        toas[flaginds] -= model

        if plot:
            xplot = np.linspace(min(x), max(x), 1000)
            modelplot = -polynomial(results.params, xplot,
                                    np.zeros(np.shape(xplot)),
                                    np.ones(np.shape(xplot)))
            plt.plot(xplot, modelplot, 'k')
            plt.title(r'{0} {1}'.format(flag, val.replace('_', '-')))
            plt.xlabel('Freq')
            plt.ylabel('TOA - mean (us)')
            plt.show()

    return toas, errs


def form_avg(files, freqs, toas, errs, flag_vals):
    """
    For each file, find the weighted average freq, toa, and error for all
        sub-bands
    """
    files_avg = []
    freqs_avg = []
    toas_avg = []
    errs_avg = []
    flag_vals_avg = []

    for file in np.unique(files):
        indicies = np.argwhere(files == file)
        files_avg.append(files[indicies][0])
        flag_vals_avg.append(flag_vals[indicies][0])

        toas_avg.append(Decimal(np.average(toas[indicies],
                                weights=1/(np.array(errs[indicies])**2))))
        errs_avg.append(Decimal(np.average(errs[indicies],
                        weights=1/(np.array(errs[indicies])**2))) /
                        Decimal(np.sqrt(len(errs[indicies]))))
        freqs_avg.append(Decimal(np.average(freqs[indicies],
                                 weights=1/(np.array(errs[indicies])**2))))

    files = np.array(files_avg).squeeze()
    freqs = np.array(freqs_avg).squeeze()
    toas = np.array(toas_avg).squeeze()
    errs = np.array(errs_avg).squeeze()
    flag_vals = np.array(flag_vals_avg).squeeze()

    return files, freqs, toas, errs, flag_vals


def apply_ecorr(params, errs, flag_vals):
    """
    Applies ecorrs (TNECORR) from .par file to toa errors
    """

    for key in params.keys():
        if "TNECORR" in key and not ("ERR" in key or "TYPE" in key):
            ecorr = Decimal(params[key])  # equad value
            group = '_'.join(key.split('_')[2:])

            print('Applying ECORR of {0}us to {1}'.format(round(ecorr, 3),
                  group))

            inds = np.argwhere(flag_vals == group)
            errs[inds] = np.sqrt(errs[inds]**2 + ecorr**2)
    print(" ")

    return errs


def write_timfile(outfile, tim, files, freqs, toas, errs, flag_vals):
    """
    Output new averaged data to a new .tim file
    """
    with open(outfile, "w+") as f:
        f.write("FORMAT 1\n")
        f.write("MODE 1\n")
        for ii in range(0, len(files)):
            ifile = files[ii]
            for line in tim:
                lsplit = line.split()
                if lsplit[0] == ifile:
                    writeline = ' '.join(lsplit[4:])
            f.write("{0} {1} {2} {3} {4}\n".
                    format(ifile, round(freqs[ii], 8), toas[ii],
                           round(errs[ii], 6), writeline))


def average_timfile(parfile, timfile, outfile=None, fd_correct=True,
                    fd_systems=False, white_flag='-group', fd_flag='-h',
                    plot=True, max_polyorder=5, downweight=True,
                    nsig=5):
    """
    Reads a .par file and .tim file, and frequency-averages the TOAs,
        saving output to outfile
    Optionally, perform a weighted fit of a linear term as a function of
        log(frequency)per system, before averaging. Remove linear trend if
        slope is inconsistent with 0
    """

    if outfile is None:
        outfile = timfile.replace('.tim', '.avg.tim')

    # Read .par file
    params = read_par(parfile)
    # Read tim file
    tim = read_tim(timfile)
    # Get necessary data from .tim file
    files, freqs, toas, errs, flag_vals = get_data(tim, flag=white_flag)
    # Apply EFACs and EQUADs to subbanded data
    errs = apply_efac_equad(params, errs, flag_vals)
    # De-disperse
    toas = dedisperse(toas, freqs, params)
    # Correct TOAs for FD parameters
    if fd_correct:
        toas = apply_fd(params, toas, freqs)
    if fd_systems:
        toas, errs = correct_systems_fd(tim, toas, files, freqs, errs,
                                        flag=fd_flag, plot=plot,
                                        max_polyorder=max_polyorder,
                                        downweight_outliers=downweight,
                                        nsig=nsig)

    files, freqs, toas, errs, flag_vals = form_avg(files, freqs, toas, errs,
                                                   flag_vals)
    errs = apply_ecorr(params, errs, flag_vals)

    # De-dedisperse
    toas = dedisperse(toas, freqs, params, reverse=True)

    write_timfile(outfile, tim, files, freqs, toas, errs, flag_vals)

    return


datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/'
# parfile = datadir + 'J0437/new.dr2.par'
# timfile = datadir + 'J0437/J0437-4715.tim'

parfiles = sorted(glob.glob(datadir + 'final/averaged/*.par'))
timfiles = sorted(glob.glob(datadir + 'final/averaged/*.tim'))

for ii in range(0, len(timfiles)):
    parfile = parfiles[ii]
    timfile = timfiles[ii]
    print(parfile)
    print(timfile)
    average_timfile(parfile=parfile, timfile=timfile, outfile=None,
                    fd_correct=True, fd_systems=False, white_flag='-group',
                    fd_flag='-group')
