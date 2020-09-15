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
import time


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


def fitter(model, params, args, mcmc=False):
    # Do fit
    if not mcmc:
        method = 'leastsq'
    else:
        method = 'emcee'
    func = Minimizer(model, params, fcn_args=args)
    results = func.minimize(method=method)
    return results


def polynomial(params, x, y, w):
    model = params['p0'] + params['p1']*x + params['p2']*x**2 + \
            params['p3']*x**3 + params['p4']*x**4 + params['p5']*x**5 + \
            params['p6']*x**6 + params['p7']*x**7 + params['p8']*x**8
    return (y - model) * w


def dm_model(params, freqs, toas, files, pars, w):
    """
    Determine the best DM
    """
    dm = params['dm']
    toas = dedisperse(toas, freqs, pars, dm=dm)

    # Now remove all time-structure by subtracting mean TOA
    for file in np.unique(files):
        inds = np.argwhere(files == file).squeeze()
        if np.size(np.array(inds)) > 1:
            avg = np.average(toas[inds], weights=w[inds]**2)
            toas[inds] -= avg
        else:
            toas[inds] = 0

    toas = shift_turns(pars, toas, files)

    # convert to float array for lmfit
    toas = np.array([float(t) for t in toas]).squeeze()
    weights = np.copy(w)
    w = np.array([float(w) for w in weights]).squeeze()
    return (toas) * w


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


def dedisperse(toas, freqs, params, dm=None, ref_freq=1000, reverse=False,
               plot=False):
    """
    Add or subtract DM delay from TOAs
    """
    if dm is None:
        dm = Decimal(params['DM'])
    else:
        dm = Decimal(float(dm))
    dm_const = Decimal(2.41e-4)
    dm_delay = dm / dm_const / (freqs / Decimal(ref_freq))**2

    period_us = 10**6/Decimal(params['F0'])
    while np.any(dm_delay >= period_us):
        dm_delay[dm_delay >= period_us] -= period_us

    if plot:
        plt.scatter(freqs, dm_delay, color='mediumblue')
        plt.xlabel('Frequency')
        plt.ylabel('Dispersion delay (us)')
        plt.show()

    delay_mjd = dm_delay / Decimal(86400) / Decimal(1e6)

    if reverse:
        print('De-dedispersion with DM={0}'.format(dm))
        return toas + delay_mjd
    else:
        # print('De-dispersion with DM={0}'.format(dm))
        return toas - delay_mjd


def apply_efac_equad(params, errs, flag_vals, outfile=None):
    """
    Applies efacs (TNEF) and equads (TNEQ) from .par file to toa errors

        errs_new = sqrt( (TNEF*errs)**2 + TNEQ )
    """
    # Apply EFACs first
    for key in params.keys():
        if "TNEF" in key and not ("ERR" in key or "TYPE" in key):
            efac = Decimal(params[key])  # efac value
            group = '_'.join(key.split('_')[2:])

            if outfile is not None:
                with open(outfile, "a+") as f:
                    f.write('# Applied TNEF of {0} to {1} \n'.
                            format(round(efac, 3), group))
            print('Applying EFAC of {0} to {1}'.format(round(efac, 3), group))

            inds = np.argwhere(flag_vals == group)
            errs[inds] *= efac
    print(" ")

    # Now apply EQUADs
    for key in params.keys():
        if "TNEQ" in key and not ("ERR" in key or "TYPE" in key):
            equad = Decimal(10**(params[key] + 6))  # equad value
            group = '_'.join(key.split('_')[2:])

            if outfile is not None:
                with open(outfile, "a+") as f:
                    f.write('# Applied TNEQ of {0}us to {1}\n'.
                            format(round(equad, 3), group))
            print('Applying EQUAD of {0}us to {1}'.format(round(equad, 3),
                  group))

            inds = np.argwhere(flag_vals == group)
            errs[inds] = np.sqrt(errs[inds]**2 + equad**2)
    print(" ")

    if outfile is not None:
        with open(outfile, "a+") as f:
            f.write('# \n')

    return errs


def apply_fd(params, toas, freqs, ref_freq=1000, plot=True, outfile=None):
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

    if outfile is not None:
        with open(outfile, "a+") as f:
            f.write('# Applied FD correction using: \n')
            for ii in range(0, polyorder.size):
                f.write("# FD{0} = {1} \n".format(polyorder[ii], fd_vals[ii]))
            f.write('# \n')

    print("Applying FD parameters")
    print("FD polynomial orders:", polyorder)
    print("Corresponding FD values:", fd_vals)

    correction = np.zeros(np.shape(freqs))
    freqs = np.array([float(f) for f in freqs]).squeeze()
    for ii in range(0, polyorder.size):
        p = polyorder[ii]
        c = fd_vals[ii]
        correction += c * np.log(freqs/ref_freq)**p * 10**6

    period_us = 10**6/Decimal(params['F0'])
    while np.any(correction > period_us):
        correction[correction > period_us] -= period_us

    if plot:
        plt.scatter(freqs, correction, color='mediumblue')
        plt.xlabel('Frequency')
        plt.ylabel('FD correction (us)')
        plt.show()

    mx = round(max(correction), 3)
    mn = round(min(correction), 3)
    print("Correction ranges from {0}us to {1}us".format(mx, mn))
    print(" ")

    # convert correction into Decimal days to apply to TOAs in MJD
    correction = np.array([Decimal(c) for c in correction]).squeeze()
    correction /= (86400 * 10**6)

    return toas - correction


def correct_systems_fd(tim, toas, files, freqs, errs, flag='-h', plot=True,
                       max_polyorder=8, downweight_outliers=True, nsig=5,
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
        try:
            flag_vals.append(line[int(ind + 1)])
        except TypeError:
            flag_vals.append('None')
    flag_vals = np.array(flag_vals).squeeze()

    for val in np.unique(flag_vals):
        if val == 'None':
            continue
        print("FD system correction for {0} {1}".format(flag, val))
        flaginds = np.argwhere(flag_vals == val).squeeze()
        if flaginds.size < 20:
            print('Too few TOAs')
            continue
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
                plt.errorbar(ifreqs[inds], itoas[inds] * 86400 * 10**6,
                             yerr=ierrs[inds], fmt='.', alpha=0.5)

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
            params.add('p6', value=0, vary=False, min=-np.inf, max=np.inf)
            if o >= 6:
                params['p6'].vary = True
            params.add('p7', value=0, vary=False, min=-np.inf, max=np.inf)
            if o >= 7:
                params['p7'].vary = True
            params.add('p8', value=0, vary=False, min=-np.inf, max=np.inf)
            if o >= 8:
                params['p8'].vary = True

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


def apply_ecorr(params, errs, flag_vals, outfile=None):
    """
    Applies ecorrs (TNECORR) from .par file to toa errors
    """

    for key in params.keys():
        if "TNECORR" in key and not ("ERR" in key or "TYPE" in key):
            ecorr = Decimal(params[key])  # equad value
            group = '_'.join(key.split('_')[2:])

            if outfile is not None:
                with open(outfile, "a+") as f:
                    f.write('# Applied TNECORR of {0}us to {1} \n'.
                            format(round(ecorr, 3), group))

            print('Applying ECORR of {0}us to {1}'.format(round(ecorr, 3),
                  group))

            inds = np.argwhere(flag_vals == group)
            errs[inds] = np.sqrt(errs[inds]**2 + ecorr**2)
    print(" ")

    if outfile is not None:
        with open(outfile, "a+") as f:
            f.write("# \n")

    return errs


def write_timfile(outfile, tim, files, freqs, toas, errs, flag_vals):
    """
    Output new averaged data to a new .tim file
    """
    with open(outfile, "a+") as f:
        f.write("FORMAT 1\n")
        f.write("MODE 1\n")
        sortind = np.argsort(toas)
        for ii in range(0, len(files)):
            ifile = files[sortind][ii]
            for line in tim:
                lsplit = line.split()
                if lsplit[0] == ifile:
                    writeline = ' '.join(lsplit[4:])
            f.write("{0} {1} {2} {3} {4}\n".
                    format(ifile, round(freqs[sortind][ii], 8),
                           toas[sortind][ii],
                           round(errs[sortind][ii], 6), writeline))
    return


def shift_turns(params, toas, files, plot=False):
    """
    For each group of sub-bands, shift TOAs to the same pulse number
    """

    period_us = 10**6/Decimal(params['F0'])
    period_mjd = Decimal(period_us / 10**6 / 86400)

    for file in np.unique(files):
        inds = np.argwhere(files == file).squeeze()
        if np.size(np.array(inds)) > 1:
            dt = np.diff(toas[inds] - np.median(toas[inds]))

            while np.any(dt > period_mjd/2):
                turns = np.zeros(np.shape(toas[inds]))
                turns = np.array([Decimal(t) for t in turns]).squeeze()
                turnind = np.argwhere(dt > period_mjd/2)
                turns[turnind] += period_mjd
                toas[inds] += turns
                dt = np.diff(toas[inds])
            while np.any(dt < -period_mjd/2):
                turns = np.zeros(np.shape(toas[inds]))
                turns = np.array([Decimal(t) for t in turns]).squeeze()
                turnind = np.argwhere(dt < -period_mjd/2)
                turns[turnind] -= period_mjd
                toas[inds] += turns
                dt = np.diff(toas[inds])
    return toas


def average_timfile(parfile, timfile, outfile=None, fd_correct=True,
                    fd_systems=False, white_flag='-group', fd_flag='-h',
                    plot=False, max_polyorder=8, downweight=True,
                    nsig=5, dm=None, fit_dm=False, ndm=100):
    """
    Reads a .par file and .tim file, and frequency-averages the TOAs,
        saving output to outfile
    Optionally, perform a weighted fit of a linear term as a function of
        log(frequency)per system, before averaging. Remove linear trend if
        slope is inconsistent with 0
    """

    if outfile is None:
        outfile = timfile.replace('.tim', '.avg.tim')

    with open(outfile, "w+") as f:
        now = time.gmtime()
        f.write("# Averaged .tim file prepared at {0}:{1}:{2} (UTC, 12h) on {3}/{4}/{5} \n".
                format(now.tm_hour,
                       now.tm_min,
                       now.tm_sec,
                       now.tm_mday,
                       now.tm_mon,
                       now.tm_year))
        f.write("# Using Daniel Reardon's 'frequency_average.py' for PPTA-DR2 \n")
        f.write("# \n")
    with open(outfile, "a+") as f:
        f.write("# White noise flag is {0} \n".format(white_flag))
    # Read .par file
    params = read_par(parfile)
    # Read tim file
    tim = read_tim(timfile)
    # Get necessary data from .tim file
    files, freqs, toas, errs, flag_vals = get_data(tim, flag=white_flag)
    # Apply EFACs and EQUADs to subbanded data
    errs = apply_efac_equad(params, errs, flag_vals, outfile=outfile)
    # Correct TOAs for FD parameters
    if fd_correct:
        toas = apply_fd(params, toas, freqs, plot=plot, outfile=outfile)

    if dm is None:
        dm = params['DM']

    if fit_dm:
        print("Original DM = {0}. Now fitting".format(dm))
        # reverse-engineer folding DM value
        dm_params = Parameters()
        dm_params.add('dm', value=dm, vary=True,
                      min=-np.inf, max=np.inf)
        results = fitter(dm_model, dm_params,
                         (freqs, toas, files, params, 1/errs), mcmc=False)
        dm = results.params['dm'].value
        print("Fitted DM".format(dm))

    with open(outfile, "a+") as f:
        f.write("# Dedispersed with DM = {0} \n".format(round(dm, 6)))
        f.write("# \n")

    # De-disperse
    toas = dedisperse(toas, freqs, params, plot=plot, dm=dm)
    # Shift TOAs by pulse period if needed
    toas = shift_turns(params, toas, files, plot=plot)

    # Now remove all time-structure by subtracting mean TOA
    if plot:
        for file in np.unique(files):
            toas_plot = np.copy(toas)
            inds = np.argwhere(files == file).squeeze()
            if np.size(np.array(inds)) > 1:
                avg = np.average(toas_plot[inds], weights=1/errs[inds]**2)
                toas_plot[inds] -= avg
            else:
                toas_plot[inds] = 0
            if plot:
                plt.errorbar(freqs[inds], toas_plot[inds] * 86400 * 10**6,
                             yerr=errs[inds], fmt='.', alpha=0.5)

        plt.ylabel('Residual')
        plt.xlabel('Frequency')
        plt.show()

    if fd_systems:
        toas, errs = correct_systems_fd(tim, toas, files, freqs, errs,
                                        flag=fd_flag, plot=plot,
                                        max_polyorder=max_polyorder,
                                        downweight_outliers=downweight,
                                        nsig=nsig)

    files, freqs, toas, errs, flag_vals = form_avg(files, freqs, toas, errs,
                                                   flag_vals)
    errs = apply_ecorr(params, errs, flag_vals, outfile)

    # De-dedisperse
    toas = dedisperse(toas, freqs, params, reverse=True, dm=dm)

    write_timfile(outfile, tim, files, freqs, toas, errs, flag_vals)

    return


datadir = '/Users/dreardon/Dropbox/Git/ppta_dr2_ephemerides/'
parfiles = [datadir + 'J0437/J0437-4715.dr2e.par']
timfiles = [datadir + 'J0437/J0437-4715.dr2e.tim']

#parfiles = sorted(glob.glob(datadir +
#                  'final/tempo2/J1909*.par'))
#timfiles = sorted(glob.glob(datadir +
#                  'final/tempo2/J1909*.tim'))

for ii in range(0, len(timfiles)):
    parfile = parfiles[ii]
    timfile = timfiles[ii]
    print(parfile)
    print(timfile)
    average_timfile(parfile=parfile, timfile=timfile, outfile=None,
                    fd_correct=True, fd_systems=False, white_flag='-group',
                    fd_flag='-h', dm=None, fit_dm=True, plot=False,
                    downweight=False)
