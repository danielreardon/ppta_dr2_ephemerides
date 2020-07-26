#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:31:56 2019

@author: dreardon

make_plots.py: Makes residual plots for general2 outputs
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
from astropy.time import Time
import sys


def read_general2(filename):
    """
    Reads a general2 output into a numpy array
    """
    with open(filename, "r") as file:
        data = []
        files = []
        header = True
        for line in file:
            if not header and not ('Finish' in line):
                files.append([(line.split('\t')[0])])
                data.append([float(x) for x in
                             (line.replace('\n', '').split('\t')[1:])])
            elif 'Starting general2 plugin' in line:
                header = False
    return np.array(data), files


def wrms(data, weights):
    """
    Given some data and weights, return the weighted rms
    """
    return np.sqrt(np.cov(np.array(data).squeeze(),
                          aweights=np.array(weights).squeeze()))


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


def makeplot(psrname, plotnum, num_tot, count, xdata, ydata, errs, freqs, files,
             average=False, cadence=False, alpha=0.1, bestband=False, scale_dm=False):
    """
    Generates pretty plots for dr2!

        If colour inputs are given, assume we're making the cadence plot
    """

    band_change_date = 2009.5

    if average:
        xdata, ydata, errs, freqs = average_subbands(xdata, ydata, errs, freqs,
                                                     files)

    if scale_dm:
        ydata = ydata*freqs**2
        errs = errs*freqs**2
        ydata = ydata-np.mean(ydata)

    plt.figure(plotnum, figsize=(8.27, 11.69))
    plt.subplot(num_tot, 1, count)

    if bestband:
        indicies = np.argwhere((freqs > 1000)*(freqs < 2000))
        try:
            wrms_20 = wrms(ydata[indicies], 1/errs[indicies]**2)
            if len(ydata[indicies])< 0.05*len(ydata) or np.ptp(xdata[indicies])<0.5*np.ptp(xdata):
                wrms_20 = np.inf
        except ZeroDivisionError:
            wrms_20 = np.inf
        indicies = np.argwhere(freqs > 2000)
        try:
            wrms_10 = wrms(ydata[indicies], 1/errs[indicies]**2)
            if len(ydata[indicies])< 0.05*len(ydata) or np.ptp(xdata[indicies])<0.5*np.ptp(xdata):
                wrms_10 = np.inf
        except ZeroDivisionError:
            wrms_10 = np.inf
        indicies = np.argwhere(freqs < 1000)
        try:
            wrms_4050 = wrms(ydata[indicies], 1/errs[indicies]**2)
            if len(ydata[indicies])< 0.05*len(ydata) or np.ptp(xdata[indicies])<0.5*np.ptp(xdata):
                wrms_4050 = np.inf
        except ZeroDivisionError:
            wrms_4050 = np.inf
        best = np.argmin([wrms_10, wrms_20, wrms_4050])

    # 10CM data
    indicies = np.argwhere(freqs > 2000)
    if not cadence:
        if not bestband or best == 0:
            xplt = xdata[indicies]
            plt.errorbar(xdata[indicies], ydata[indicies],
                         yerr=errs[indicies], fmt='none', alpha=alpha, zorder=2,
                         capsize=1, elinewidth=0.7, ecolor='mediumblue')
    else:
        plt.scatter(xdata[indicies], ydata[indicies],
                    marker='.', alpha=alpha, zorder=2,
                    color='mediumblue')
                    #color=colour[indicies, :].squeeze())

    # 40/50CM data
    indicies = np.argwhere(freqs < 1000)
    if not cadence:
        if not bestband or best == 2:
            xplt = xdata[indicies]
            plt.errorbar(xdata[indicies], ydata[indicies],
                         yerr=errs[indicies], fmt='none', alpha=alpha, zorder=0,
                         capsize=1, elinewidth=0.7, ecolor='crimson')
    else:
        plt.scatter(xdata[indicies], ydata[indicies],
                    marker='.', alpha=alpha, zorder=2,
                    color='crimson')
                    #color=colour[indicies, :].squeeze())

    # 20CM data
    indicies = np.argwhere((freqs > 1000)*(freqs < 2000))
    if not cadence:
        if not bestband or best == 1:
            xplt = xdata[indicies]
            plt.errorbar(xdata[indicies], ydata[indicies],
                         yerr=errs[indicies], fmt='none', alpha=alpha*0.8, zorder=1,
                         capsize=1, elinewidth=0.7, ecolor='darkcyan')
    else:
        plt.scatter(xdata[indicies], ydata[indicies],
                    marker='.', alpha=alpha*0.8, zorder=2,
                    color='darkcyan')
                    #color=colour[indicies, :].squeeze())


    if not bestband:
        tot_wrms = wrms(ydata, 1/errs**2)
    else:
        tot_wrms = np.min([wrms_10, wrms_20, wrms_4050])
    if not cadence:
        plt.tight_layout(pad=0)
#        if not bestband:
#            plt.ylim(min(ydata), max(ydata))
#        else:
#            plt.ylim(min(ydata[indicies]), max(ydata[indicies]))
        bottom, top = plt.ylim()  # return the current ylim
        left, right = plt.xlim()  # return the current xlim
        plt_range = top-bottom
        plt.xlim(2002.1, 2019)
        plt.text(2001.8, -0.08*plt_range,
                 r"{0:2.2f} $\mu$s".format(round(tot_wrms, 2)), fontsize=9, color='dimgray',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
        #plt.text(min(xdata)-1.4, 0,
        #         r"{0} $\mu$s".format(round(tot_wrms, 2)), fontsize=9)
        plt.plot([2003, 2018.5], [0, 0], 'k:', alpha=0.3, zorder=-1)
        #plt.plot([band_change_date, band_change_date], [bottom+0.15*plt_range, top-0.15*plt_range], 'k', alpha=0.3, zorder=-2)
    else:
        if psrname not in ['J1017-7156', 'J1125-6014', 'J1446-4701', 'J1545-4550', 'J1832-0836', 'J2241-5236']:
            legacy = np.loadtxt('/Users/dreardon/Desktop/Misc/finish_dr2/general2_output/dr1e/{0}.tim.output'.format(psrname))
            t = Time(legacy, format='mjd')
            legacy = t.byear  # observation year
            plt.scatter(legacy[legacy<2004], np.zeros(np.shape(legacy[legacy<2004])), marker='.',
                        alpha=alpha, zorder=-1, color='black')
        plt.tight_layout(pad=0)
        bottom, top = plt.ylim()  # return the current ylim
        left, right = plt.xlim()  # return the current xlim
        plt.xlim(1994, 2018.95)
        plt.xticks([1994, 1997, 2000, 2003, 2006, 2009, 2012, 2015, 2018])
        plt.ylim(-1, 1)
    plt.yticks([])
    if count < num_tot:
        plt.axis('off')
    else:
        plt.box(on=None)
        plt.xlabel('Year', fontsize=10.25)
        plt.xticks(fontsize=10)
    plt.text(2018.8, 0.2*plt.ylim()[0], psrname, fontsize=10)
    if not cadence:
        plt.subplots_adjust(wspace=None, hspace=0)


"""
Start of script
"""


data_dir = '/Users/dreardon/Desktop/Misc/finish_dr2/'
filenames = sorted(glob.glob(data_dir + 'general2_output/J*.output'))

count = 1
for file in filenames:

    # Read in data and change it somewhat
    psrname = file.split('/')[-1].split('.')[0]
    data, files = read_general2(file)

    # Create a column for year
    n = len(filenames)
    t = Time(data[:, 1], format='mjd')
    data[:, 1] = t.byear  # observation year
    # Fix posttn
    data[:, 6] = data[:, 6] - np.mean(data[:, 6])  # temponest postfit

    band_change_date = Time(55319, format='mjd').byear


    """
    Define data for use
    """

    sat = data[:, 0].squeeze()
    year = data[:, 1].squeeze()
    pre = data[:, 2].squeeze()*10**6
    post = data[:, 3].squeeze()*10**6
    errs = data[:, 4].squeeze()
    freqs = data[:, 5].squeeze()
    if count == 1:
        minfreq = min(freqs)
    posttn = data[:, 6].squeeze()*10**6
    tnred = data[:, 7].squeeze()
    tndm = data[:, 8].squeeze()
    cadence_data = np.array([0 if (f < 2000)*(f > 1000) else
                             -0.35*np.sign(f-1400) for f in freqs])

    """
    Now make some plots!
    """
    # Cadence
    makeplot(psrname, 1, n, count, year, cadence_data, errs, freqs, files,
             average=True, cadence=True)
    # Postfit residuals
    makeplot(psrname, 2, n, count, year, post, errs, freqs, files,
             average=True, alpha=0.5)
    # Noise-subtracted
    makeplot(psrname, 3, n, count, year, posttn, errs, freqs, files,
             average=True, alpha=0.5)
    # DM subtracted
    makeplot(psrname, 4, n, count, year, post-tndm*10**6, errs, freqs, files,
             average=True, alpha=0.5)
    # Red subtracted
    makeplot(psrname, 5, n, count, year, post-tnred*10**6, errs, freqs, files,
             average=True, alpha=0.5)
    # Postfit residuals in best band
    makeplot(psrname, 6, n, count, year, post, errs, freqs, files,
             average=True, alpha=0.8, bestband=True)
    # Postfit residuals in best band with noise subtracted
    makeplot(psrname, 7, n, count, year, posttn, errs, freqs, files,
             average=True, alpha=0.8, bestband=True)

    count += 1


"""
Save plots
"""

plt.figure(1)
plt.savefig(data_dir + 'plots/cadence.pdf', papertype='a4', bbox_inches='tight')
plt.savefig(data_dir + 'plots/cadence.png', papertype='a4', bbox_inches='tight')

plt.figure(2)
plt.savefig(data_dir + 'plots/postfit_avg.pdf', papertype='a4', bbox_inches='tight')
plt.savefig(data_dir + 'plots/postfit_avg.png', papertype='a4', bbox_inches='tight')


plt.figure(3)
plt.savefig(data_dir + 'plots/postfit_tn_avg.pdf', papertype='a4', bbox_inches='tight')
plt.savefig(data_dir + 'plots/postfit_tn_avg.png', papertype='a4', bbox_inches='tight')


plt.figure(4)
plt.savefig(data_dir + 'plots/postfit_tnsubdm_avg.pdf', papertype='a4',
            bbox_inches='tight')
plt.savefig(data_dir + 'plots/postfit_tnsubdm_avg.png', papertype='a4',
            bbox_inches='tight')

plt.figure(5)
plt.savefig(data_dir + 'plots/postfit_tnsubred_avg.pdf', papertype='a4',
            bbox_inches='tight')
plt.savefig(data_dir + 'plots/postfit_tnsubred_avg.png', papertype='a4',
            bbox_inches='tight')

plt.figure(6)
plt.savefig(data_dir + 'plots/postfit_avg_bestband.pdf', papertype='a4',
            bbox_inches='tight')
plt.savefig(data_dir + 'plots/postfit_avg_bestband.png', papertype='a4',
            bbox_inches='tight')

plt.figure(7)
plt.savefig(data_dir + 'plots/postfit_tn_avg_bestband.pdf', papertype='a4',
            bbox_inches='tight')
plt.savefig(data_dir + 'plots/postfit_tn_avg_bestband.png', papertype='a4',
            bbox_inches='tight')
