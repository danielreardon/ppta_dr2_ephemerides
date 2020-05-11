#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:53:00 2019

@author: dreardon

Runs basic white, red, and DM noise model for all pulsars in datadir
"""

import numpy as np
from enterprise_extensions import models, model_utils
import matplotlib.pyplot as plt
import glob
from ppta_utils import *

datadir = "./Datasets/dr2_best7/"
parfiles = sorted(glob.glob(datadir + '/*.par'))
timfiles = sorted(glob.glob(datadir + '/*.tim'))

psrs = load_pulsars(datadir, ephem='DE436')

for psr in psrs:
    print("Running single psr analysis for: ", psr.name)
    # Exponential DM dip for J1713
    dm_expdip = True if psr.name == 'J1713+0747' else False
    refit = False if psr.name == 'J0437-4715' else True
    dm_expdip_tmin = 54500
    dm_expdip_tmax = 54900
    plt.errorbar(psr.toas/86400, psr.residuals*1e6, psr.toaerrs*1e6, fmt='.')
    plt.show()
    # Create a single pulsar model
    pta = models.model_singlepsr_noise(psr, psd='powerlaw', red_var=True,
                                       white_vary=True, dm_var=True,
                                       dm_psd='powerlaw', dm_annual=False,
                                       dm_expdip=dm_expdip,
                                       dm_expdip_tmin=dm_expdip_tmin,
                                       dm_expdip_tmax=dm_expdip_tmax)
    print(pta.params)

    outdir = datadir + "/chains/singlePsrNoise/" + psr.name
    pp = default_sampler(pta, N=3e5, outdir=outdir, plot=True, refit=refit)

    # Now, save noise files
    make_noise_files(psr.name, pp.chain, pp.pars,
                     outdir=datadir+'/noiseFiles/')
