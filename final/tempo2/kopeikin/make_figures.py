#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer
import glob


# Path and shortname of temponest results, plus parameters of interest for plotting
#data_shortname = 'results/J1022-DM_annual-J1022+1001-'
#data_shortname = 'results/J1643-DM_annual-J1643-1224-'
data_shortname = 'results/J1909-DM_annual-J1909-3744-'
plot_params = ["KOM", "KIN"]  # plot_params = None  (If plotting everything)
#plot_params = ["OMDOT", "PX", "M2", "KIN", "PBDOT"]
derive_params = False  # option to derive posteriors for Mp and D for J0437-4715
par_file = data_shortname

# Define the temponest output files
chain_file = data_shortname + 'post_equal_weights.dat'
figure_file = data_shortname + 'triplot.pdf'
rescale_file = data_shortname + 'T2scaling.txt'
param_file = data_shortname + '.paramnames'

# Get parameter names
names = pd.read_csv(param_file, delim_whitespace=True, names=("sno", "param"), dtype={'sno': np.str, 'param': np.str})
cnames = names['param'].tolist()
cnames.append("Likelihood")
print(cnames)
new_names = []
new_plot_params = []
for name in cnames:  # rename parameters for plotting
    if name in plot_params:
        if name == "M2":
            new_plot_params.append("$M_c$")
            new_names.append("$M_c$")
        elif name == "KIN":
            new_plot_params.append("$i$")
            new_names.append("$i$")
        elif name == "KOM":
            new_plot_params.append("$\Omega$")
            new_names.append("$\Omega$")
        else:
            new_plot_params.append(name)
    else:
        new_names.append(name)

# Read data and rescale tempo2 parameters
posteriors = pd.read_csv(chain_file, names=cnames, delim_whitespace=True)
rescale_vals = pd.read_csv(rescale_file, delim_whitespace=True, names=("param", "mean", "std"), dtype={'param':np.str, 'mean':np.str, 'std':np.str})
for name in rescale_vals['param'].values:
    try:
        posteriors[name] = float(rescale_vals.loc[rescale_vals['param'] == name, 'mean']) + float(rescale_vals.loc[rescale_vals['param'] == name, 'std'])*posteriors[name]
    except:
        continue

kom = posteriors["KOM"].values*360
kin = posteriors["KIN"].values*180

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
x = kin
y = kom
# Define the borders
xmin = 0
xmax = 180
ymin = 0
ymax = 360
# Create meshgrid
for ii in [0, 1]:
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    if ii == 0:
        x = kin[np.argwhere(kin < 90).squeeze()]
        y = kom[np.argwhere(kin < 90).squeeze()]
    else:
        x = kin[np.argwhere(kin > 90).squeeze()]
        y = kom[np.argwhere(kin > 90).squeeze()]
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values, bw_method=0.1)
    f = np.reshape(kernel(positions).T, xx.shape)
    fig = plt.figure(1, figsize=(8,8))
    ax = fig.gca()
    #f, xedges, yedges = np.histogram2d(x, y, bins=20, range=[[xmin, xmax],[ymin, ymax]], density=True)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.scatter(x, y, marker='.', color='crimson')
    cset = ax.contour(xx, yy, f, colors = 'mediumblue', levels=[np.percentile(f, q=68.27), np.percentile(f, q=95.45), np.percentile(f, q=99.73)])
    #cset = ax.contour(xx, yy, f, colors = 'mediumblue', levels=[np.percentile(f, q=0.27), np.percentile(f, q=4.55), np.percentile(f, q=31.73)])
    #ax.scatter()
    #ax.imshow(np.rot90(f), extent=[xmin, xmax, ymin, ymax])
    #cset = ax.contour(xx, yy, f, colors='k')
    #ax.clabel(cset, fontsize=10)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.title('2D Gaussian Kernel density estimation')


#plt.figure(1, figsize=(9, 9))
#plt.scatter(kin, kom, marker='.', color='mediumblue')

plt.xlabel('i')
plt.ylabel('om')
plt.figure(2, figsize=(9, 9))
plt.hist(kom, bins=100, density=True)
plt.xlabel('kom')
plt.figure(3, figsize=(9, 9))
plt.hist(kin, bins=100, density=True)
plt.xlabel('kin')


plt.figure(1)
plt.savefig("kinkom.pdf")
plt.close()

plt.figure(2)
plt.savefig("kom.pdf")
plt.close()

plt.figure(3)
plt.savefig("kin.pdf")
plt.close()
exit()


# Add chain and plot
c = ChainConsumer()
c.add_chain(posteriors.values, parameters=new_names, bins=100, name="J1909-3744")
c.configure(diagonal_tick_labels=False, flip=False, sigma2d=True, kde=False, sigmas=[1,2,3], shade=True, shade_alpha=0.1, label_font_size=16, tick_font_size=16, max_ticks=3, summary=True)
c.plotter.plot(filename=figure_file, figsize=(16,16), legend=False, parameters=new_plot_params)
#plt.show()
plt.savefig("triplot.pdf")
#plt.close()

#exit()

# Setup for deriving parameters of interest for J0437 
if not derive_params:
    exit()
else:
    import libstempo as T
    import GalDynPsr
    import scipy.constants as sc
    import matplotlib.pyplot as plt
    import scipy.stats as stats

# Load pulsar
psr = T.tempopulsar(parfile=par_file, timfile=tim_file, maxobs=30000)

# Set parameters from par file
pbdot = psr["PBDOT"].val
pbdot_err = psr["PBDOT"].err
mc = psr["M2"].val
mc_err = psr["M2"].err
i = psr["KIN"].val
i_err = psr["KIN"].err
asini = psr["A1"].val
asini_err = psr["A1"].err
pb = psr["PB"].val
pb_err = psr["PB"].err
pmra = psr["PMRA"].val
pmra_err = psr["PMRA"].err
pmdec = psr["PMDEC"].val
pmdec_err = psr["PMDEC"].err
# From Reardon et al. (2016)
dkpc = 0.15679
sigd = 0.00025
ldeg = 253.39
sigl = 0
bdeg = -41.96
sigb = 0

# Define other useful constants 
M_sun = 1.98847542e+30  # kg
mass_function = (4*np.pi**2/sc.G)*(asini*sc.c)**3/(pb*86400)**2/M_sun
rad_to_mas = 180*3600*1000/np.pi
parsec_to_m = 3.08567758e+16
sec_per_year = 86400*365.2425
n_samples = len(posteriors["PBDOT"].values)

# Pulsar mass
mc_posterior = posteriors["M2"].values
i_posterior = posteriors["KIN"].values
f2_posterior = posteriors["F2"].values
f0_posterior = posteriors["F0"].values
f1_posterior = posteriors["F1"].values

mtot_posterior = np.sqrt(np.abs((np.power((np.multiply(mc_posterior, np.sin(i_posterior*np.pi/180))), 3)/mass_function)))
mp_posterior = mtot_posterior - mc_posterior
mp = np.mean(mp_posterior)
mp_err = np.std(mp_posterior)
f2 = np.mean(f2_posterior)
f2_err = np.std(f2_posterior)
f0 = np.mean(f0_posterior)
f0_err = np.std(f0_posterior)
f1 = np.mean(f1_posterior)
f1_err = np.std(f1_posterior)


if display:
    plt.hist(mp_posterior, bins=50)
    plt.xlabel("$M_p$ ($M_\odot$)")
    plt.title("$M_p={0}\pm{1}$ ($M_\odot$)".format(round(mp, 3), round(mp_err, 3)))
    plt.show()
    plt.savefig("Mp.pdf")    
    plt.close()

# Pulsar distance from Shklovskii effect
pbdot_posterior = posteriors["PBDOT"].values # observed
pbdot_grav = -3.2e-16

# Expected Shklovskii and Galactic terms
Ex_pl =  GalDynPsr.modelLb.Expl(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term parallel to the Galactic plane
Ex_z =  GalDynPsr.modelLb.Exz(ldeg, sigl, bdeg, sigb, dkpc, sigd) # excess term perpendicular to the Galactic plane
errpl = np.abs(0.03*Ex_pl) # 3% for circular Bovy et al. (2012)
errz = np.abs(0.1*Ex_z) # 10% for vertical Holmberg and Flynn (2004)
pbdot_gal = GalDynPsr.pdotint.PdotGal(Ex_pl, Ex_z, pb*86400)
pbdot_gal_err =GalDynPsr.pdotint.ErrPdotGal(Ex_pl, errpl, Ex_z, errz, pb*86400, pb_err*86400)

print("Observed Pbdot = ", np.mean(pbdot_posterior), " +/- ", np.std(pbdot_posterior))
print("Galactic Pbdot contribution = ", pbdot_gal, " +/- ", pbdot_gal_err) 
pbdot_gal_posterior = np.random.normal(loc=pbdot_gal, scale=pbdot_gal_err, size=n_samples)  # sample randomly from pbdot_gal distribution 
pbdot_shklovskii = np.subtract(pbdot_posterior, pbdot_gal_posterior) - pbdot_grav
print("Observed Shklovskii contribution = ", np.mean(pbdot_shklovskii), " +/- ", np.std(pbdot_shklovskii))
pm = np.sqrt(pmra**2 + pmdec**2)/(sec_per_year*rad_to_mas)
D_posterior = sc.c*pbdot_shklovskii/(pm**2*pb*86400)/parsec_to_m
D = np.mean(D_posterior)
D_err = np.std(D_posterior)

if display:
    plt.hist(pbdot_posterior, bins=50, label="Observed")
    plt.hist(pbdot_shklovskii, bins=50, label="Kinematic")
    plt.xlabel("$\dot{P}_b$")
    plt.legend(loc='upper center')
    plt.savefig("pbdot1.pdf")
    plt.close()
    plt.hist(D_posterior, bins=50)
    plt.xlabel("$D$ (pc)")
    plt.title("$D={0}\pm{1}$ pc".format(round(D, 2), round(D_err, 2)))
    plt.savefig("pbdot2.pdf")
    plt.close()
    
print("======================================")
print("FINAL NUMBERS")
print("======================================")
print("Pulsar mass = ", mp, " +/- ",  mp_err)
print("Pulsar distance = ", D, " +/- ",  D_err)
print("Orbital Inclination = ", np.mean(i_posterior), " +/- ", np.std(i_posterior))


print("======================================")

# Radial velocity from F2!
#f2 = -6.15114143991082e-29
p0_obs = 1/f0 
p0_int = p0_obs/(1 - 75000/sc.c) # approximate Doppler correction
p1_obs = -f1/f0**2 
p1_int = p1_obs - GalDynPsr.pdotint.PdotGal(Ex_pl, Ex_z, p0_int) # Shklovskii correction
p2_obs = f1*(2*f1/f0**3) - (1/f0**2)*(f2)
v_r = (2*p1_int*D*parsec_to_m - p2_obs*(sc.c/pm**2))/(3*p0_int)/1000 # km/s ... assuming p2_int=0
print("Observed P2 = ", p2_obs)
print("Radial velocity = ", v_r, " km/s")

# Make shapiro delay figure
"""
First, if shapiro.out does not exist, make it!
Use:
tempo2 -output general2 -f data/J0437-4715.temponest.par data/J0437-4715.new.tim -nobs 30000 -nofit -set M2 0 -s "{bat}\t{post}\t{err}\t{freq}\t{binphase}\n" > shapiro.out
with a par file that has TNSubtractRed 1 and  TNSubtractDM 1, and then cut off the leading junk to just leave the post-fit residuals
"""
data = np.loadtxt('shapiro.out')
post_res = data[:, 1]
toa_err = data[:, 2]/10**6
binphase = data[:, 4]
# bin the data
n_bins=80
phases = np.linspace(0, 1, n_bins)
dphase = phases[1] - phases[0]
average_residual = []
average_error = []
for iphase in phases:
    mask = ((binphase > (iphase - dphase)) & (binphase < (iphase + dphase)))
    weights = np.divide(1, np.power(toa_err[mask], 1))
    average_residual.append(np.average(post_res[mask], weights=weights)*10**6)
    average_error.append(np.sqrt(np.sum(np.multiply(weights, (post_res[mask] - average_residual[-1]/10**6)**2))/np.sum(weights))*10**6/4)
    #print(iphase, average_residual[-1], average_error[-1])
plt.errorbar(phases*2*np.pi, average_residual, yerr=average_error)
plt.ylabel("Residual (us)")
plt.xlabel("Binary phase (rad)")
plt.savefig("shapiro.pdf")
plt.close()
