"""
 parameter labels for binary and solitarly pulsar tables
"""

def getParamLabels():

  parameterLabels = {
  'NTOA'                 : r'Number of TOAs\dotfill',
  'NEPOCH'               : r'Number of observations\dotfill',
  'MJDRange'             : r'MJD range\dotfill',
  'RAJ'                  : r'Right ascension (RA), $\alpha$ (hh:mm:ss)\dotfill',
  'DECJ'                 : r'Declination (DEC), $\delta$ (dd:mm:ss)\dotfill',
  'F0'                   : r'Spin frequency, $f$ (${\rm s}^{-1}$)\dotfill',
  'F1'                   : r'First spin frequency derivative, ${\dot{f}}$ (${\rm s}^{-2}$)\dotfill',
  'F2'                   : r'Second spin frequency derivative, ${\ddot{f}}$ (${\rm s}^{-3}$)\dotfill',
  'DM'                   : r'Dispersion measure, DM (${\rm cm}^{-3}\,{\rm pc}$)\dotfill',
  'PMRA'                 : r'Proper motion in RA, $\mu_\alpha \cos\delta$ (${\rm mas}\,{\rm yr}^{-1}$)\dotfill',
  'PMDEC'                : r'Proper motion in DEC, $\mu_\delta$ (${\rm mas}\,{\rm yr}^{-1}$)\dotfill',
  'PX'                   : r'Parallax, $\pi$ (${\rm mas}$)\dotfill',
  'PB'                   : r'Obital period, $P_{\mathrm{b}}$ ($\mathrm{d}$)\dotfill',
  'PBDOT'                : r'First time derivative of orbital period, ${\dot P}_{\mathrm{b}}$ \dotfill',
  'A1'                   : r'Projected semimajor axis, $x$ (lt-s)\dotfill',
  'XDOT'                 : r'Rate of change of projected semi-major axis, ${\dot x}$ \dotfill',
  'TASC'                 : r'Epoch of ascending node, $T_{\mathrm{ASC}}$ (MJD)\dotfill',
  'T0'                   : r'Epoch of periastron, $T_0$ (MJD)\dotfill',
  'OM'                   : r'Longitude of periastron, $\omega$ (deg)\dotfill',
  'ECC'                  : r'Eccentricity of orbit, $e$\dotfill',
  'OMDOT'                : r'Rate of advance of periastron, ${\dot \omega}$ (deg\,yr$^{-1}$)\dotfill',
  'EPS1'                 : r'First Laplace parameter, $\epsilon_1 = e \sin \omega$\dotfill',
  'EPS2'                 : r'Second Laplace parameter, $\epsilon_2 = e \cos \omega$\dotfill',
  'EPS1DOT'              : r'Time derivative of $\epsilon_1$, $\dot{\epsilon}_1$ (${\rm s}^{-1}$)\dotfill',
  'EPS2DOT'              : r'Time derivative of $\epsilon_2$, $\dot{\epsilon}_2$ (${\rm s}^{-1}$)\dotfill',
  'SINI'                 : r'Sine of inclination angle, $\sin i$\dotfill',
  'M2'                   : r'Companion mass, $M_{\mathrm{c}}$ ($M_{\odot}$)\dotfill',
  'H3'                   : r'Orthometric Shapiro delay parameter, $h_3$ ($\mu\,$s)\dotfill',
  'H4'                   : r'Orthometric Shapiro delay parameter, $h_4$ ($\mu\,$s)\dotfill',
  'STIG'                 : r'Orthometric Shapiro delay parameter, $\zeta = h_4 / h_3$\dotfill',
  'KOM'                  : r'Longitude of ascending node, $\Omega$ (deg.)\dotfill',
  'KIN'                  : r'Inclination angle, $i$ (deg.)\dotfill',
  'ELAT'                 : r'Ecliptic latitude $\beta$ (deg.)\dotfill',
  'ELONG'                : r'Ecliptic longitude $\lambda$ (deg.)\dotfill',
  'PMELAT'               : r'Proper motion in ecliptic latitude, $\mu_\beta$ (${\rm mas}\,{\rm yr}^{-1}$)\dotfill',
  'PMELONG'              : r'Proper motion in ecliptic longitude, $\mu_\lambda \cos\beta$ (${\rm mas}\,{\rm yr}^{-1}$)\dotfill',
  'PX_LKB(med/16th/84th)': r'Parallax (L-K bias corrected), $\pi_{lk}$ (${\rm mas}$)\dotfill',
  'D_LKB(med/16th/84th)' : r'Parallax distance (L-K bias corrected), $D_{\pi,lk}$ (kpc)\dotfill',
  'D_PX(med/16th/84th)'  : r'Parallax distance, $D_\pi$ (kpc)\dotfill',
  'D_SHK(med/16th/84th)' : r'Shklovskii distance $D_{\rm shk}$ (kpc)\dotfill',
  'INC_LIM(med/std)'     : r'Upper limit on $i$ (deg.)\dotfill',
  'INC(med/16th/84th)'   : r'Inclination angle, $i$ (deg.)\dotfill',
  'M2(med/16th/84th)'    : r'Companion mass, $M_{\mathrm{c}}$ ($M_{\odot}$)\dotfill',
  'MP(med/16th/84th)'    : r'Pulsar mass, $M_{\mathrm{P}}$ ($M_{\odot}$) \dotfill',
  'MTOT(med/16th/84th)'  : r'Total system mass $M_{\rm tot}$ ($M_{\odot}$)\dotfill',
  'OMDOT_GR'             : r'$\dot{\omega}_{\rm{GR}}$(Fix name)\dotfill'
  }


  return parameterLabels
