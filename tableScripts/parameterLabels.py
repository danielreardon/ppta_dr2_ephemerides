"""
 parameter labels for binary and solitarly pulsar tables
"""

def getParamLabels(): 

  parameterLabels = {
  'NTOA'                 : r'Number of TOAs',
  'MJDRange'             : r'MJD range',
  'RAJ'                  : r'Right ascension (RA), $\alpha$ (hh:mm:ss)',
  'DECJ'                 : r'Declination (DEC), $\delta$ (dd:mm:ss)',
  'F0'                   : r'Pulse frequency, $\nu$ (${\rm s}^{-1}$)',
  'F1'                   : r'First time frequency derivative, ${\dot \nu}$ (${\rm s}^{-2}$)',
  'DM'                   : r'Dispersion measure, DM (${\rm cm}^{-3}\,{\rm pc}$)',
  'PMRA'                 : r'Proper motion in RA (${\rm mas}\,{\rm yr}^{-1}$)',
  'PMDEC'                : r'Proper motion in DEC (${\rm mas}\,{\rm yr}^{-1}$)',
  'PX'                   : r'Parallax, $\pi$ (${\rm mas}$)',
  'PB'                   : r'Obital period, $P_{\mathrm{b}}$ ($\mathrm{d}$)',
  'PBDOT'                : r'First time derivative of orbital period, ${\dot P}_{\mathrm{b}}$ ($10^{-12}$)',
  'A1'                   : r'Projected semimajor axis, $x$ (lt-s)',
  'A1DOT'                : r'Rate of change of projected semi-major axis, ${\dot x}$ ($10^{-12}$)',
  'TASC'                 : r'Epoch of ascending node, $T_{\mathrm{ASC}}$ (MJD)',
  'TO'                   : r'Epoch of periastron, $T_0$ (MJD)',
  'OM'                   : r'Longitude of periastron, $\omega_0$ (deg)',
  'ECC'                  : r'Eccentricity of orbit, $e$',
  'ECCDOT'               : r'Rate of change of eccentricity, ${\dot e}$ (deg/yr)',
  'OMDOT'                : r'Rate of advance of periastron, ${\dot \omega_0}$ (deg/yr)',
  'EPS1'                 : r'EPS1, $e \sin \omega_0$',
  'EPS2'                 : r'EPS2, $e \cos \omega_0$',
  'EPS1DOT'              : r'First derivative of EPS1',
  'EPS2DOT'              : r'First derivative of EPS2',
  'SINI'                 : r'Sine of inclination angle, $\sin i$',
  'M2'                   : r'Companion mass, $M_{\mathrm{c}}$ ($M_{\odot}$)',
  'H3'                   : r'\textbf{fix} H3',
  'H4'                   : r'\textbf{fix} H4',
  'STIG'                 : r'\textbf{fix} STIG',
  'KOM'                  : r'\textbf{fix} KOM',
  'KIN'                  : r'\textbf{fix} KIN',
  'ELAT'                 : r'Ecliptic latitude $XX$ (deg.)',
  'ELONG'                : r'Ecliptic longitude $XX$ (deg.)',
  'PMELAT'               : r'Proper motion in ecliptic latitude $XX$ (unit)',
  'PMELONG'              : r'Proper motion in ecliptic longitude $XX$ (unit)', 
  'MASS_FUNC'            : r'\textbf{fix} MASS\_FUNC', 
  'D_PX(med/16th/84th)'  : r'\textbf{fix} D\_PX ($\pm$ central $68\%$ range)',
  'D_SHK(med/16th/84th)' : r'\textbf{fix} D\_SHK ($\pm$ central $68\%$ range)',
  'INC_LIM(med/std)'     : r'\textbf{fix} INC\_LIM(med/std)',
  'INC(med/16th/84th)'   : r'\textbf{fix} INC ($\pm$ central $68\%$ range)',
  'M2(med/16th/84th)'    : r'Companion mass $M_2$ ($M_{\odot}$) ($\pm$ central $68\%$ range)',
  'MP(med/16th/84th)'    : r'Pulsar mass, $M_{\mathrm{P}}$ ($M_{\odot}$) ($\pm$ central $68\%$ range)',
  'MTOT(med/16th/84th)'  : r'Total mass $M_{\rm tot}$ ($M_{\odot}$) ($\pm$ central $68\%$ range)',
  'OMDOT_GR'             : r'\textbf{fix} OMDOT\_GR(Fix name)'
  }


  return parameterLabels
