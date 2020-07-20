"""
 parameter labels for binary and solitarly pulsar tables
"""

def getParamLabels(): 

  parameterLabels = {
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
  'OM'                   : r'Longitude of periastron, $XX$ (deg)',
  'ECC'                  : r'Eccentricity of orbit, $XX$',
  'ECCDOT'               : r'Rate of change of eccentricity, $XX$ (deg/yr)',
  'OMDOT'                : r'Rate of advance of periastron, $XX$ (deg/yr)',
  'EPS1'                 : r'EPS1(Fix name)',
  'EPS2'                 : r'EPS2(Fix name)',
  'EPS1DOT'              : r'EPS1DOT(Fix name)',
  'EPS2DOT'              : r'EPS2DOT(Fix name)',
  'SINI'                 : r'Sine of inclination angle, $XX$',
  'M2'                   : r'Companion mass, $M_{\mathrm{c}}$ ($M_{\odot}$)',
  'H3'                   : r'H3(Fix name)',
  'H4'                   : r'H4(Fix name)',
  'STIG'                 : r'STIG(Fix name)',
  'KOM'                  : r'KOM(Fix name)',
  'KIN'                  : r'KIN(Fix name)',
  'ELAT'                 : r'Ecliptic latitude $XX$ (deg.)',
  'ELONG'                : r'Ecliptic longitude $XX$ (deg.)',
  'PMELAT'               : r'Proper motion in ecliptic latitude $XX$ (unit)',
  'PMELONG'              : r'Proper motion in ecliptic longitude $XX$ (unit)', 
  'MASS_FUNC'            : r'MASS\_FUNC(Fix name)', 
  'D_PX(med/16th/84th)'  : r'D\_PX(med/16/84)(Fix name)',
  'D_SHK(med/16th/84th)' : r'D\_SHK(med/16th/84th)(Fix name)',
  'INC_LIM(med/std)'     : r'INC\_LIM(med/std)(Fix name)',
  'INC(med/16th/84th)'   : r'INC(med/16th/84th)(Fix name)',
  'M2(med/16th/84th)'    : r'M2(med/16th/84th)(Fix name)',
  'MP(med/16th/84th)'    : r'MP(med/16th/84th)(Fix name)',
  'MTOT(med/16th/84th)'  : r'MTOT(med/16th/84th)(Fix name)',
  'OMDOT_GR'             : r'OMDOT\_GR(Fix name)'
  }


  return parameterLabels
