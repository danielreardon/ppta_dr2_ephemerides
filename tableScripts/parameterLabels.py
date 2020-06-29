"""
 parameter labels for binary and solitarly pulsar tables
"""

def getParamLabels(): 

  parameterLabels = {
  'F0'      : r'Pulse frequency, $\nu$ (${\rm s}^{-1}$)',
  'F1'      : r'First time frequency derivative, ${\dot \nu}$ (${\rm s}^{-2}$)',
  'DM'      : r'Dispersion measure, DM (${\rm cm}^{-3}\,{\rm pc}$)',
  'PMRA'    : r'Proper motion in RA (${\rm mas}\,{\rm yr}^{-1}$)',
  'PMDEC'   : r'Proper motion in DEC (${\rm mas}\,{\rm yr}^{-1}$)',
  'PX'      : r'Parallax, $\pi$ (${\rm mas}$)',
  'PB'      : r'Obital period, $P_{\mathrm{b}}$ ($\mathrm{d}$)',
  'PBDOT'   : r'First time derivative of orbital period, ${\dot P}_{\mathrm{b}}$ ($10^{-12}$)',
  'A1'      : r'Projected semimajor axis, $x$ (lt-s)',
  'A1DOT'   : r'Rate of change of projected semi-major axis, ${\dot x}$ ($10^{-12}$)',
  'TASC'    : r'Epoch of ascending node, $T_{\mathrm{ASC}}$ (MJD)',
  'TO'      : r'Epoch of periastron, $T_0$ (MJD)',
  'OM'      : r'Longitude of periastron, $XX$ (deg)',
  'ECC'     : r'Eccentricity of orbit, $XX$',
  'ECCDOT'  : r'Rate of change of eccentricity, $XX$ (deg/yr)',
  'OMDOT'   : r'Rate of advance of periastron, $XX$ (deg/yr)',
  'EPS1'    : r'EPS1 ...',
  'EPS2'    : r'EPS2 ...',
  'EPS1DOT' : r'EPS1DOT ...',
  'EPS2DOT' : r'EPS2DOT ...',
  'SINI'    : r'Sine of inclination angle, $XX$',
  'M2'      : r'Companion mass, $M_{\mathrm{c}}$ ($M_{\odot}$)',
  'H3'      : r'H3 ...',
  'H4'      : r'H4 ...',
  'STIG'    : r'STIG ....',
  'KOM'     : r'KOM ...',
  'KIN'     : r'KIN ...'
  }


  return parameterLabels
