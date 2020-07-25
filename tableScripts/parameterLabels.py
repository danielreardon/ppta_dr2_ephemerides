"""
 parameter labels for binary and solitarly pulsar tables
"""

def getParamLabels():

  parameterLabels = {
  'NTOA'                 : r'Number of TOAs',
  'MJDRange'             : r'MJD range',
  'RAJ'                  : r'Right ascension (RA), $\alpha$ (hh:mm:ss)',
  'DECJ'                 : r'Declination (DEC), $\delta$ (dd:mm:ss)',
  'F0'                   : r'Spin frequency, $\nu$ (${\rm s}^{-1}$)',
  'F1'                   : r'First spin frequency derivative, ${\dot \nu}$ (${\rm s}^{-2}$)',
  'F2'                   : r'Second spin frequency derivative, ${\dotdot \nu}$ (${\rm s}^{-3}$)',
  'DM'                   : r'Dispersion measure, DM (${\rm cm}^{-3}\,{\rm pc}$)',
  'PMRA'                 : r'Proper motion in RA, $\mu_\alpha$\cos\delta (${\rm mas}\,{\rm yr}^{-1}$)',
  'PMDEC'                : r'Proper motion in DEC, $\mu_\delta$ (${\rm mas}\,{\rm yr}^{-1}$)',
  'PX'                   : r'Parallax, $\pi$ (${\rm mas}$)',
  'PB'                   : r'Obital period, $P_{\mathrm{b}}$ ($\mathrm{d}$)',
  'PBDOT'                : r'First time derivative of orbital period, ${\dot P}_{\mathrm{b}}$ ($10^{-12}$)',
  'A1'                   : r'Projected semimajor axis, $x$ (lt-s)',
  'XDOT'                 : r'Rate of change of projected semi-major axis, ${\dot x}$ ($10^{-12}$)',
  'TASC'                 : r'Epoch of ascending node, $T_{\mathrm{ASC}}$ (MJD)',
  'T0'                   : r'Epoch of periastron, $T_0$ (MJD)',
  'OM'                   : r'Longitude of periastron, $\omega$ (deg)',
  'ECC'                  : r'Eccentricity of orbit, $e$',
  'OMDOT'                : r'Rate of advance of periastron, ${\dot \omega}$ (deg/yr)',
  'EPS1'                 : r'First Laplace parameter, $\epsilon_1 = e \sin \omega$',
  'EPS2'                 : r'Second Laplace parameter, $\epsilon_2 = e \cos \omega$',
  'EPS1DOT'              : r'Time derivative of $\epsilon_1$, $\dot{\epsilon}_1$ (${\rm s}^{-1}$)',
  'EPS2DOT'              : r'Time derivative of $\epsilon_2$, $\dot{\epsilon}_2$ (${\rm s}^{-1}$)',
  'SINI'                 : r'Sine of inclination angle, $\sin i$',
  'M2'                   : r'Companion mass, $M_{\mathrm{c}}$ ($M_{\odot}$)',
  'H3'                   : r'Orthometric Shapiro delay parameter, $h_3$ ($\mu\,$s)',
  'H4'                   : r'Orthometric Shapiro delay parameter, $h_4$ ($\mu\,$s)',
  'STIG'                 : r'Orthometric Shapiro delay parameter, $\zeta = h_4 / h_3$',
  'KOM'                  : r'Longitude of ascending node, $\Omega$ (deg.)',
  'KIN'                  : r'Inclination angle, $i$ (deg.)',
  'ELAT'                 : r'Ecliptic latitude $\beta$ (deg.)',
  'ELONG'                : r'Ecliptic longitude $\lambda$ (deg.)',
  'PMELAT'               : r'Proper motion in ecliptic latitude, $\mu_\beta$ (${\rm mas}\,{\rm yr}^{-1}$)',
  'PMELONG'              : r'Proper motion in ecliptic longitude, $\mu_\lambda$ (${\rm mas}\,{\rm yr}^{-1}$)',
  'D_PX(med/16th/84th)'  : r'$D_\pi$ ($\pm$ central $68\%$ range)',
  'D_SHK(med/16th/84th)' : r'$D_\rm{shk}$ ($\pm$ central $68\%$ range)',
  'INC_LIM(med/std)'     : r'95\%-confidence upper limit on $i$ (deg.)',
  'INC(med/16th/84th)'   : r'Inclination angle, $i$ (deg.)',
  'M2(med/16th/84th)'    : r'Companion mass, $M_{\mathrm{c}}$ ($M_{\odot}$)',
  'MP(med/16th/84th)'    : r'Pulsar mass, $M_{\mathrm{P}}$ ($M_{\odot}$) ($\pm$ central $68\%$ range)',
  'MTOT(med/16th/84th)'  : r'Total mass $M_{\rm tot}$ ($M_{\odot}$) ($\pm$ central $68\%$ range)',
  'OMDOT_GR'             : r'$\dot{\omega}_{\rm{GR}}$(Fix name)'
  }


  return parameterLabels
