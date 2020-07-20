"""
Read in par files and make latex table for solitary pulsars
"""

import os, json
import numpy as np
from astropy import units as u
import astropy
from uncertainties import ufloat
import uncertainties
import decimal
import string

from astropy import units as u
from astropy.coordinates import SkyCoord

# read in par file - copied from Daniel's script
import readParFile
import parameterLabels



"""
Converts to short format undertainties  from:  
https://pythonhosted.org/uncertainties/user_guide.html 
"""
class ShorthandFormatter(string.Formatter):

    def format_field(self, value, format_spec):
        if isinstance(value, uncertainties.UFloat):
            return value.format(format_spec+'S')  # Shorthand option added
        # Special formatting for other types can be added here (floats, etc.)
        else:
            # Usual formatting:
            return super(ShorthandFormatter, self).format_field(
                value, format_spec)






"""
Write the parameter line in latex table format
"""
def writeLine(parameters,tableFile,parameterName):
    table = open(tableFile,'a')
    table.write(parameterName)

    frmtr = ShorthandFormatter()

    for p in parameters:

        try: 
            shortFormat = frmtr.format("{0:.1u}",p)

            # does the string contain "e" ?
            if (shortFormat.find("e"))!=-1:
                shortFormat = shortFormat.replace("e", "\\times 10^{")+"}"
            else: pass

        except:
            shortFormat = p


        table.write('\t & \t ${}$'.format(shortFormat))

    table.write('\\\\ \n')
    table.close()

    return None




"""
Writes the MJD range to the table
Maybe some missing from the par files
To do - check par files for start/finish
"""
def writeMJDRange(psrDeets,tabFile):

    # row title
    table = open(tabFile,'a')
    table.write('MJD range')

    # if can't get value, write none for now
    for i in range(len(psrDeets)):
    
      try: 
         start = psrDeets[i]['START']
         table.write('\t & \t ${:.1f}$'.format(start))
      except:
         table.write('\t & \t none?')
      try:
         finish = psrDeets[i]['FINISH']
         table.write('--${:.1f}$'.format(finish))
      except: 
         table.write('--none?')
    

    table.write('\\\\ \n')
    table.close()

    return None



"""
Writing the RA and DEC rows
"""
def writeSkyPos(psrDeets,tabFile):

    # row title
    table = open(tabFile,'a')

    ipsr = 0

    ras  = []
    decs = []

    for i in range(len(psrDeets)):

        # get sky position 
        pos = SkyCoord(psrDetails[ipsr]['RAJ'] + ' ' +str(psrDetails[ipsr]['DECJ']),unit=(u.hourangle,u.deg))

        # get RA as string with error
        raSecondsAndErr = ufloat(pos.ra.hms.s,psrDeets[ipsr]['RAJ_ERR'])
        frmtr = ShorthandFormatter()
        shortFormat = frmtr.format("{0:.1u}",raSecondsAndErr)
        ras.append('$'+str(int(pos.ra.hms.h))+'$:$'+str(int(pos.ra.hms.m))+'$:$'+str(shortFormat)+'$')

        # get DEC as string with error
        decSecondsAndErr = ufloat(abs(pos.dec.dms.s),psrDeets[ipsr]['DECJ_ERR'])
        frmtr = ShorthandFormatter()
        shortFormat = frmtr.format("{0:.1u}",decSecondsAndErr)
        decs.append('$'+str(int(pos.dec.dms.d))+'$:$'+str(int(abs(pos.dec.dms.m)))+'$:$'+str(shortFormat)+'$')


    # write ra
    table.write('Right ascension (RA), hh:mm:ss')
    for ra in ras:
        table.write('\t & \t {}'.format(ra))
    table.write('\\\\ \n')
  
    # write dec
    table.write('Declination (DEC), dd:mm:ss')
    for dec in decs:
        table.write('\t & \t {}'.format(dec))
    table.write('\\\\ \n')

    return None    



# point to list of pulsars to include in the table
# (split into two groups as table too wide)
psrNames = np.genfromtxt('psrLists/psrListSolitary-all.list',dtype=str)


# read in pulsar details
psrDetails = []
for psr in psrNames:

    parLoc = '/fred/oz002/hmiddleton/ppta_ephemeris/repositories/ppta_dr2_ephemerides/publish_collection/dr2/{}.par'.format(psr)
    #parLoc = '/fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr2_boris/new_params_ver1/{}.par'.format(psr)

    psrPars = readParFile.read_par(parLoc)
    psrDetails.append(psrPars)





# place to save the table
# clearing file and making table top matter
tableFile='localTexTables/solitary-table.tex'
table = open(tableFile,'w')
table.write("""
\\begin{table*}
\\begin{tabular}{llllllll}
\\hline\\hline \\\\\
""")

## write pulsar name row
table.write('Pulsar Name ')
for p in psrNames:
  table.write('\t & \t {}'.format(p))
table.write(' \\\\ \n \\\\ \\hline \\\\ \n')
table.close()


# parameters 

"""
## write number of toas row
names = [ psrDetails[i]['NTOA'] for i in range(len(psrNames)) ]
writeLine(names,tableFile,'Number of TOAs')


## write MJD range row
writeMJDRange(psrDetails,tableFile)

## write sky position row (RA & DEC)
writeSkyPos(psrDetails,tableFile)

"""
# parameters that can all be treated the same 

params = (['F0','F1','DM','PMRA','PMDEC','PX'])
parameterNames = parameterLabels.getParamLabels()
"""
parameterNames = ([r'Pulse frequency, $\nu$ (${\rm s}^{-1}$)',
                   r'First frequency derivative, $\nu$ (${\rm s}^{-2}$)',
                   r'Dispersion measure, DM (${\rm cm}^{-3}\,{\rm pc}$)',
                   r'Proper motion in RA (${\rm mas}\,{\rm yr}^{-1}$)',
                   r'Proper motion in DEC (${\rm mas}\,{\rm yr}^{-1}$)',
                   r'Parallax, $\pi$ (${\rm mas}$)'])
"""

# write parameters from list 
for ipar, par in enumerate(params):

  print ('\n ',par) 

  paramList = []

  for ipsr, psr in enumerate(psrNames):

      try: 
        parameter = ufloat(psrDetails[ipsr][par], psrDetails[ipsr][str(par+'_ERR')])
        paramList.append(parameter)
      except:
        print('no parameter!')
        # if doesn't exist, write 0.0 for now
        #paramList.append(ufloat(0.0,0.0))
        paramList.append('-')

  # write parameter line
  writeLine(paramList,tableFile,parameterNames[par])


# end table stuff
table=open(tableFile,'a')
table.write("""
\\\\ \\hline\\hline
\end{tabular}
\end{table*}
""")
table.close()






