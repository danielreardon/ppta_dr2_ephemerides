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
import argparse

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
def writeLine(parameters,tableFile,parameterName,fitOrDer,parLabel=None):
    table = open(tableFile,'a')
    table.write(parameterName)

    frmtr = ShorthandFormatter()

    # fix for PMELAT, PMELONG, ELAT, ELONG, PMRA, PMDEC
    if parLabel == 'PMELAT' or parLabel == 'PMELONG' or parLabel == 'ELAT' or parLabel == 'ELONG':
        fitOrDer = 0
    else: pass
    if parLabel == 'PMRA' or parLabel == 'PMDEC':
        fitOrDer = 1
    else: pass



    for p in parameters:

        try: 
            shortFormat = frmtr.format("{0:.1u}",p)

            # does the string contain "e" ?
            if (shortFormat.find("e"))!=-1:
                shortFormat = shortFormat.replace("e", "\\times 10^{")+"}"
            else: pass

        except:
            shortFormat = p

        if fitOrDer == 0:
           table.write('\t & \t $\\mathbf{{ {} }}$'.format(shortFormat))
        elif fitOrDer == 1:
           table.write('\t & \t ${}$'.format(shortFormat))

    table.write('\\\\ \n')
    table.close()

    return None




"""
Writes the MJD range to the table
Maybe some missing from the par files
To do - check par files for start/finish
"""
def writeMJDRange(psrDeets,tabFile,label):

    # row title
    table = open(tabFile,'a')
    table.write(label)

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
def writeSkyPos(psrDeets,tabFile,parLabels):

    # row title
    table = open(tabFile,'a')

    ipsr = 0

    ras  = []
    decs = []

    for ipsr in range(len(psrDeets)):

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

        # move on to next pulsar
        #ipsr+=1


    # write ra
    table.write(parLabels['RAJ'])
    for ra in ras:
        table.write('\t & \t {}'.format(ra))
    table.write('\\\\ \n')
  
    # write dec
    table.write(parLabels['DECJ'])
    for dec in decs:
        table.write('\t & \t {}'.format(dec))
    table.write('\\\\ \n')

    return None    




def get_parameters_for_table(solitaryOrBinary):


    # paratmeters from the par files
    fromParFiles  = (['F0','F1','DM','PMRA','PMDEC','PX',\
                      'PB','A1','TASC',\
                      'PBDOT','A1DOT',\
                      'TO','OM','OMDOT','ECC','ECCDOT',\
                      'EPS1','EPS2','EPS1DOT','EPS2DOT',\
                      'M2','SINI',\
                      'H3','H4','STIG',\
                      'KOM','KIN'])
 
    # parameters form derived_parameters.txt 
    fromDerivedParams = (['ELAT','ELONG','PMELONG','PMELAT','PMELONG','MASS_FUNC',\
                          'D_PX(med/16th/84th)', 'D_SHK(med/16th/84th)',\
                          'INC(med/16th/84th)',\
                          'M2(med/16th/84th)', 'MP(med/16th/84th)', 'MTOT(med/16th/84th)',\
                          'OMDOT_GR'])

    # keeping track of which params are from where 
    fittedOrDerived = {}
    for p in fromParFiles:
        fittedOrDerived[p] = 0
    for p in fromDerivedParams:
        fittedOrDerived[p] = 1


    # parameters in the order that we want them to appear in the table
    # some tricky parameters to think about: INC_Lim (automatically display <) and M2 (can be fitted or derived?)

    if solitaryOrBinary == 'solitary':

        params = (['ELAT','ELONG',\
                   'F0', 'F1',\
                   'DM',\
                   'PMRA', 'PMDEC', 'PMELAT', 'PMELONG',\
                   'D_PX(med/16th/84th)'])


    elif solitaryOrBinary == 'binary':

        params = (['ELAT','ELONG',\
                   'F0', 'F1',\
                   'DM',\
                   'PMRA', 'PMDEC', 'PMELAT', 'PMELONG',\
                   'PB', 'PBDOT', 'A1', 'A1DOT', 'TASC',\
                   'TO', 'OM', 'OMDOT', 'OMDOT_GR',\
                   'ECC', 'ECCDOT', \
                   'EPS1', 'EPS2', 'EPS1DOT', 'EPS2DOT',\
                   'SINI',\
                   'D_PX(med/16th/84th)', 'D_SHK(med/16th/84th)',\
                   'INC(med/16th/84th)',\
                   'MASS_FUNC',\
                   'MP(med/16th/84th)', 'MTOT(med/16th/84th)',\
                   'H3', 'H4', 'STIG', \
                   'KOM', 'KIN' ])

    return fittedOrDerived, params



# point to list of pulsars to include in the table
# (split binary psrs into groups as table too wide)
parser = argparse.ArgumentParser(description='which group of pulsars to make the table for')
parser.add_argument('--solOrBin', type=str, help='choose solitary or binary')
parser.add_argument('--groupNo', type=int, help='integer (1,2,3,4) for group of psrs',default=0)
args = parser.parse_args()

solBin = str(args.solOrBin)
whichGroup = int(args.groupNo)

if solBin=='solitary':
    if whichGroup==1 or whichGroup==2:
        psrNames = np.genfromtxt('psrLists/psrListSolitary-group{}.list'.format(whichGroup),dtype=str)
        tableFile='localTexTables/solitaryTable-group{}.tex'.format(whichGroup)
    else: 
        print('Error: need to choose a group from [1,2] for solitary pulsar table')
        exit()

elif solBin=='binary':
    if whichGroup==1 or whichGroup==2 or whichGroup==3 or whichGroup==4:
        psrNames = np.genfromtxt('psrLists/psrListBinary-group{}.list'.format(whichGroup),dtype=str)
        tableFile='localTexTables/binaryTable-group{}.tex'.format(whichGroup)
    else: 
        print('Error: need to choose a group from [1,2,3,4] for binary pulsar table')
        exit()




# read in pulsar details from par files
psrDetails = []
for psr in psrNames:

    if psr == 'J1713+0747' or psr == 'J1909-3744':
        parLoc = '/fred/oz002/hmiddleton/ppta_ephemeris/repositories/ppta_dr2_ephemerides/publish_collection/dr2/{}.kop.par'.format(psr)
    else: 
        parLoc = '/fred/oz002/hmiddleton/ppta_ephemeris/repositories/ppta_dr2_ephemerides/publish_collection/dr2/{}.par'.format(psr)
    #parLoc = '/fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr2_boris/new_params_ver1/{}.par'.format(psr)

    psrPars = readParFile.read_par(parLoc)
    psrDetails.append(psrPars)




# read in derived pulsar details
psrDerived = []
for psr in psrNames:

    if psr == 'J1713+0747' or psr == 'J1909-3744':
        parName = '{}.kop.par'.format(psr)
    else: 
        parname = '{}.par'.format(psr)

    psrDer = readParFile.get_derived_params(parname)
    psrDerived.append(psrDer)





# a place to save the table
# clearing file and making table top matter
#tableFile='localTexTables/binaryTable-group{}.tex'.format(whichGroup)
table = open(tableFile,'w')
table.write("""
\\begin{table*}
\\begin{adjustbox}{angle=90}
\\footnotesize
\\begin{tabular}{llllllll}
\\hline\\hline \\\\\
""")


## write pulsar name row
table.write('Pulsar Name ')
for p in psrNames:
  table.write('\t & \t {}'.format(p))
table.write(' \\\\ \n \\\\ \\hline \\\\ \n')
table.close()


# getting the parameter labels / names for the table headings
parameterNames = parameterLabels.getParamLabels()


## write number of toas row
names = [ psrDetails[i]['NTOA'] for i in range(len(psrNames)) ]
writeLine(names,tableFile,parameterNames['NTOA'],1)

## write MJD range row
writeMJDRange(psrDetails,tableFile,parameterNames['MJDRange'])

## write sky position row (RA & DEC)
writeSkyPos(psrDetails,tableFile,parameterNames)



fittedOrDerived, params = get_parameters_for_table(solBin)

# getting the parameter labels / names for the table headings
parameterNames = parameterLabels.getParamLabels()


# write parameters from list 
for ipar, par in enumerate(params):

  print ('\n ',par,type(par)) 

  paramList = []

  for ipsr, psr in enumerate(psrNames):

    print(fittedOrDerived[par])    
    if fittedOrDerived[par]==0: 
      try: 
        parameter = ufloat(psrDetails[ipsr][par], psrDetails[ipsr][str(par+'_ERR')])
        paramList.append(parameter)
      except:
        """
        # dealing with different names - not needed I think  
        if par == 'ECC':
          try: 
            parE = 'E'
            parameter = ufloat(psrDetails[ipsr][parE], psrDetails[ipsr][str(parE+'_ERR')])
            paramList.append(parameter)
            print ('got E')
          except:
            pass
        elif par == 'ECCDOT': 
          try:
            parEDOT = 'EDOT'
            parameter = ufloat(psrDetails[ipsr][parEDOT], psrDetails[ipsr][str(parEDOT+'_ERR')])
            paramList.append(parameter)
            print ('got EDOT')
          except:
            pass
        elif par == 'A1DOT':
          try: 
            parXDOT = 'XDOT'
            parameter = ufloat(psrDetails[ipsr][parXDOT], psrDetails[ipsr][str(parXDOT+'_ERR')])
            paramList.append(parameter)
            print ('got XDOT')
          except:
            pass
        else:
        """
        print('no parameter!')
        paramList.append('-')

    else:
       try: 
         parameter =  psrDerived[ipsr][par]
         #print('here ', psrDerived[ipsr][par+'_16th'])
         high = float(psrDerived[ipsr][par+'_84th']) - float(psrDerived[ipsr][par])
         low  = float(psrDerived[ipsr][par]) - float(psrDerived[ipsr][par+'_16th'])
         parameterStr =  '{0:.3}^{{ +{1:.3} }}_{{ -{2:.3} }}'.format(float(psrDerived[ipsr][par]), high, low)
         paramList.append(parameterStr)
         print('parSTRING ', parameterStr)
       except: 
         paramList.append('-')

  # write parameter line
  writeLine(paramList,tableFile,parameterNames[par],fittedOrDerived[par],parLabel=par)


# end table stuff
table=open(tableFile,'a')
table.write("""
\\\\ \\hline\\hline
\\end{tabular}\\hfill\\
\\end{adjustbox}
\\caption{\\label{tab:XXXXX}
Placeholder caption.....
}
\\end{table*}
""")
table.close()






