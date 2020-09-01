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
Rounds float to number of sig figs
"""
from math import log10, floor
def round_sig(x, sig=1, small_value=1.0e-9):
    if sig==1 and str(x*10**10)[0] == '1':
        sig+=1  # First digit, add a sig fig
    value = round(x, sig - int(floor(log10(max(abs(x), abs(small_value))))) - 1)
    return value




def read_general2(filename, header=False):
    """
    Reads a general2 output into a numpy array
    """
    with open(filename, "r") as file:
        data = []
        files = []
        for line in file:
            if not header and not ('Finish' in line):
                files.append([(line.split('\t')[0])])
                data.append([float(x) for x in
                             (line.replace('\n', '').split('\t')[1:])])
            elif 'Starting general2 plugin' in line:
                header = False
            elif 'Finish' in line:
                return np.array(data), files
    return np.array(data), files


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
            if (shortFormat.find("e0"))!=-1:
                shortFormat = shortFormat.replace("e0", "\\times 10^{")+"}"
            elif (shortFormat.find("e-0"))!=-1:
                shortFormat = shortFormat.replace("e-0", "e-")
                shortFormat = shortFormat.replace("e", "\\times 10^{")+"}"
            elif (shortFormat.find("e"))!=-1:
                shortFormat = shortFormat.replace("e", "\\times 10^{")+"}"
            else: pass

        except:
            shortFormat = p

        if fitOrDer == 0:
            # get bolding right for derived M2
            if par=='M2' or par=='KIN' and shortFormat.find("+")!=-1:
               table.write('\t & \t ${}$'.format(shortFormat))
            else:
               table.write('\t & \t $\\mathbf{{ {} }}$'.format(shortFormat))
        elif fitOrDer == 1:
           table.write('\t & \t ${}$'.format(shortFormat))

    table.write('\\\\ \n')
    table.close()

    return None


def formatDerivedParams(psrDerived,ipsr,par):

    """ 
    Formatting for derived params with ^{+...}_{-...}
    """

    parameter =  psrDerived[ipsr][par]
    high = float(psrDerived[ipsr][par+'_84th']) - float(psrDerived[ipsr][par])
    low  = float(psrDerived[ipsr][par]) - float(psrDerived[ipsr][par+'_16th'])
    high = round_sig(high)
    low = round_sig(low)
    if high%1 == 0:
         high = int(high)
    if low%1 == 0:
        low = int(low)
    parameter = float(parameter)

    # if the high and low uncertainties are equal print as value(error)
    if high==low:
        parameterToWrite = ufloat(parameter,high)
        return parameterToWrite
    else: pass 

    if high>2 and low>2:
        parameterToWrite =  '{0}^{{ +{1} }}_{{ -{2} }}'.format(round(parameter), high, low)
    else:
        digit = np.max([len(str(high).split('.')[-1]), len(str(low).split('.')[-1])])
        parameterToWrite =  '{0}^{{ +{1} }}_{{ -{2} }}'.format(round(parameter, int(digit)), high, low)
    return parameterToWrite




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
         table.write('\t & \t ${}$'.format(int(np.floor(start))))
      except:
         table.write('\t & \t none?')
      try:
         finish = psrDeets[i]['FINISH']
         table.write(' -- ${}$'.format(int(np.floor(finish))))
      except:
         table.write('--none?')


    table.write('\\\\ \n')
    table.close()

    return None



def writeSkyPos(psrDeets,tabFile,parLabels):
    """
    Writitng RA and DEC rows
    """

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
    fromParFiles  = (['F0','F1','F2','DM','PMRA','PMDEC','PX',\
                      'PB','A1','TASC',\
                      'PBDOT','XDOT',\
                      'T0','OM','OMDOT','ECC',\
                      'EPS1','EPS2','EPS1DOT','EPS2DOT',\
                      'M2','SINI',\
                      'H3','H4','STIG',\
                      'KOM','KIN'])

    # parameters form derived_parameters.txt
    fromDerivedParams = (['ELAT','ELONG','PMELONG','PMELAT','PMELONG','MASS_FUNC',\
                          'D_PX(med/16th/84th)', 'D_SHK(med/16th/84th)',\
                          'INC(med/16th/84th)', 'INC_LIM(med/std)',\
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
                   'F0', 'F1','F2',\
                   'DM',\
                   'PMRA', 'PMDEC', 'PMELAT', 'PMELONG', 'PX',\
                   'D_PX(med/16th/84th)'])


    elif solitaryOrBinary == 'binary':

        params = (['ELAT','ELONG',\
                   'F0', 'F1',\
                   'DM',\
                   'PMRA', 'PMDEC', 'PMELAT', 'PMELONG',
                   'PX',\
                   'PB', 'A1',\
                   'T0', 'OM', 'ECC',\
                   'TASC', 'EPS1', 'EPS2',\
                   'PBDOT', 'OMDOT', 'XDOT',\
                   'SINI', 'M2', 'H3', 'H4', 'STIG',\
                   'KOM', 'KIN',\
                   'INC_LIM(med/std)',\
                   'D_PX(med/16th/84th)', 'D_SHK(med/16th/84th)',\
                   'MP(med/16th/84th)', 'MTOT(med/16th/84th)'])

    return fittedOrDerived, params



# point to list of pulsars to include in the table
# (split binary psrs into groups as table too wide)
parser = argparse.ArgumentParser(description='which group of pulsars to make the table for')
parser.add_argument('--solOrBin', type=str, help='choose solitary or binary',default='binary')
parser.add_argument('--groupNo', type=int, help='integer (1,2,3,4) for group of psrs',default=1)
args = parser.parse_args()

solBin = str(args.solOrBin)
whichGroup = int(args.groupNo)

datadir = '/Users/dreardon/Dropbox/Git/'
#datadir = '/fred/oz002/hmiddleton/ppta_ephemeris/repositories'

for solBin in ['solitary', 'binary']:

    if solBin=='solitary':
        groups = [1,2]
    else:
        groups = [1,2,3,4]

    for whichGroup in groups:

        if solBin=='solitary':
            if whichGroup==1 or whichGroup==2:
                psrNames = np.genfromtxt('psrLists/psrListSolitary-group{}.list'.format(whichGroup),dtype=str)
                tableFile=datadir + '/ppta_dr2_ephemerides/localTexTables/solitaryTable-group{}.tex'.format(whichGroup)
            else:
                print('Error: need to choose a group from [1,2] for solitary pulsar table')
                exit()

        elif solBin=='binary':
            if whichGroup==1 or whichGroup==2 or whichGroup==3 or whichGroup==4:
                psrNames = np.genfromtxt('psrLists/psrListBinary-group{}.list'.format(whichGroup),dtype=str)
                tableFile=datadir + '/ppta_dr2_ephemerides/localTexTables/binaryTable-group{}.tex'.format(whichGroup)
            else:
                print('Error: need to choose a group from [1,2,3,4] for binary pulsar table')
                exit()



        # read in pulsar details from par files
        psrDetails = []
        for psr in psrNames:

            #if psr == 'J1713+0747' or psr == 'J1909-3744':
            #    parLoc = datadir + '/ppta_dr2_ephemerides/publish_collection/dr2e/{}.kop.par'.format(psr)
            #    outLoc = datadir + '/ppta_dr2_ephemerides/publish_collection/dr2e/output/{}.kop.par.out'.format(psr)
            #else:
            #    parLoc = datadir + '/ppta_dr2_ephemerides/publish_collection/dr2e/{}.par'.format(psr)
            #    outLoc = datadir + '/ppta_dr2_ephemerides/publish_collection/dr2e/output/{}.par.out'.format(psr)
            parLoc = datadir + '/ppta_dr2_ephemerides/final/tempo2/{}.par'.format(psr)
            outLoc = datadir + '/ppta_dr2_ephemerides/final/tempo2/{}.out'.format(psr)
            #parLoc = '/fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr2_boris/new_params_ver1/{}.par'.format(psr)

            try:
                psrPars = readParFile.read_par(parLoc)
            except FileNotFoundError:
                parLoc = parLoc.replace('dr2e','dr2')
                outLoc = outLoc.replace('dr2e','dr2')
                psrPars = readParFile.read_par(parLoc)
            data, files = read_general2(outLoc, header=True)
            psrPars['START'] = np.min(data[:,0])
            psrPars['FINISH'] = np.max(data[:,0])
            psrPars['NEPOCH'] = len(np.unique(files))
            psrDetails.append(psrPars)




        # read in derived pulsar details
        psrDerived = []
        for psr in psrNames:

            if psr == 'J1713+0747' or psr == 'J1909-3744':
                parname = '{}.kop.par'.format(psr)
            else:
                parname = '{}.par'.format(psr)

            psrDer = readParFile.get_derived_params(parname)
            psrDerived.append(psrDer)





        # a place to save the table
        # clearing file and making table top matter
        #tableFile='localTexTables/binaryTable-group{}.tex'.format(whichGroup)
        table = open(tableFile,'w')
        table.write("""
        \\begin{table}
        \\footnotesize
        \\begin{tabular}{llllllll}
        \\hline\\hline \\noalign{\\vskip 1.5mm}
        """)


        ## write pulsar name row
        table.write('Pulsar Name ')
        for p in psrNames:
          p = p.replace('-','$-$')
          table.write('\t & \t {}'.format(p))
        table.write(' \n \\\\ \\hline \\noalign{\\vskip 1.5mm} \n')
        table.close()


        # getting the parameter labels / names for the table headings
        parameterNames = parameterLabels.getParamLabels()


        ## write number of toas row
        names = [ psrDetails[i]['NTOA'] for i in range(len(psrNames)) ]
        writeLine(names,tableFile,parameterNames['NTOA'],1)

        names = [ psrDetails[i]['NEPOCH'] for i in range(len(psrNames)) ]
        writeLine(names,tableFile,parameterNames['NEPOCH'],1)

        ## write MJD range row
        writeMJDRange(psrDetails,tableFile,parameterNames['MJDRange'])


        # read in pulsar details from par files
        psrDetails = []
        for psr in psrNames:

            parLoc = datadir + '/ppta_dr2_ephemerides/publish_collection_refit/dr2/{}.par'.format(psr)
            outLoc = datadir + '/ppta_dr2_ephemerides/final/tempo2/{}.out'.format(psr)

            try:
                psrPars = readParFile.read_par(parLoc)
            except FileNotFoundError:
                parLoc = parLoc.replace('dr2e','dr2')
                outLoc = outLoc.replace('dr2e','dr2')
                psrPars = readParFile.read_par(parLoc)
            data, files = read_general2(outLoc, header=True)
            psrPars['START'] = np.min(data[:,0])
            psrPars['FINISH'] = np.max(data[:,0])
            psrPars['NEPOCH'] = len(np.unique(files))
            psrDetails.append(psrPars)

        ## write sky position row (RA & DEC)
        writeSkyPos(psrDetails,tableFile,parameterNames)

        # read in pulsar details from par files
        psrDetails = []
        for psr in psrNames:

            #if psr == 'J1713+0747' or psr == 'J1909-3744':
            #    parLoc = datadir + '/ppta_dr2_ephemerides/publish_collection/dr2e/{}.kop.par'.format(psr)
            #    outLoc = datadir + '/ppta_dr2_ephemerides/publish_collection/dr2e/output/{}.kop.par.out'.format(psr)
            #else:
            #    parLoc = datadir + '/ppta_dr2_ephemerides/publish_collection/dr2e/{}.par'.format(psr)
            #    outLoc = datadir + '/ppta_dr2_ephemerides/publish_collection/dr2e/output/{}.par.out'.format(psr)
            parLoc = datadir + '/ppta_dr2_ephemerides/final/tempo2/{}.par'.format(psr)
            outLoc = datadir + '/ppta_dr2_ephemerides/final/tempo2/{}.out'.format(psr)
            #parLoc = '/fred/oz002/dreardon/ppta_dr2_ephemerides/partim/dr2_boris/new_params_ver1/{}.par'.format(psr)

            try:
                psrPars = readParFile.read_par(parLoc)
            except FileNotFoundError:
                parLoc = parLoc.replace('dr2e','dr2')
                outLoc = outLoc.replace('dr2e','dr2')
                psrPars = readParFile.read_par(parLoc)
            data, files = read_general2(outLoc, header=True)
            psrPars['START'] = np.min(data[:,0])
            psrPars['FINISH'] = np.max(data[:,0])
            psrPars['NEPOCH'] = len(np.unique(files))
            psrDetails.append(psrPars)



        fittedOrDerived, params = get_parameters_for_table(solBin)

        # getting the parameter labels / names for the table headings
        parameterNames = parameterLabels.getParamLabels()


        # write parameters from list
        nparam = 5
        for ipar, par in enumerate(params):
          if par == 'MASS_FUNC':
              continue

          print ('\n ',par,type(par))

          paramList = []

          for ipsr, psr in enumerate(psrNames):

            print(fittedOrDerived[par])
            if fittedOrDerived[par]==0:
              try:
                parameter = ufloat(psrDetails[ipsr][par], psrDetails[ipsr][str(par+'_ERR')])
                paramList.append(parameter)
              except:
                # is there a derived parameter for M2?
                if par=='M2':
                  try:
                    paraString = formatDerivedParams(psrDerived,ipsr,'M2(med/16th/84th)')
                    paramList.append(paraString)
                  except:
                    paramList.append('-')
                elif par=='KIN':
                  try:
                    paraString = formatDerivedParams(psrDerived,ipsr,'INC(med/16th/84th)')
                    paramList.append(paraString)
                  except:
                    paramList.append('-')
                else:
                  paramList.append('-')

            else:
                try:
                  # check about including errors for these params and how many dp to quote here.
                  if par=='ELAT' or par=='ELONG' or par=='PMELAT' or par=='PMELONG' or par=='MASS_FUNC' or par=='OMDOT_GR':
                    parameter = psrDerived[ipsr][par]
                    paramList.append('{0:.5f}'.format(float(parameter)))
                  elif par=='INC_LIM(med/std)':
                    parameter = psrDerived[ipsr][par]
                    parameter = parameter.replace('<','')
                    parameter = float(parameter) #+ float(psrDerived[ipsr][par+'_16th'])
                    paramList.append('<{0}'.format(round(parameter)))
                  else:
                    ''' -> moved to function formatDerivedParams(...)
                    parameter =  psrDerived[ipsr][par]
                    #print('here ', psrDerived[ipsr][par+'_16th'])
                    high = float(psrDerived[ipsr][par+'_84th']) - float(psrDerived[ipsr][par])
                    low  = float(psrDerived[ipsr][par]) - float(psrDerived[ipsr][par+'_16th'])
                    high = round_sig(high)
                    low = round_sig(low)
                    if high%1 == 0:
                        high = int(high)
                    if low%1 == 0:
                        low = int(low)
                    parameter = float(parameter)
                    if high>2 and low>2:
                        parameterStr =  '{0}^{{ +{1} }}_{{ -{2} }}'.format(round(parameter), high, low)
                    else:
                        digit = np.max([len(str(high).split('.')[-1]), len(str(low).split('.')[-1])])
                        parameterStr =  '{0}^{{ +{1} }}_{{ -{2} }}'.format(round(parameter, int(digit)), high, low)
                    '''
                    parameterStr = formatDerivedParams(psrDerived,ipsr,par)
                    paramList.append(parameterStr)
                    print('parSTRING ', parameterStr)
                except:
                  paramList.append('-')

          # write parameter line
          if nparam % 5 == 0:
              table=open(tableFile,'a')
              table.write('\n \\noalign{\\vskip 1.5mm} \n')
              table.close()
          writeLine(paramList,tableFile,parameterNames[par],fittedOrDerived[par],parLabel=par)
          nparam += 1

        # end table stuff
        table=open(tableFile,'a')
        table.write("""
        \\noalign{\\vskip 1.5mm}
        \\hline\\hline
        \\end{tabular}\\hfill\\
        \\caption{\\label{tab:XXXXX}
        Placeholder caption.....
        }
        \\end{table}
        """)
        table.close()






