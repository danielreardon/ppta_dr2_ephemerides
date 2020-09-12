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



def replaceEWithTimes10(value):

    # does the string contain "e" ?
    if (value.find("e0"))!=-1:
        value = value.replace("e0", "\\times 10^{")+"}"
    elif (value.find("e-0"))!=-1:
        value = value.replace("e-0", "e-")
        value = value.replace("e", "\\times 10^{")+"}"
    elif (value.find("e+0"))!=-1:
        value = value.replace("e+0","e+")
        value = value.replace("e", "\\times 10^{")+"}"
    elif (value.find("e"))!=-1:
        value = value.replace("e", "\\times 10^{")+"}"
    else: pass

    return value





def writeLine(parameters,tableFile,parameterName,keepTrackFitOrDer,parLabel=None):

    """
    Write the parameter line in latex table format
    will put fitted params in bold and params from the derived_param file not-bold
    """


    table = open(tableFile,'a')
    table.write(parameterName)

    frmtr = ShorthandFormatter()

    """
    no longer needed (I think...)
    # fix for PMELAT, PMELONG, ELAT, ELONG, PMRA, PMDEC
    if parLabel == 'PMELAT' or parLabel == 'PMELONG' or parLabel == 'ELAT' or parLabel == 'ELONG':
        fitOrDer = 0
    else: pass
    if parLabel == 'PMRA' or parLabel == 'PMDEC':
        fitOrDer = 1
    else: pass
    """

    # loop over the parameter values for each pulsar in this batch
    # using "keepTrackFitOrDer" to remember whether the parameters came from the
    # par files or the derived params file
    for i,p in enumerate(parameters):

        try:
            #does the error start with a 1?
            errorRounded = round_sig(p.std_dev)
            if str(p.std_dev*10**10)[0] == '1' and errorRounded!=0:
                p.std_dev = errorRounded
                shortFormat = frmtr.format("{0}",p)
            else: 
                shortFormat = frmtr.format("{0:.1u}",p)
            #old version shortFormat = frmtr.format("{0:.1u}",p)

            # does the string contain "e" ?
            shortFormat = replaceEWithTimes10(shortFormat)
        except:
            shortFormat = p
  
        # fitted in bold
        if keepTrackFitOrDer[i]=='f':
            table.write('\t & \t $\\mathbf{{ {} }}$'.format(shortFormat))

        # derived not bold
        elif keepTrackFitOrDer[i]=='d':
            table.write('\t & \t ${}$'.format(shortFormat))

        # the rest, probably writing '-' for no parameter
        else:
            table.write('\t & \t ${}$'.format(shortFormat))

    table.write('\\\\ \n')
    table.close()

    return None


def formatDerivedParams(psrDerived,ipsr,par):

    """
    Formatting for derived params with ^{+...}_{-...}
    EXCEPT for ECC(med/std) which is written as value(error)
    These are all not bold..
    """

    parameter =  psrDerived[ipsr][par]

    if par=='ECC(med/std)':
        print(psrDerived[ipsr][par],psrDerived[ipsr][par+'_ERR'])


    high = float(psrDerived[ipsr][par+'_84th']) - float(psrDerived[ipsr][par])
    low  = float(psrDerived[ipsr][par]) - float(psrDerived[ipsr][par+'_16th'])
    #print(high,low)
    high = round_sig(high)
    low = round_sig(low)
    #print(high,low)
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
        value = replaceEWithTimes10(str(round(parameter,int(digit))))
        low   = replaceEWithTimes10(str(low))
        high  = replaceEWithTimes10(str(high))

        parameterToWrite = '{{ {0} }} ^{{ +{1} }}_{{ -{2} }}'.format(value,high,low)

        #print(parameterToWrite)

    return parameterToWrite





def writeMJDRange(psrDeets,tabFile,label):

    """
    Writes the MJD range to the table
    """

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



def writeSkyPos(psrDets,tabFile,parLabels):
    """
    Writitng RA and DEC rows
    """

    # row title
    table = open(tabFile,'a')

    ipsr = 0

    ras  = []
    decs = []
    pmras  = []
    pmdecs = []

    for ipsr in range(len(psrDets)):

        # get sky position
        pos = SkyCoord(psrDets[ipsr]['RAJ'] + ' ' +str(psrDets[ipsr]['DECJ']),unit=(u.hourangle,u.deg))

        # get RA as string with error
        raSecondsAndErr = ufloat(pos.ra.hms.s,psrDets[ipsr]['RAJ_ERR'])
        frmtr = ShorthandFormatter()
        shortFormat = frmtr.format("{0:.1u}",raSecondsAndErr)
        ras.append('$'+str(int(pos.ra.hms.h))+'$:$'+str(int(pos.ra.hms.m))+'$:$'+str(shortFormat)+'$')

        # get DEC as string with error
        decSecondsAndErr = ufloat(abs(pos.dec.dms.s),psrDets[ipsr]['DECJ_ERR'])
        frmtr = ShorthandFormatter()
        shortFormat = frmtr.format("{0:.1u}",decSecondsAndErr)
        decs.append('$'+str(int(pos.dec.dms.d))+'$:$'+str(int(abs(pos.dec.dms.m)))+'$:$'+str(shortFormat)+'$')

        # PMRA and PMDEC - old version
        #pmras.append( frmtr.format("{0:.1u}",ufloat(psrDets[ipsr]['PMRA'], psrDets[ipsr]['PMRA_ERR'])))
        #pmdecs.append(frmtr.format("{0:.1u}",ufloat(psrDets[ipsr]['PMDEC'],psrDets[ipsr]['PMDEC_ERR'])))
        pmras.append(ufloat(psrDets[ipsr]['PMRA'], psrDets[ipsr]['PMRA_ERR']))
        pmdecs.append(ufloat(psrDets[ipsr]['PMDEC'],psrDets[ipsr]['PMDEC_ERR']))
  

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
    table.close()
   
    fitOrDerived = ['f' for i in range(len(pmras))]
    writeLine(pmras,tabFile,parLabels['PMRA'],fitOrDerived,parLabel='PMRA')
    writeLine(pmdecs,tabFile,parLabels['PMDEC'],fitOrDerived,parLabel='PMDEC')
  
    return None




def get_parameters_for_table(solitaryOrBinary):


    # paratmeters from the par files
    fromParFiles  = (['F0','F1','F2','DM','ELAT','ELONG','PMRA','PMDEC','PMELAT','PMELONG','PX',\
                      'PB','A1','TASC',\
                      'PBDOT','XDOT',\
                      'T0','OM','OMDOT','ECC',\
                      'EPS1','EPS2','EPS1DOT','EPS2DOT',\
                      'M2','SINI',\
                      'H3','H4','STIG',\
                      'KOM','KIN'])

    # parameters form derived_parameters.txt
    #fromDerivedParams = (['ELAT','ELONG','PMELONG','PMELAT','PMELONG','MASS_FUNC',\
    fromDerivedParams = (['MASS_FUNC',\
                          'ECC(med/std)', 'OM(med/std)', 'T0(med/std)',\
                          'D_PX(med/16th/84th)', 'D_SHK(med/16th/84th)',\
                          'INC(med/16th/84th)', 'INC_LIM(med/std)',\
                          'M2(med/16th/84th)', 'MP(med/16th/84th)', 'MTOT(med/16th/84th)',\
                          'MTOT_GR(med/16th/84th)',\
                          'OMDOT_GR', 'VT(med/16th/84th)'])

    # keeping track of which params are from where
    fittedOrDerived = {}
    for p in fromParFiles:
        fittedOrDerived[p] = 0
    for p in fromDerivedParams:
        fittedOrDerived[p] = 1


    # parameters in the order that we want them to appear in the table
    # some tricky parameters to think about: INC_Lim (automatically display <) and M2 (can be fitted or derived?)

    if solitaryOrBinary == 'solitary':

        params = (['F0', 'F1',\
                   'DM',\
                   'ELAT','ELONG',\
                   'PMELAT', 'PMELONG', 'PX',\
                   'D_PX(med/16th/84th)', 'VT(med/16th/84th)'])


    elif solitaryOrBinary == 'binary':

        params = (['F0', 'F1',\
                   'DM',\
                   'ELAT','ELONG',\
                   'PMELAT', 'PMELONG',
                   'PX',\
                   'PB', 'A1',\
                   'T0', 'OM', 'ECC',\
                   'TASC', 'EPS1', 'EPS2',\
                   'PBDOT', 'OMDOT', 'XDOT',\
                   'SINI', 'M2', 'H3', 'H4', 'STIG',\
                   'KOM', 'KIN',\
                   'INC_LIM(med/std)',\
                   'D_PX(med/16th/84th)', 'D_SHK(med/16th/84th)',\
                   'MP(med/16th/84th)', 'MTOT(med/16th/84th)', 'MTOT_GR(med/16th/84th)', 'VT(med/16th/84th)'])

    return fittedOrDerived, params



# point to list of pulsars to include in the table
# (split binary psrs into groups as table too wide)
parser = argparse.ArgumentParser(description='which group of pulsars to make the table for')
parser.add_argument('--solOrBin', type=str, help='choose solitary or binary',default='binary')
parser.add_argument('--groupNo', type=int, help='integer (1,2,3,4) for group of psrs',default=1)
args = parser.parse_args()

solBin = str(args.solOrBin)
whichGroup = int(args.groupNo)

#datadir = '/Users/dreardon/Dropbox/Git/'
datadir = '/fred/oz002/hmiddleton/ppta_ephemeris/repositories'

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
                #parname = '{}.kop.par'.format(psr)
                parname = '{}.par'.format(psr)
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
        print('Number of TOAS')
        names = [ psrDetails[i]['NTOA'] for i in range(len(psrNames)) ]
        noTOABolding = ['n' for i in range(len(psrNames))]
        writeLine(names,tableFile,parameterNames['NTOA'],noTOABolding)

        ## write number of observations row
        print('Number of observations')
        names = [ psrDetails[i]['NEPOCH'] for i in range(len(psrNames)) ]
        noObservationsBolding = ['n' for i in range(len(psrNames))]
        writeLine(names,tableFile,parameterNames['NEPOCH'],noObservationsBolding)

        ## write MJD range row
        print('MJD range')
        writeMJDRange(psrDetails,tableFile,parameterNames['MJDRange'])


        # read in pulsar details from par files
        psrDetails = []
        for psr in psrNames:

            # Use non-ecliptic for sky positions
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
        print('sky position (RA,DEC)')
        writeSkyPos(psrDetails,tableFile,parameterNames)



        # Now go back to ecliptic
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

        """   ###### I think this just repeats the lines above #####
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
        """

        # this is used for J1713  for the ecliptic parameters only
        parLocJ1713Ecliptic = datadir + '/ppta_dr2_ephemerides/publish_collection_refit/dr2e/ecliptic/J1713+0747.kop_ecliptic.par'
        J1713EclipticParameters=readParFile.read_par(parLocJ1713Ecliptic)



        fittedOrDerived, params = get_parameters_for_table(solBin)

        # getting the parameter labels / names for the table headings
        parameterNames = parameterLabels.getParamLabels()


        # write parameters from list
        nparam = 5
        for ipar, par in enumerate(params):
          if par == 'MASS_FUNC':
              continue

          print ('\n ',par,type(par))
          #print('\t\t ', fittedOrDerived[par])
          paramList = []


          # keep track of where parameters came from with values:
          # n: 'neither', f:fitted, d:derived
          keepingTrackFitDerived = ['n' for i in range(len(psrNames))]


          for ipsr, psr in enumerate(psrNames):

            #print(fittedOrDerived[par])
            if fittedOrDerived[par]==0:
              try:
                if par!='DM':
                  parameter = ufloat(psrDetails[ipsr][par], psrDetails[ipsr][str(par+'_ERR')])
                  paramList.append(parameter)
                  keepingTrackFitDerived[ipsr] = 'f'
                else: 
                  paramList.append(psrDetails[ipsr][par])
                  
                """
                parameter = ufloat(psrDetails[ipsr][par], psrDetails[ipsr][str(par+'_ERR')])
                paramList.append(parameter)
                keepingTrackFitDerived[ipsr] = 'f'
                """
              except:
                # is there a derived parameter for some parameters and gets ecliptic values for J1713 only
                if par=='M2':
                  try:
                    paraString = formatDerivedParams(psrDerived,ipsr,'M2(med/16th/84th)')
                    paramList.append(paraString)
                    keepingTrackFitDerived[ipsr] = 'd'
                  except:
                    paramList.append('-')
                elif par=='KIN':
                  try:
                    paraString = formatDerivedParams(psrDerived,ipsr,'INC(med/16th/84th)')
                    paramList.append(paraString)
                    keepingTrackFitDerived[ipsr] = 'd'
                  except:
                    paramList.append('-')
                elif par=='ECC':
                  try:
                    paraString = ufloat(psrDerived[ipsr]['ECC(med/std)'],psrDerived[ipsr]['ECC(med/std)_ERR'])
                    paramList.append(paraString)
                    keepingTrackFitDerived[ipsr] = 'd'
                  except:
                    paramList.append('-')
                elif par=='OM':
                  try:
                    paraString = ufloat(psrDerived[ipsr]['OM(med/std)'],psrDerived[ipsr]['OM(med/std)_ERR'])
                    paramList.append(paraString)
                    keepingTrackFitDerived[ipsr] = 'd'
                  except:
                    paramList.append('-')
                elif par=='T0':
                  try:
                    paraString = ufloat(psrDerived[ipsr]['T0(med/std)'],psrDerived[ipsr]['T0(med/std)_ERR'])
                    paramList.append(paraString)
                    keepingTrackFitDerived[ipsr] = 'd'
                  except:
                    paramList.append('-')
                elif    par=='ELAT'    and psr=='J1713+0747' \
                     or par=='ELONG'   and psr=='J1713+0747' \
                     or par=='PMELAT'  and psr=='J1713+0747' \
                     or par=='PMELONG' and psr=='J1713+0747':
                  try:
                    paraString = ufloat(J1713EclipticParameters[par],J1713EclipticParameters[par+'_ERR']) 
                    paramList.append(paraString)
                    keepingTrackFitDerived[ipsr] = 'f'
                  except:
                    paramList.append('-')
                else:
                  paramList.append('-')

            else:
                try:
                  # check about including errors for these params and how many dp to quote here.
                  # is this first if used? 
                  if par=='ELAT' or par=='ELONG' or par=='PMELAT' or par=='PMELONG' or par=='MASS_FUNC' or par=='OMDOT_GR':
                    print('test ', ipsr,psr,par,psrDerived[ipsr][par])
                    parameter = psrDerived[ipsr][par]
                    paramList.append('{0:.5f}'.format(float(parameter)))
                    keepingTrackFitDerived[ipsr] = 'd'
                  elif par=='INC_LIM(med/std)':
                    parameter = psrDerived[ipsr][par]
                    parameter = parameter.replace('<','')
                    parameter = float(parameter) #+ float(psrDerived[ipsr][par+'_16th'])
                    paramList.append('<{0}'.format(round(parameter)))
                    keepingTrackFitDerived[ipsr] = 'd'
                  else:
                    parameterStr = formatDerivedParams(psrDerived,ipsr,par)
                    paramList.append(parameterStr)
                    #print('parSTRING ', parameterStr)
                    keepingTrackFitDerived[ipsr] = 'd'
                except:
                  paramList.append('-')

          #print(keepingTrackFitDerived)

          # write parameter line
          if nparam % 5 == 0:
              table=open(tableFile,'a')
              table.write('\n \\noalign{\\vskip 1.5mm} \n')
              table.close()
          #writeLine(paramList,tableFile,parameterNames[par],fittedOrDerived[par],parLabel=par)
          writeLine(paramList,tableFile,parameterNames[par],keepingTrackFitDerived,parLabel=par)
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






