"""
Read in par files and compare old and new values
It needs a list of files to compare which are in
./comparisonFileLists
"""

import os, json
import numpy as np
from astropy import units as u
import astropy
from uncertainties import ufloat
import uncertainties
import decimal
import string

import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord


import sys
sys.path.append('../tableScripts')
# read in par file - copied from Daniel's script
import readParFile



# do the n Sigma range overlap?
def insideOneSigma(par1,err1,par2,err2,nSig=1):

    min1 = float(par1)-nSig*float(err1)
    max1 = float(par1)+nSig*float(err1)
    min2 = float(par2)-nSig*float(err2)
    max2 = float(par2)+nSig*float(err2)

    if   min1>max2: overlaps = 'False'
    elif min2>max1: overlaps = 'False'
    else: overlaps = 'True'
    #print (overlaps)

    # check if the parameter is zero
    if any([par1,par2,err1,err2])==0.0:
      overlaps='n/a'
    return overlaps




def getLabels(path):


   if ('kop' in path)==True:
       if ('kop_alt' in path)==True:
           pathLabel = 'kop_alt'
       else: 
           pathLabel = 'kop'

   elif ('basic' in path)==True:
       pathLabel = 'basic' 
   

   elif ('t2' in path)==True:
       pathLabel = 't2'

   else: pathLabel = 'default'


   return pathLabel



# not used at the moment
def plotThisParamPSR(new,old,param):

   oldVal=old[param]
   newVal=new[param]

   oldErr=old[param+'_ERR']
   newErr=new[param+'_ERR']
   print(oldVal,newVal,'\n',oldErr,newErr)

   """
   plt.scatter(1,newVal)
   plt.errorbar(1,newVal,yerr=newErr,label='new')
   plt.scatter(2,oldVal)
   plt.errorbar(2,oldVal,yerr=oldErr,label='old')
   plt.legend()
   """
   plt.scatter(oldVal,newVal)
   plt.errorbar(oldVal,newVal,xerr=oldErr,yerr=newErr)

  
   minimum = min(oldVal-oldErr,newVal-newErr)
   maximum = max(oldVal+oldErr,newVal+newErr)
   difference = maximum - minimum

   #plt.ylim(minimum-(0.1*difference),maximum+(0.1*difference))
   #plt.show()







#which data set to compare - choose from 
# ppta15
# eptaDR1
# nanograv11yr
# nanograv12.5yrNB
# nanograv12.5yrWB 	
whichDatasetToCompare = 'nanograv12.5yrWB'

filePaths = np.genfromtxt('comparisonFileLists/{}CommonFiles.list'.format(whichDatasetToCompare),dtype='str')
psrList = filePaths[:,0]
pptaDR2Paths = filePaths[:,1]
comparisonPaths = filePaths[:,2]


# get pulsar name right for different data releases
if whichDatasetToCompare == 'ppta15' or whichDatasetToCompare=='eptaDr1':
   psrLabel = 'PSRJ'
else: 
   psrLabel = 'PSR'


# what parameters do we want to compare? 
params = (['M2'])




for par in params: 


    print("""

    Comparing PPTADR2 to {}
    Results for parameter {}
    
    Do the uncertainties overlap?

    """.format(whichDatasetToCompare,par))


    print('PSR Name\tDR2Label\tComp.Label\t1sig\t2sig\t3sig')

    for psr, newParPath, oldParPath in zip(psrList,pptaDR2Paths,comparisonPaths): 

        newResult = readParFile.read_par(newParPath)
        oldResult = readParFile.read_par(oldParPath)

        newLabel,oldLabel = getLabels(newParPath), getLabels(oldParPath)

        #print('alt' in newParPath)

        try: 
            compare1Sig = insideOneSigma(newResult[par],newResult[par+'_ERR'],oldResult[par],oldResult[par+'_ERR'],nSig=1)
            compare2Sig = insideOneSigma(newResult[par],newResult[par+'_ERR'],oldResult[par],oldResult[par+'_ERR'],nSig=2)
            compare3Sig = insideOneSigma(newResult[par],newResult[par+'_ERR'],oldResult[par],oldResult[par+'_ERR'],nSig=3)
            print('{}\t{}\t\t{}\t\t{}\t{}\t{}'.format(psr,newLabel,oldLabel,compare1Sig,compare2Sig,compare3Sig))
        except:
            print('{}\t{}\t\t{}\t\tn/a\tn/a\tn/a'.format(psr,newLabel,oldLabel))








"""
old version for writing to files

          absDiff = abs(newResult[par]-oldResult[par])
          percentageDiff = (absDiff/abs(oldResult[par]))*100.       

          compare.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(oldResult[psrLabel],
                                                                                  newResult[par],
                                                                                  newResult[par+'_ERR'],
                                                                                  oldResult[par],
                                                                                  oldResult[par+'_ERR'],
                                                                                  compare1Sig,
                                                                                  compare2Sig,
                                                                                  compare3Sig,
                                                                                  absDiff,
                                                                                  percentageDiff,
                                                                                  absDiff/abs(newResult[par+'_ERR']),
                                                                                  absDiff/abs(oldResult[par+'_ERR']),
                                                                                  absDiff/(abs(newResult[par+'_ERR'])+abs(oldResult[par+'_ERR']))))

          print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(oldResult[psrLabel],compare1Sig,compare2Sig,compare3Sig,
                                                    absDiff/abs(newResult[par+'_ERR']),
                                                    absDiff/abs(oldResult[par+'_ERR']),
                                                    absDiff/(abs(newResult[par+'_ERR'])+abs(oldResult[par+'_ERR']))))

        except:
          compare.write('{}\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\n'.format(oldResult[psrLabel]))
          print('{}\t n/a\tn/a\tn/a\tn/a\tn/a\tn/a'.format(oldResult[psrLabel]))



    compare.close()

"""
