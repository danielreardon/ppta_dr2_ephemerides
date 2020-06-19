"""
This script finds which pulsars are common to two datasets
and makes a list of file paths to par files which we'll 
compare in ../psrComparison.py 

All the file paths point to places on OzSTAR in oz002
"""
import numpy as np
import glob


# pick one of:
#   ppta15
#   eptaDR1
#   nanograv11yr
#   nanograv12.5yrNB
#   nanograv12.5yrWB
whichDataset = 'nanograv12.5yrWB'

if whichDataset   == 'ppta15':
    pathToComparisonDataset = '/fred/oz002/hmiddleton/ppta_ephemeris/dataSets/{}'.format(whichDataset)
elif whichDataset == 'eptaDR1':
    pathToComparisonDataset = '/fred/oz002/hmiddleton/ppta_ephemeris/dataSets/eptaDR1/EPTA_v2.2_git'
elif whichDataset == 'nanograv11yr':
    pathToComparisonDataset = '/fred/oz002/hmiddleton/ppta_ephemeris/dataSets/NANOGrav_11y/par/'
elif whichDataset == 'nanograv12.5yrNB':
    pathToComparisonDataset = '/fred/oz002/hmiddleton/ppta_ephemeris/dataSets/NANOGrav_12yv3/narrowband/par/'
elif whichDataset == 'nanograv12.5yrWB':
    pathToComparisonDataset = '/fred/oz002/hmiddleton/ppta_ephemeris/dataSets/NANOGrav_12yv3/wideband/par/'
else: 
    print('{} dataset not found!'.format(whichDataset))
    exit()

# the list of all PPTA DR2 pulsars
allPSRList = np.genfromtxt('allPPTADR2PSRs.list',dtype=str)

# a place to save the files for common pulsars (one for each comparison dataset)
commonFiles = open('{}CommonFiles.list'.format(whichDataset),'w')



for psr in allPSRList: 

    # this is a list of files for the PPTA DR2 dataset for each pulsar
    pptaDR2FileNames = glob.glob('/fred/oz002/hmiddleton/ppta_ephemeris/repositories/ppta_dr2_ephemerides/publish_collection/dr2/{}*.par'.format(psr))

    # does this pulsar exist in the other dataset? 
    if whichDataset == 'ppta15':
        comparisonFileNames = glob.glob('{}/{}*.par'.format(pathToComparisonDataset,psr))
    elif whichDataset == 'eptaDR1':
        comparisonFileNames = glob.glob('{0}/{1}/{1}*.par'.format(pathToComparisonDataset,psr))
    elif whichDataset == 'nanograv11yr':
        comparisonFileNames = glob.glob('{}/{}*.par'.format(pathToComparisonDataset,psr))
    elif whichDataset == 'nanograv12.5yrNB':
        comparisonFileNames = glob.glob('{}/{}*.par'.format(pathToComparisonDataset,psr))
    elif whichDataset == 'nanograv12.5yrWB':
        comparisonFileNames = glob.glob('{}/{}*.par'.format(pathToComparisonDataset,psr))
    else: 
        print(' ')
        exit()


    # skip if the pular is not in the other dataset 
    if len(comparisonFileNames)==0:
        print('skipping {}'.format(psr))
        pass

    # if the pulsar common to both datasets, save the paths compare to a file
    else:
        for new in pptaDR2FileNames: 
            for old in comparisonFileNames:
                commonFiles.write('{}\t{}\t{}\n'.format(psr,new,old))

commonFiles.close()

