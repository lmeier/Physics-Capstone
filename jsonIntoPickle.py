import pynbody
import sys
import numpy as np
import pylab as plt
import pynbody.analysis.profile as profile
import glob, os, pickle
import math
import json

fileLocation = '/scratch/database/nihao/gasoline2.1/'
fileLocation = '/scratch/database/nihao/gasoline2.1/'
#fileLocationAGN = '/scratch/mb6605/NIHAO_BH/'
fileLocationAGN = '/scratch/mb6605/NIHAO_ELL_BH/'

datafiles = sorted(glob.glob(fileLocationAGN + 'g?.??e??/g?.??e??' +  '.0001?'))
galaxies = []
for i in datafiles:
    galaxies.append(i[-14:-6])




galDict = dict()
for gal in galaxies:
    redshift = []
    ms = []
    sfr = []
    time = []

    datafiles = sorted(glob.glob(fileLocationAGN +'/' +gal+'/' + gal + '.0????'))
    jsondatafiles = sorted(glob.glob(fileLocationAGN + '/' + gal + '/0????/allplotdata0????.json'))
    print(gal)
    print(len(datafiles))
    print(len(jsondatafiles))
    print(datafiles)
    print(jsondatafiles)
    for df, jdf in zip(datafiles, jsondatafiles):

        f2 = open(jdf, 'r')
        data_input = json.load(f2)
        try:
            ms_cur = data_input[0]['mstar']
            sfr_cur = data_input[0]['sfr']
        except:
            continue
        sim = pynbody.load(df)
        redshift.append(sim.properties['z'])
        time.append(sim.properties['time'].in_units('Gyr'))
        ms.append(ms_cur)
        sfr.append(sfr_cur)
    print(sfr)
    print(ms)
    propDict = {'sfr': sfr, 'ms': ms, 'z': redshift, 'time': time}
    galDict[gal] = propDict

with open('pickles/jsonAGN_ELL_NIHAOproperties.pkl', 'wb') as fp:
    pickle.dump(galDict, fp)


exit()


