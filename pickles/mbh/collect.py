import pynbody
import numpy as np
import pylab as plt 
import pynbody.analysis.profile as profile
import glob, os, pickle
import math

'''
mbh-mstar-g1.05e11-NA.pkl  mbh-mstar-g2.19e11-NA.pkl  mbh-mstar-g8.26e11-NA.pkl
mbh-mstar-g1.12e12-NA.pkl  mbh-mstar-g2.79e12-NA.pkl
mbh-mstar-g1.44e10-NA.pkl  mbh-mstar-g7.55e11-NA.pkl
'''


#datafiles = glob.glob('/home/lem507/2018/pickles/mbh/mbh*.pkl')

datafiles = glob.glob('/home/lem507/2018/nadine//mbh*.pkl')

galDict = dict()
for d in datafiles:
    print(d)
    #/home/lem507/2018/pickles/mbh/mbh-mstar-g1.44e10.pkl
    #/home/lem507/2018/nadine/mbh-mstar-g8.26e11-redshift.pkl

    galaxy = d[-21:-13]
    print(galaxy)
    currentDict =  pickle.load(open(d, 'rb'))
    print(currentDict.keys())
    galDict[galaxy] = currentDict


#agnQuenchingProperties.pkl
with open('/home/lem507/2018/pickles/agnQuenchingBHonly.pkl', 'wb') as fp:
    pickle.dump(galDict, fp)
