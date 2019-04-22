import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob, os, pynbody# matplotlib
import pickle
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import copy
from scipy.optimize import curve_fit
fileLocationBH = '/scratch/mb6605/gasol_32_180423/'

#z0NIHAO_BHproperties.pkl
#sfrArray, mStarArray = pickle.load(open("pickles/sfrAndMstarArray.pkl", "rb"))
#galDictEllBH = pickle.load(open("pickles/AGN_ELL_NIHAOproperties.pkl", "rb"))
galDictEllBH = pickle.load(open("pickles/jsonAGN_ELL_NIHAOproperties.pkl", "rb"))
galDictEll = pickle.load(open("pickles/ELL_WOBH_NIHAOproperties.pkl", "rb"))
#NIHAO non ellipical
galDict = pickle.load(open("pickles/ALLjanNEW_NIHAOproperties.pkl", "rb"))
#galDictBH = pickle.load(open("pickles/NIHAO_BHproperties.pkl", "rb"))
galDictBH = pickle.load(open("pickles/jsonNIHAO_BHproperties.pkl", "rb"))
print(galDictEll.keys())
print(galDictEllBH.keys())
print(galDict.keys())
print(galDictBH.keys())

sfr = []
mstar = []
sfrEll = []
mstarEll = []
galsTot = 0
labeled = False
fig = plt.figure(figsize=(6,6))
for g in galDictBH.keys():
    try:
        propDict = galDict[g]
        propDictBH = galDictBH[g]
    except:
        continue
    try:
        ms = propDict['ms'][-1]
        sfr = propDict['sfr'][-1]
        z = propDict['z'][-1]
        msBH = propDictBH['ms'][-1]
        sfrBH = propDictBH['sfr'][-1]
        msBHlast = propDictBH['ms'][-1]
        sfrBHlast = propDictBH['sfr'][-1]
    except:
        continue
    print("MS")
    print(np.log10(propDictBH['ms']))
    print("SFR")
    print(np.log10(propDictBH['sfr']))
    print("Last MS", np.log10(msBHlast))
    print("Last SFR", np.log10(sfrBHlast))

    if sfr < 10**-6:
        sfr = 10**-6
    if sfrBH < 10**-6:
        sfrBH = 10**-6

    zBH = propDictBH['z'][-1]
    if zBH > 0.03:
       continue
    print(g)
    print(z)
    print(zBH)
    plt.plot(np.log10([ms, msBH]), np.log10([sfr, sfrBH]), color='k',  alpha=0.2)
    if labeled == False:
        plt.scatter(np.log10(ms), np.log10(sfr),  color ='r', marker = 'o', alpha=0.5, label = 'Without AGN')
        plt.scatter(np.log10(msBH), np.log10(sfrBH),  color ='r', alpha=0.5, marker="^", label = 'With AGN')
        labeled = True
    else:
        plt.plot(np.log10(ms), np.log10(sfr),  color ='r', marker = 'o', alpha=0.5,)
        plt.plot(np.log10(msBH), np.log10(sfrBH),  color ='r', alpha=0.5, marker="^", )
    galsTot += 1

for g in galDictEllBH.keys():
    try:
        propDict = galDictEll[g]
        propDictBH = galDictEllBH[g]
    except:
       continue
    ms = propDict['ms'][-1]
    sfr = propDict['sfr'][-1]
    z = propDict['z'][-1]
    zBH = propDictBH['z'][-1]
    msBH = propDictBH['ms'][-1]
    sfrBH = propDictBH['sfr'][-1]
    msBHlast = propDictBH['ms'][-1]
    sfrBHlast = propDictBH['sfr'][-1]
    print(propDictBH['ms'])
    print('gal', g)
    print('zBH', zBH)
    if sfr < 10**-6:
        sfr = 10**-6
    if sfrBH < 10**-6:
        sfrBH = 10**-6
    plt.plot(np.log10([ms, msBH]), np.log10([sfr, sfrBH]), color='k', alpha=0.2)
    plt.plot(np.log10(ms), np.log10(sfr),  color ='r', marker='o',  alpha=0.5, )
    plt.plot(np.log10(msBH), np.log10(sfrBH),  color ='r', alpha=0.5, marker="^")


    #plt.plot(([ms, msBH]), ([sfr, sfrBH]), color='k', alpha=0.2)
    #plt.plot((ms), (sfr),  color ='r', marker='o',  alpha=0.5, )
    #plt.plot((msBH), (sfrBH),  color ='r', alpha=0.5, marker="^")
    galsTot += 1
print(galsTot)
'''
for g in galDictEllBH.keys():
    propDict = galDictEll[g]
    z = propDict['z']
    if len(z) == 0:
        continue
    sfrEll.append(propDict['sfr'][-1])
    mstarEll.append(propDict['ms'][-1])

sfrAll = sfr + sfrEll
mstarAll = mstar + mstarEll

def func(x, SFR0, slope):
    return SFR0*((x/10**10)**slope)

popt, pcov = curve_fit(func, mstarAll, sfrAll)
SFR0 = popt[0]
slope = popt[1]
'''

x = [10.62, 10.88, 10.375, 10.1, 10.62]
y = [-0.83, -1.13, -1.6, -1.3, -0.83]
plt.plot(x,y, '-', color='k', linewidth=2, alpha =0.5)
#draw high density triangle for blue cloud
x = [ 10.05, 10.4, 8.7, 8.6, 8.1, 10.05] #3rd to last 10**8.25
y = [ 0.5, 0.1, -1.5, -1.4, -.95, 0.5] # 10**-1.25
plt.plot(x,y, '-', color='k', linewidth=2, alpha = 0.5)
#these are redshift 0.02 - 0.085



'''
plt.scatter(mstar, sfr, color='r', alpha=0.5, label="NIHAO")
plt.scatter(mstarEll, sfrEll, color ='r', alpha=0.5, marker="^", label="NIHAO Ellipticals")
fit_vals = np.asarray([10.0**6.0, 10.0**12])
fit = SFR0*((fit_vals/10**10)**slope)
plt.plot(fit_vals, fit, 'r-', label="fit")
plt.legend(loc="upper left")
'''
#ax.set_facecolor('w')
plt.ylabel("log SFR [M${_\odot}$/ year]", fontsize=18)
plt.xlabel("log M$_{\star}$ [M${_\odot}$]", fontsize=18)
plt.legend(loc="upper left")
#plt.xlim(10**8, 10**12.5)
#plt.ylim(10**-2.5, 10**1.8)
plt.xlim(8, 12.5)
plt.ylim(-2.5, 1.8)
#plt.xscale('log')
#plt.yscale('log')
plt.tight_layout()
plt.savefig('plots/2_mssfrzWandWOBH.png', dpi=300)
plt.show()
