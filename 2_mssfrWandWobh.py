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


#sfrArray, mStarArray = pickle.load(open("pickles/sfrAndMstarArray.pkl", "rb"))
galDictEllBH = pickle.load(open("pickles/AGN_ELL_NIHAOproperties.pkl", "rb"))
galDictEll = pickle.load(open("pickles/ELL_WOBH_NIHAOproperties.pkl", "rb"))
#NIHAO non ellipical
galDict = pickle.load(open("pickles/ALLjanNEW_NIHAOproperties.pkl", "rb"))
galDictBH = pickle.load(open("pickles/NIHAO_BHproperties.pkl", "rb"))



sfr = []
mstar = []
sfrEll = []
mstarEll = []
galsTot = 0
for g in galDictBH.keys():
    try:
        propDict = galDict[g]
        propDictBH = galDictBH[g]
    except:
        continue
    ms = propDict['ms'][-1]
    sfr = propDict['sfr'][-1]
    z = propDict['z'][-1]
    msBH = propDictBH['ms'][-1]
    sfrBH = propDictBH['sfr'][-1]
    zBH = propDictBH['z'][-1]
    print(g)
    print(z)
    print(zBH)
    plt.plot([ms, msBH], [sfr, sfrBH], color='k',  alpha=0.2)
    plt.plot(ms, sfr,  color ='r', marker = 'o', alpha=0.5, )
    plt.plot(msBH, sfrBH,  color ='r', alpha=0.5, marker="^")
    galsTot += 1
for g in galDictEllBH.keys():
    try:
        propDict = galDictEll[g]
        propDictBH = galDictEllBH[g]
    except:
       continue
    ms = propDict['ms'][-1]
    sfr = propDict['sfr'][-1]
    msBH = propDictBH['ms'][-1]
    sfrBH = propDictBH['sfr'][-1]
    plt.plot([ms, msBH], [sfr, sfrBH], color='k', alpha=0.2)
    plt.plot(ms, sfr,  color ='r', marker='o',  alpha=0.5, )
    plt.plot(msBH, sfrBH,  color ='r', alpha=0.5, marker="^")
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



#draw high density rectangle for red sequence
x = [10**10.6, 10**10.7, 10**10.4, 10**10.3, 10**10.6]
y = [10**-1.025, 10**-1.125, 10**-1.45, 10**-1.35, 10**-1.025]
plt.plot(x,y, '-', color='r', linewidth=2)
#draw high density triangle for blue cloud
x = [ 10**10.15, 10**10.4, 10**8.7, 10**8.6, 10**8.1, 10**10.15] #3rd to last 10**8.25
y = [ 10**0.3, 10**0.1, 10**-1.5, 10**-1.4, 10**-.95, 10**0.3] # 10**-1.25
plt.plot(x,y, '-', color='r', linewidth=2)
#these are redshift 0.02 - 0.085

#plt.scatter(mStarArray, sfrArray, color='b', alpha=0.5)

'''
plt.scatter(mstar, sfr, color='r', alpha=0.5, label="NIHAO")
plt.scatter(mstarEll, sfrEll, color ='r', alpha=0.5, marker="^", label="NIHAO Ellipticals")
fit_vals = np.asarray([10.0**6.0, 10.0**12])
fit = SFR0*((fit_vals/10**10)**slope)
plt.plot(fit_vals, fit, 'r-', label="fit")
plt.legend(loc="upper left")
'''
#ax.set_facecolor('w')
plt.ylabel("SFR [M${_\odot}$/ Year]", fontsize=18)
plt.xlabel("M$_{\star}$ [M${_\odot}$]", fontsize=18)
#plt.xlim(10**8, 10**12.5)
#plt.ylim(10**-2.5, 10**1.8)
plt.xscale('log')
plt.yscale('log')
plt.savefig('plots/2_mssfrzWandWOBH.png', dpi=300)
plt.show()
