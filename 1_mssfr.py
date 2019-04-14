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
galDict = pickle.load(open("pickles/ALLjanNEW_NIHAOproperties.pkl", "rb"))
galDictEll = pickle.load(open("pickles/AGN_ELL_NIHAOproperties.pkl", "rb"))
sfr = []
mstar = []
sfrEll = []
mstarEll = []

print(len(galDict.keys()))
for g in galDict.keys():
    propDict = galDict[g]
    z = propDict['z']
    sfr.append(propDict['sfr'][-1])
    mstar.append(propDict['ms'][-1])


for g in galDictEll.keys():
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
'''
#bottom right
# tl, tr, br, bl tl
x = [10**10.62, 10**10.88, 10**10.375, 10**10.1, 10**10.62]
y = [10**-0.83, 10**-1.13, 10**-1.6, 10**-1.3, 10**-0.83]
plt.plot(x,y, '-', color='r', linewidth=2, alpha =0.5)
#draw high density triangle for blue cloud
x = [ 10**10.05, 10**10.4, 10**8.7, 10**8.6, 10**8.1, 10**10.05] #3rd to last 10**8.25
y = [ 10**0.5, 10**0.1, 10**-1.5, 10**-1.4, 10**-.95, 10**0.5] # 10**-1.25
plt.plot(x,y, '-', color='r', linewidth=2, alpha = 0.5)
#these are redshift 0.02 - 0.085

#plt.scatter(mStarArray, sfrArray, color='b', alpha=0.5)

plt.scatter(mstar, sfr, color='r', alpha=0.5, label="NIHAO")
plt.scatter(mstarEll, sfrEll, color ='r', alpha=0.5, marker="^", label="NIHAO Ellipticals")

fit_vals = np.asarray([10.0**6.0, 10.0**13])
fit = SFR0*((fit_vals/10**10)**slope)
plt.plot(fit_vals, fit, 'r-', label="fit")

plt.legend(loc="upper left")
#ax.set_facecolor('w')
plt.ylabel("SFR [M${_\odot}$/ Year]", fontsize=18)
plt.xlabel("M$_{\star}$ [M${_\odot}$]", fontsize=18)
plt.xlim(10**8, 10**12.5)
plt.ylim(10**-2.5, 10**1.8)
plt.xscale('log')
plt.yscale('log')
plt.savefig('plots/1_mssfrz0.png', dpi=300)
plt.show()
