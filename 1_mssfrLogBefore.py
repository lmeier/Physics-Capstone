import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob, os, pynbody# matplotlib
import pickle
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import copy
from scipy.optimize import curve_fit
from scipy import stats
fileLocationBH = '/scratch/mb6605/gasol_32_180423/'

import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 1.3
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['legend.scatterpoints'] = 1
mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.handlelength'] = 1.4
mpl.rcParams['legend.handletextpad'] = 0.5
#mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'cmr10'
mpl.rcParams['mathtext.fontset'] = 'cm'
#mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.unicode_minus']=False

#sfrArray, mStarArray = pickle.load(open("pickles/sfrAndMstarArray.pkl", "rb"))
galDict = pickle.load(open("pickles/ALLjanNEW_NIHAOproperties.pkl", "rb"))
galDictEll = pickle.load(open("pickles/jsonAGN_ELL_NIHAOproperties.pkl", "rb"))
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

logsfr = []
logmstar = []
for i, j in zip(sfr, mstar):
    if i > 0:
        logsfr.append(np.log10(i))
        logmstar.append(np.log10(j))
slope, intercept, r_value, p_value, std_err = stats.linregress(logmstar, logsfr)
print(slope)
print(intercept)



'''
def func(x, SFR0, slope):
    return SFR0*((x/10**10)**slope)

popt, pcov = curve_fit(func, mstar, sfr)
SFR0 = popt[0]
slope = popt[1]
'''

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


fig = plt.figure(figsize=(6,6))
#bottom right
# tl, tr, br, bl tl
x = [10.8, 10.98, 10.38, 10.2, 10.8]
y = [-.95, -1.20, -1.6, -1.35, -.95]
plt.plot(x,y, '-', color='k', linewidth=2, alpha =0.5)

#top left box
x = [ 10.15, 10.5, 8.8, 8.1, 10.15] #3rd to last 10**8.25
y = [ 0.25, -0.1, -1.6, -.95, 0.25] # 10**-1.25
plt.plot(x,y, '-', color='k', linewidth=2, alpha = 0.5)


#plt.scatter(mStarArray, sfrArray, color='b', alpha=0.5)
plt.scatter(logmstar, logsfr, color='r', alpha=0.5, label="NIHAO")
plt.scatter(np.log10(mstarEll), np.log10(sfrEll), color ='k', alpha=0.5, marker="^", label="NIHAO LARGE w/ BH")

fit_vals = np.asarray([6.0, 13])
#fit = SFR0*((fit_vals/10**10)**slope)
fit = slope*fit_vals + intercept
plt.plot(fit_vals, fit, 'r-', label="fit")

plt.legend(loc="upper left", numpoints=1)
#ax.set_facecolor('w')
plt.ylabel("log SFR [M${_\odot}$/ year]", fontsize=18)
plt.xlabel("log M$_{\star}$ [M${_\odot}$]", fontsize=18)
plt.xlim(8, 12.5)
plt.ylim(-2.5, 1.8)
plt.title('z = 0')
#plt.xscale('log')
#plt.yscale('log')
plt.tight_layout()
plt.savefig('plots/1_mssfrz0logbeforeNoRen.png', dpi=300)
plt.show()


fig = plt.figure(figsize=(6,6))
#bottom right
# tl, tr, br, bl tl

x = [10.8, 10.98, 10.38, 10.2, 10.8]
y = [-.95, -1.20, -1.6, -1.35, -.95]
plt.plot(x,y, '-', color='k', linewidth=2, alpha =0.5)

#top left box
x = [ 10.15, 10.5, 8.8, 8.1, 10.15] #3rd to last 10**8.25
y = [ 0.25, -0.1, -1.6, -.95, 0.25] # 10**-1.25
plt.plot(x,y, '-', color='k', linewidth=2, alpha = 0.5)





#plt.scatter(mStarArray, sfrArray, color='b', alpha=0.5)
plt.scatter(logmstar, logsfr, color='r', alpha=0.5, label="NIHAO")


#plt.legend(loc="upper left")
#ax.set_facecolor('w')
plt.ylabel("log SFR [M${_\odot}$/ year]", fontsize=18)
plt.xlabel("log M$_{\star}$ [M${_\odot}$]", fontsize=18)
plt.xlim(8, 12.5)
plt.ylim(-2.5, 1.8)
plt.title('z = 0')
#plt.xscale('log')
#plt.yscale('log')
plt.tight_layout()
plt.savefig('plots/1_mssfrz0NoElllogbeforeNoRen.png', dpi=300)
plt.show()
