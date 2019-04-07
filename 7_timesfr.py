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



galDictEllBH = pickle.load(open("pickles/AGN_ELL_NIHAOproperties.pkl", "rb"))
galDictEll = pickle.load(open("pickles/ELL_WOBH_NIHAOproperties.pkl", "rb"))
#NIHAO non ellipical
galDict = pickle.load(open("pickles/ALLjanNEW_NIHAOproperties.pkl", "rb"))
galDictBH = pickle.load(open("pickles/NIHAO_BHproperties.pkl", "rb"))
sfr = []
mstar = []
sfrEll = []
mstarEll = []

for g in galDict.keys():
    propDict = galDict[g]
    z = propDict['z']
    sfr.append(propDict['sfr'][-1])
    mstar.append(propDict['ms'][-1])


def plot_ga




plt.plot(x,y, '-', color='r', linewidth=2)




plt.legend(loc="upper left")
plt.ylabel("SFR [M${_\odot}$/ Year]", fontsize=18)
plt.xlabel("Time]", fontsize=18)
plt.yscale('log')
plt.savefig('plots/7_timesfr.png', dpi=300)
plt.show()
