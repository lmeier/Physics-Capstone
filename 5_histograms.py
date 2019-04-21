import pynbody
import numpy as np
import pylab as plt
import pynbody.analysis.profile as profile
import pynbody.filt as filt
import glob, os, pickle
import math
import pickle
from sets import Set
import os.path
import json
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import linregress
#from redshifts import checkRedshifts
fileLocation = '/scratch/database/nihao/gasoline2.1/'
fileLocationELL = '/scratch/mb6605/NIHAO_ELL_BH/'
fileLocationNIHAO = '/scratch/mb6605/NIHAO_BH/'
fileLoc = '/scratch/mb6605/'

def myExpFunc(x, a, b):
    return a * np.power(x, b)

def linearFunc(x, m, b):
    return m

def norm(x, minimum, maximum):
    return ((x- minimum) / (maximum - minimum))

def linear_to_log_func(x, m, b):
    intermediate = m*x + b
    return np.power(10, intermediate)

def make_fits(galaxies):
    timeMsDict = dict()
    timeSfrDict = dict()
    for g in galaxies.keys():
        propDict = galaxies[g]
        ms = propDict['ms']
        redshift = propDict['z']
        timeNew = propDict['time']
        sfr = propDict['sfr']
        print(sfr)
        for idx, t in enumerate(timeNew):
            try:
                ms_list = timeMsDict[t]
                sfr_list = timeSfrDict[t]
                ms_list.append(ms[idx])
                sfr_list.append(sfr[idx])
                timeMsDict[t] = ms_list
                timeSfrDict[t] = sfr_list
            except:
                ms_list = []
                ms_list.append(ms[idx])
                sfr_list = []
                sfr_list.append(sfr[idx])
                timeMsDict[t] = ms_list
                timeSfrDict[t] = sfr_list
    avg_fit = dict()
    for t, ms in timeMsDict.items():
        print('time')
        sfr = timeSfrDict[t]
        print(t)
        print('length of ms')
        print(len(ms))
        print(len(sfr))
        if len(sfr) < 15 or len(ms) < 15:
            print('continuing')
            continue
        try:
            sfr_filt, ms_filt = zip(*((id, other) for id, other in zip(sfr, ms) if (id  >  10**-4)))
        except:
            print(sfr)
            print(ms)
            continue
        if len(ms_filt) < 15:
            continue
        log_ms = np.log10(ms_filt)
        log_sfr = np.log10(sfr_filt)
        m, b = linregress(log_ms, log_sfr)[0:2]
        avg_fit[t] = (m, b)
    return avg_fit

#def plot_hist(deltas=all_delta_after_log, concs=allConc, minConc = 8, maxConc =11, binCount=23, plotTitle='NIHAO'):
mins = []
maxs = []
def plot_hist(deltas, concs, pointToRemove=[-1], minConc = 8, maxConc =11, binCount=23, plotTitle='NIHAO'):
    all_delta_percentages = np.array(deltas)
    all_delta_percentages = all_delta_percentages[np.logical_not(np.isnan(all_delta_percentages))]
    n, bins, patches = plt.hist(all_delta_percentages, bins=binCount) #15 is good, as is 23
    A = np.vstack((np.digitize(all_delta_percentages, bins), concs)).T
    res_lst = np.array([np.mean(A[A[:, 0] == i, 1]) for i in range(len(bins))])
    res_lst[np.isnan(res_lst)] = min(concentrations)
    res_lst = np.array(res_lst)
    maxC = np.average(res_lst)
    minC = np.average(res_lst)
    for i, j in zip(n, res_lst):
        if i < 7:
            continue
        if j < minC:
            minC = j
        if j > maxC:
            maxC = j
    print(maxC)
    maxs.append(maxC)
    print(minC)
    mins.append(minC)
    print("n", n)
    print('bins', bins)
    print("A", A)
    print("res_lst", res_lst)
    minC = 6.5
    maxC = 12.1
    cmap = plt.get_cmap('plasma')
    #sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(minConc, maxConc))
    sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(minC, maxC))
    #sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(min(res_lst), max(res_lst)))
    sm._A = []

    for i in range(len(patches)):
       if i in pointToRemove:
           patches[i].set_edgecolor('white')
           patches[i].set_facecolor('white')
       else:
           #patches[i].set_facecolor(cmap(norm(res_lst[i],minConc, maxConc) ))
           patches[i].set_facecolor(cmap(norm(res_lst[i],minC, maxC) ))
           #patches[i].set_facecolor(cmap(norm(res_lst[i],min(res_lst), max(res_lst) )))
    cbar = plt.colorbar(sm,)
    cbar.set_label('c', rotation=270, fontsize=18)
    plt.ylabel("Number of galaxy datapoints", fontsize=18)
    plt.xlabel("Delta Log(SFR) [M$_\odot$/ Year]", fontsize=18)
    plt.xlim(-3, 3)
    plt.yscale('log')
    plt.title(plotTitle)
    plt.savefig(str("plots/9_hist_" + plotTitle.replace(" ", "") + ".png"), dpi=300)
    plt.show()

def make_delta_conc_lists(gals, galDict, concDict, avgFitDict):
    all_delta_percentages = []
    all_delta_after_log = []
    sfrEndAtTime = []
    sfrEndConc = []
    allConc = []
    allMstar = []
    minMaxSetYet = False
    for num, g in enumerate(gals):
        if g not in concDict.keys():
            galsNotIncluded.append(g)
            continue
        propDict = galDict[g]
        conc0 = concDict[g]
        if minMaxSetYet == False:
            cur_min = conc0
            cur_max = conc0
            minMaxSetYet = True
        else:
            if conc0 < cur_min:
                cur_min = conc0
            if conc0 > cur_max:
                cur_max = conc0
    for g in gals:
        if g not in concDict.keys():
            continue
        propDict = galDict[g]
        conc0 = concDict[g]
        #concentration0 = "n/a"
        ms = propDict['ms']
        redshift = propDict['z']
        timeNew = propDict['time']
        sfr = propDict['sfr']
        delta_sfr = []
        delta_sfr_percentage = []
        delta_after_log = []
        concList = []
        msList = []
        time = timeNew[:]
        timeNew = []
        time.reverse()
        ms.reverse()
        sfr.reverse()
        alreadyEnded = False
        try:
            prev_m = ms[0]
            prev_t = time[0]
            prev_sfr = sfr[0]
        except:
            pass
        for ms_a, s, t in zip(ms, sfr, time):
            try:
                m_cur, b_cur = avgFitDict[t]
                #trying adding 1 to break the divide by 0 error
                #undid
                avg = linear_to_log_func(np.log10(ms_a), m_cur, b_cur)
                #if alreadyEnded == True:
                #    continue
                if s <= 0:
                   sfrEndAtTime.append(t)
                   sfrEndConc.append(conc0)
                   alreadyEnded = True
                   continue
                delta = s - avg
                delta_sfr.append(delta)
                delta_sfr_percentage.append(delta/avg)
                delta_after_log.append( np.log10(s) - np.log10(avg))
                concList.append(conc0)
                timeNew.append(t)
                msList.append(ms_a)
            except:
                continue
        timeNew.reverse()
        delta_sfr.reverse()
        delta_sfr_percentage.reverse()
        delta_after_log.reverse()
        all_delta_percentages.extend(delta_sfr)
        all_delta_after_log.extend(delta_after_log)
        allConc.extend(concList)
        allMstar.extend(msList)

    return [all_delta_after_log, allConc]

#def make_hist(galDict, galDictKeys, concDict, title):


#computing average fit

names =["Total_NIHAO_Classic", "NIHAO_Classic_High_Mass",
"NIHAO_Classic_Low_Mass", "NIHAO_Ellipticals_All", "NIHAO_With_AGN", "NIHAO_Without_AGN"]

galDictEllBH = pickle.load(open("pickles/AGN_ELL_NIHAOproperties.pkl", "rb"))
galDictEll = pickle.load(open("pickles/ELL_WOBH_NIHAOproperties.pkl", "rb"))
#NIHAO non ellipical
galDict = pickle.load(open("pickles/ALLjanNEW_NIHAOproperties.pkl", "rb"))
galDictBH = pickle.load(open("pickles/NIHAO_BHproperties.pkl", "rb"))

#DMO2 because thats the right conc's
concDictNIHAO = pickle.load(open('/home/lem507/2018/pickles/concentrationsDMO2.pkl', 'rb'))
concDictEll = pickle.load(open('/home/lem507/2018/pickles/concentrationsDMO_AGN_ELL2.pkl', 'rb'))
allConc = concDictNIHAO.copy()
allConc.update(concDictEll)
print(allConc.keys())
galsNotIncluded = []

current_galDict = pickle.load(open('/home/lem507/2018/pickles/ALLjanNEW_NIHAOproperties.pkl', 'rb'))
concentrations = allConc.values()
current_gals = [g for g in current_galDict.keys()]

highMassGals = []
lowMassGals = []
for g in galDict.keys():
    if int(g[-2:]) > 11 or (int(g[-2:]) == 11 and int(g[1]) >= 5):
       highMassGals.append(g)
    else:
        lowMassGals.append(g)
#NIHAO Classic
avgFitDict = make_fits(galDict)
delta_list, conc_list = make_delta_conc_lists(galDict.keys(), galDict, allConc, avgFitDict)
plot_hist(delta_list, conc_list,pointToRemove=[2], minConc=5, maxConc=11.5,  plotTitle='NIHAO Classic')

#NIHAO Classis High Mass
avgFitDict = make_fits(galDict)
delta_list, conc_list = make_delta_conc_lists(highMassGals, galDict, allConc, avgFitDict)
plot_hist(delta_list, conc_list, minConc=6, maxConc=9.6, plotTitle = 'NIHAO Classic High Mass')

#NIHAO Classis Low Mass
avgFitDict = make_fits(galDict)
delta_list, conc_list = make_delta_conc_lists(lowMassGals, galDict, allConc, avgFitDict)
plot_hist(delta_list, conc_list, minConc=6, maxConc=15, plotTitle='NIHAO Classic Low Mass')

#NIHAO with Ellipticals
nihaoWithEll = galDict.copy()
nihaoWithEll.update(galDictEllBH)
avgFitDict = make_fits(nihaoWithEll)
delta_list, conc_list = make_delta_conc_lists(nihaoWithEll.keys(), nihaoWithEll, allConc, avgFitDict)
plot_hist(delta_list, conc_list, pointToRemove = [0], minConc = 6, maxConc = 12, plotTitle="NIHAO with Ellipticals")

#With AGN
withAgn = galDictBH.copy()
withAgn.update(galDictEllBH)
avgFitDict = make_fits(withAgn)
delta_list, conc_list = make_delta_conc_lists(withAgn.keys(), withAgn, allConc, avgFitDict)
plot_hist(delta_list, conc_list, pointToRemove=[1], minConc=5, maxConc=9.6,  plotTitle="Galaxies with AGN")

#Without AGN
withoutAgn = galDict.copy()
withAgn.update(galDictEll)
avgFitDict = make_fits(withoutAgn)
delta_list, conc_list = make_delta_conc_lists(withoutAgn.keys(), withoutAgn, allConc, avgFitDict)
plot_hist(delta_list, conc_list, minConc = 5, maxConc = 15, plotTitle="Galaxies without AGN")

print(galsNotIncluded)


print(mins)
print(maxs)




#def plot_hist(deltas=all_delta_after_log, concs=allConc, minConc = 8, maxConc =11, binCount=23, plotTitle='NIHAO'):
