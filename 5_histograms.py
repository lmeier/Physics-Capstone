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
def plot_hist(deltas, concs, minConc = 8, maxConc =11, binCount=23, plotTitle='NIHAO'):
    all_delta_percentages = np.array(deltas)
    all_delta_percentages = all_delta_percentages[np.logical_not(np.isnan(all_delta_percentages))]
    n, bins, patches = plt.hist(all_delta_percentages, bins=binCount) #15 is good, as is 23
    A = np.vstack((np.digitize(all_delta_percentages, bins), concs)).T
    res_lst = np.array([np.mean(A[A[:, 0] == i, 1]) for i in range(len(bins))])
    res_lst[np.isnan(res_lst)] = min(concentrations)
    res_lst = np.array(res_lst)
    cmap = plt.get_cmap('plasma')
    sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(minConc, maxConc))
    sm._A = []

    for i in range(len(patches)):
       patches[i].set_facecolor(cmap(norm(res_lst[i],min(res_lst), max(res_lst) )))

    plt.colorbar(sm)
    plt.ylabel("Number of galaxy datapoints", fontsize=18)
    plt.xlabel("Delta Log(SFR) [M$_\odot$/ Year]", fontsize=18)
    plt.xlim(-3, 3)
    plt.yscale('log')
    plt.title(plotTitle)
    #plt.savefig(str("del_sfr_ELL_plots/NIHAO_hist.png"), dpi=300)
    plt.show()

def make_delta_conc_lists(gals, galDict, concDict, avgFitDict):
    all_delta_percentages = []
    all_delta_after_log = []
    sfrEndAtTime = []
    sfrEndConc = []
    allConc = []
    allMstar = []
    for num, g in enumerate(gals):
        propDict = galDict[g]
        conc0 = concDict[g]
        if num == 0:
            cur_min = conc0
            cur_max = conc0
        else:
            if conc0 < cur_min:
                cur_min = conc0
            if conc0 > cur_max:
                cur_max = conc0
    for g in gals:
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

alDictEllBH = pickle.load(open("pickles/AGN_ELL_NIHAOproperties.pkl", "rb"))
galDictEll = pickle.load(open("pickles/ELL_WOBH_NIHAOproperties.pkl", "rb"))
#NIHAO non ellipical
galDict = pickle.load(open("pickles/ALLjanNEW_NIHAOproperties.pkl", "rb"))
galDictBH = pickle.load(open("pickles/NIHAO_BHproperties.pkl", "rb"))


current_galDict = pickle.load(open('/home/lem507/2018/pickles/ALLjanNEW_NIHAOproperties.pkl', 'rb'))
current_concDict = pickle.load(open('/home/lem507/2018/pickles/concentrations.pkl', 'rb'))
concentrations = current_concDict.values()
current_gals = [g for g in current_galDict.keys()]


current_avgFitDict = make_fits(current_galDict)
delta_list, conc_list = make_delta_conc_lists(current_gals, current_galDict, current_concDict, current_avgFitDict)
plot_hist(delta_list, conc_list)
#def plot_hist(deltas=all_delta_after_log, concs=allConc, minConc = 8, maxConc =11, binCount=23, plotTitle='NIHAO'):
