import pynbody
import numpy as np
import pylab as plt
import matplotlib.patheffects as path_effects
import pynbody.analysis.profile as profile
import pynbody.filt as filt
import glob, os, pickle
import math
import pickle
import os.path
import json
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import linregress
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


fileLocation = '/scratch/database/nihao/gasoline2.1/'
fileLocationELL = '/scratch/mb6605/NIHAO_ELL_BH/'
fileLocationNIHAO = '/scratch/mb6605/NIHAO_BH/'
fileLoc = '/scratch/mb6605/'
fileBase = '/Users/liammeier/FalconLocal'
fileBaseFalcon = '/home/lem507'

sfr_all = []
time_all = []
z_all = []
mvir_all = []
mstar_all = []
haloid_all = []

def smooth_for_plotting(x,y):
    new_a = []
    new_b = []
    for num, (a,b) in enumerate(zip(x,y)):
        if num == 0:
            prev_a = a
            prev_b = b
        else:
           new_a.append((a+prev_a)/2)
           new_b.append((b+prev_b)/2)
           prev_a = a
           prev_b = b
    return new_a, new_b

def myExpFunc(x, a, b):
    return a * np.power(x, b)

def linearFunc(x, m, b):
    return m

def linear_to_log_func(x, m, b):
    intermediate = m*x + b
    return np.power(10, intermediate)

def norm_one(x):
    return ((x - min(concentrations)) / (max(concentrations) - min(concentrations)))

def norm(x, minimum, maximum):
    return ((x- minimum) / (maximum - minimum))

mstar_nihao = []
sfr_nihao = []
time_nihao = []

#galDict = pickle.load(open(fileBase + '/2018/pickles/AGN_ELL_NIHAOproperties.pkl', 'rb'))

galDict = pickle.load(open('/home/lem507/2018/pickles/ALLjanNEW_NIHAOproperties.pkl', 'rb'))
galaxies = galDict.keys()
galaxies.sort()

timeMsDict = dict()
timeSfrDict = dict()
nihaoGalDict = pickle.load(open('/home/lem507/2018/pickles/ALLjanNEW_NIHAOproperties.pkl', 'rb'))

for g in galDict.keys():
    propDict = galDict[g]
    ms = propDict['ms']
    redshift = propDict['z']
    timeNew = propDict['time']
    sfr = propDict['sfr']
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


avg_fit_dict = dict()

for t, ms in timeMsDict.items():
    sfr = timeSfrDict[t]
    if len(sfr) < 15 or len(ms) < 15:
        print('continuing')
        continue
    try:
        sfr_filt, ms_filt = zip(*((id, other) for id, other in zip(sfr, ms) if (id  >  10**-4)))
    except:
        continue
    if len(ms_filt) < 15:
        continue
    log_ms = np.log10(ms_filt)
    log_sfr = np.log10(sfr_filt)
    m, b = linregress(log_ms, log_sfr)[0:2]
    avg_fit_dict[t] = (m, b)
    fit_vals = np.asarray([10.0**4.0, 10.0**6, 10.0**10, 10.0**14])


plt.plot([0.0,14], [0,0], color='k')
n = 1 #number of divisions

all_delta_percentages = []
all_delta_after_log = []

sfrEndAtTime = []
sfrEndConc = []

concDict = pickle.load(open(fileBaseFalcon + '/2018/pickles/concentrations.pkl', 'rb'))
concentrations = concDict.values()
allConc = []
allMstar = []


cmap = plt.get_cmap('plasma')
#colors = [cmap(i) for i in np.linspace(0, 1, number)]

#hacky code for getting a colorbar
sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(vmin=min(concentrations), vmax=max(concentrations)))
# fake up the array of the scalar mappable. Urgh...
sm._A = []
#plt.colorbar(sm)


#for a in ['g3.59e12', 'g1.05e13', 'g7.92e12']:
#    galDict.pop(a)

'''
g1 = [g for g in galDict.keys()[0:10]]
g2 = [g for g in galDict.keys()[10:20]]
g3 = [g for g in galDict.keys()[20:30]]
g4 = [g for g in galDict.keys() if g not in g1+g2+g3]
'''

g1 = ['g8.26e11', 'g5.38e11']
g2 = ['g8.26e11', 'g3.21e11']
g3 = ['g8.26e11', 'g4.90e11']
g4 = ['g8.26e11', 'g5.55e11']
g5 = ['g8.26e11']
#gal_all = ['g1.95e10', 'g4.86e10', 'g6.77e10', 'g9.59e10', 'g1.05e11', 'g1.37e11', 'g1.52e11',
# 'g2.04e11', 'g3.06e11', 'g3.21e11', 'g4.90e11', 'g5.38e11', 'g5.55e11', 'g6.96e11',  'g7.55e11',
# 'g8.26e11', 'g8.28e11', 'g1.05e11', 'g1.12e12', 'g2.79e12']
gal_arr = [g5]
#gal_arr =[gal_all]
size_names = [ "2gals",]
#size_names = ["1"]
plot = [True, True, True, True]

for i, name, p in zip(gal_arr, size_names, plot):
    for num, g in enumerate(i):
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
    for g in i:
        propDict = galDict[g]
        conc0 = concDict[g]
        ms = propDict['ms']
        redshift = propDict['z']
        timeNew = propDict['time']
        sfr = propDict['sfr']
        mstar_nihao.extend(ms)
        sfr_nihao.extend(sfr)
        delta_sfr = []
        delta_sfr_percentage = []
        delta_after_log = []
        concList = []
        msList = []
        time = timeNew[:]
        timeNew = []

        #reversing these three to help put the 10**-1 cutoff
        time.reverse()
        ms.reverse()
        sfr.reverse()

        above_cutoff_yet = False
        alreadyEnded = False
        try:
            prev_m = ms[0]
            prev_t = time[0]
            prev_sfr = sfr[0]
        except:
            pass
        for ms_a, s, t in zip(ms, sfr, time):
            try:
                m_cur, b_cur = avg_fit_dict[t]
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
        #reverse everything to correct time
        timeNew.reverse()
        delta_sfr.reverse()
        delta_sfr_percentage.reverse()
        delta_after_log.reverse()
        all_delta_percentages.extend(delta_sfr)
        all_delta_after_log.extend(delta_after_log)
        allConc.extend(concList)
        allMstar.extend(msList)
        if p == False:
            continue

        a, b = smooth_for_plotting(timeNew, delta_after_log)
        a, b = smooth_for_plotting(a,b)
        #label_long = g + ' ' +  str(conc0)[0:3]
        label_long = g
        linestyle='-'
        #color=cmap(norm(conc0, cur_min, cur_max))
        #plt.plot(a, b, color=cmap(norm(conc0, cur_min, cur_max)),  label=label_long,  ls=linestyle,path_effects=[path_effects.SimpleLineShadow(),path_effects.Normal()] )
        plt.plot(a, b, linewidth=2,  label=label_long,  ls=linestyle,path_effects=[path_effects.Stroke(linewidth=3, foreground='k'), path_effects.Normal()])
        #plt.plot(a, b, color='k',  label=label_long,  ls=linestyle, linewidth)
        #plt.plot(a, b, color=cmap(norm(conc0, cur_min, cur_max)),  label=label_long,  ls=linestyle,
        #plt.plot(a, b, color=cmap(norm(conc0, cur_min, cur_max)),linewidth =3,  label=label_long,  ls=linestyle,)

    if p == False:
        break

    sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(vmin=cur_min, vmax=cur_max))
    #make empty array to do custom color map
    sm._A = []
    #cbar = plt.colorbar(sm,)
    #cbar.set_label('c', rotation=270, fontsize=18)

    #plt.plot([0.0, 12.5], [0,0])
    plt.ylim(-0.7, 0.7)
    plt.ylabel("$\Delta$ log SFR [M$_\odot$/ year]", fontsize=18)
    plt.xlabel("Time [gyr]", fontsize=18)
    plt.title( "NIHAO Galaxies")
    plt.legend(prop={'size': 8}, loc='best')
    plt.savefig(str(fileBaseFalcon + "/capstoneFinal/plots/4_deltasfrg1.12e12"  + ".png"), dpi=300)
    plt.show()

all_delta_percentages = list(all_delta_after_log)

allMstar = np.array(allMstar)

all_delta_percentages = np.array(all_delta_percentages)
weighted = np.log10(allMstar + 1)  * all_delta_percentages
weighted_clean = weighted[np.logical_not(np.isnan(weighted))]
#this line below removes NaNs, right?
all_delta_percentages = all_delta_percentages[np.logical_not(np.isnan(all_delta_percentages))]

##### histogram ######
n, bins, patches = plt.hist(all_delta_percentages, bins=17)
A = np.vstack((np.digitize(all_delta_percentages, bins), allConc)).T
#res = {bins[int(i)]: np.mean(A[A[:, 0] == i, 1]) for i in np.unique(A[:, 0])}
res_lst = np.array([np.mean(A[A[:, 0] == i, 1]) for i in range(len(bins))])
res_lst[np.isnan(res_lst)] = min(concentrations)
res_lst = np.array(res_lst)
print('n', n)
print('bins', bins)

temp_min = 8.7
temp_max = 11
sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(temp_min, temp_max))
sm._A = []
for i in range(len(patches)):
   patches[i].set_facecolor(cmap(norm(res_lst[i],min(res_lst), max(res_lst) )))


plt.colorbar(sm)
plt.ylabel("Number of galaxy datapoints", fontsize=18)
plt.xlabel("Delta Log(SFR) [M$_\odot$/ Year]", fontsize=18)

#n, bins, patches = plt.hist(all_delta_percentages, color=cmap(norm(res_lst)), bins=51)
plt.title('Ellipticals')
plt.yscale('log')
plt.xlim(-3, 3)
#plt.xscale('symlog')
#plt.savefig(str(fileBaseFalcon + "ELL_hist.png"), dpi=300)
plt.show()

for i in avg_fit_dict:
    print(i, avg_fit_dict[i])

