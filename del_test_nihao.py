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

mstar_nihao = []
sfr_nihao = []
time_nihao = []

galDict = pickle.load(open('/home/lem507/2018/pickles/novNEW_NIHAOproperties.pkl', 'rb'))
#galDict = pickle.load(open('/home/lem507/2018/pickles/AGN_ELL_NIHAOproperties.pkl', 'rb'))
galaxies = galDict.keys()
galaxies.sort()

timeMsDict = dict()
timeSfrDict = dict()


#nihaoGalDict = pickle.load(open('/home/lem507/2018/pickles/ALLnovNEW_NIHAOproperties.pkl', 'rb'))
#nihaoGalDict = pickle.load(open('/home/lem507/2018/pickles/ALLjanNEW_NIHAOpropertiesWORKINGCOPY.pkl', 'rb'))
nihaoGalDict = pickle.load(open('/home/lem507/2018/pickles/ALLjanNEW_NIHAOproperties.pkl', 'rb'))
#nihaoGalDict = pickle.load(open('/home/lem507/2018/pickles/AGN_ELL_NIHAOproperties.pkl', 'rb'))

for g in nihaoGalDict.keys():
    propDict = nihaoGalDict[g]
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
            #print("made it here")
        except:
            print(idx)
            print(len(timeNew))
            print(len(ms))
            print(len(redshift))
            ms_list = []
            ms_list.append(ms[idx])
            sfr_list = []
            sfr_list.append(sfr[idx])
            timeMsDict[t] = ms_list
            timeSfrDict[t] = sfr_list
            print(g)
for t, ms in timeMsDict.items():
    print(t)
    print(len(ms))


avg_fit_dict = dict()

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
    avg_fit_dict[t] = (m, b)
    fit_vals = np.asarray([10.0**4.0, 10.0**6, 10.0**10, 10.0**14])
    #plt.plot(fit_vals, myExpFunc(fit_vals, *popt), color='g')
    #plt.plot(fit_vals, myExpFunc(fit_vals, m, b), color='g')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()

    #if need to plot values can do this below, commenting out for now
    '''
    plt.scatter(ms_filt, sfr_filt)
    plt.scatter(mstar_nihao,sfr_nihao,color='red', alpha =0.1)
    plt.plot(fit_vals, linear_to_log_func(np.log10(fit_vals), m, b), color = 'g')
    plt.ylabel("SFR [M$_\odot$/ Year]", fontsize=18)
    plt.xlabel("M$_{\star}$ [M$_\odot$]", fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(10**7, 10**(12.5))
    plt.ylim(10**-3.5, 10**3.0)
    #fitted_vals = 10**((a)*np.log10(fit_vals) + (c))
    plt.title(t)
    plt.show()
    plt.clf()
    '''


plt.plot([0.0, 13.7], [0,0], color='k')
n = 1 #number of divisions


#['g1.05e11', 'g1.12e12', 'g1.37e11', 'g1.52e11', 'g1.92e12', 'g1.95e10', 'g2.04e11', 'g2.79e12', 'g3.06e11', 'g3.21e11', 'g4.86e10', 'g4.90e11', 'g5.38e11', 'g5.55e11', 'g6.77e10', 'g6.96e11', 'g7.08e11', 'g7.55e11', 'g8.26e11', 'g8.28e11', 'g9.59e10']

all_delta_percentages = []
all_delta_after_log = []


sfrEndAtTime = []
sfrEndConc = []




concDict = pickle.load(open('/home/lem507/2018/pickles/concentrations.pkl', 'rb'))
#concDict = pickle.load(open('/home/lem507/2018/pickles/concentrationsDMO_AGN_ELL.pkl', 'rb'))
concentrations = concDict.values()
#normConc = plt.colors.Normalize(concentrations)
allConc = []
allMstar = []
def norm_one(x):
    return ((x - min(concentrations)) / (max(concentrations) - min(concentrations)))

def norm(x, minimum, maximum):
    return ((x- minimum) / (maximum - minimum))

cmap = plt.get_cmap('plasma')
#colors = [cmap(i) for i in np.linspace(0, 1, number)]

#hacky code for getting a colorbar
'''
Z = [[0,0],[0,0]]
levels = np.arange(min(concentrations), max(concentrations))
CS3 = plt.contourf(Z, concentrations, cmap='viridis')
plt.clf()
'''
sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(vmin=min(concentrations), vmax=max(concentrations)))
# fake up the array of the scalar mappable. Urgh...
sm._A = []
#plt.colorbar(sm)


g1 = ['g1.95e10', 'g4.86e10',  'g1.05e11', 'g1.37e11', 'g1.52e11', 'g6.77e10','g9.59e10']
g2 = ['g2.04e11', 'g3.06e11', 'g3.21e11', 'g4.90e11', 'g5.38e11', 'g5.55e11', 'g6.96e11']
g3 = [ 'g7.55e11', 'g8.26e11', 'g8.28e11', 'g1.05e11', 'g1.12e12', 'g2.79e12',  'g7.08e11']

g4 = [g for g in nihaoGalDict.keys() if g not in g1+g2+g3]
'''
#for a in ['g3.59e12', 'g1.05e13', 'g7.92e12']:
#nihaoGalDict.pop(a)

g1 = [g for g in nihaoGalDict.keys()[0:6]]
g2 = [g for g in nihaoGalDict.keys()[6:12]]
g3 = [g for g in nihaoGalDict.keys()[12:]]
g4 = [g for g in nihaoGalDict.keys() if g not in g1+g2+g3]

'''
#gal_all = ['g1.95e10', 'g4.86e10', 'g6.77e10', 'g9.59e10', 'g1.05e11', 'g1.37e11', 'g1.52e11', 'g2.04e11', 'g3.06e11', 'g3.21e11', 'g4.90e11', 'g5.38e11', 'g5.55e11', 'g6.96e11',  'g7.55e11', 'g8.26e11', 'g8.28e11', 'g1.05e11', 'g1.12e12', 'g2.79e12']
gal_arr = [g1, g2, g3, g4]
#gal_arr =[gal_all]
size_names = ["1", "2", "3",]
#size_names = ["1"]
plot = [True, True, True, False]
#change this back1!!!
plot= [False, False, False, False]

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
        #concentration0 = "n/a"
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
        if g == 'g8.28e11':
            print("g8.28e11")
            print(timeNew)
            print(sfr)
            continue
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
        #a, b = smooth_for_plotting(timeNew, delta_sfr)

        #a, b = smooth_for_plotting(timeNew, delta_sfr_percentage)
        #previously was the one above
        a, b = smooth_for_plotting(timeNew, delta_after_log)
        a, b = smooth_for_plotting(a,b)
        label_long = g + ' ' +  str(conc0)[0:3]
        #c=B,cmap=cm.jet,vmin=0.,vmax=2.
        #linestyle = '--' if conc0 < 9 else '-'
        linestyle='-'
        plt.plot(a, b, color=cmap(norm(conc0, cur_min, cur_max)),  label=label_long,  ls=linestyle)

    #plt.clim(min(concentrations), max(concentrations))
    #pcm = ax.pcolormesh(concentrations, vmin=min(concentrations), vmax=max(concentrations), cmap='RdBu_r')
    #plt.colorbar(concentrations, cmap='seismic')
    if p == False:
        break


    sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(vmin=cur_min, vmax=cur_max))
    # fake up the array of the scalar mappable. Urgh...
    sm._A = []
    #plt.colorbar(sm)

    plt.colorbar(sm)
    plt.plot([0.0, 13.7], [0,0])
    #removing symlog for now
    #plt.yscale('symlog')
    #plt.ylim(-(10**1.5), 10**1.5)
    plt.ylabel("Delta SFR [M$_\odot$/ Year]", fontsize=18)
    plt.xlabel("Time [Gyr]", fontsize=18)
    plt.title( " NIHAO Galaxies")
    plt.legend(prop={'size': 8}, loc='best')

    #plt.savefig(str("del_sfr_ELL_plots/NIHAOdeltasfr" + name + ".png"), dpi=300)
    plt.show()
    #plt.clf()

#trying this out with new format
all_delta_percentages = list(all_delta_after_log)
#print(all_delta_percentages)

allMstar = np.array(allMstar)
#print(allMstar)
all_delta_percentages = np.array(all_delta_percentages)
weighted = np.log10(allMstar + 1)  * all_delta_percentages
weighted_clean = weighted[np.logical_not(np.isnan(weighted))]
#this line below removes NaNs, right?
all_delta_percentages = all_delta_percentages[np.logical_not(np.isnan(all_delta_percentages))]

##### first histogram ######

n, bins, patches = plt.hist(all_delta_percentages, bins=23) #15 is good, as is 23
A = np.vstack((np.digitize(all_delta_percentages, bins), allConc)).T
#res = {bins[int(i)]: np.mean(A[A[:, 0] == i, 1]) for i in np.unique(A[:, 0])}
res_lst = np.array([np.mean(A[A[:, 0] == i, 1]) for i in range(len(bins))])
res_lst[np.isnan(res_lst)] = min(concentrations)
res_lst = np.array(res_lst)
print('n', n)
print('bins', bins)
'''
print('a')
#cmap = plt.set_cmap(sm)
#sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=min(concentra    tions), vmax=max(concentrations)))
sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(vmin=min(res_lst), vmax=max(res_lst)))
sm._A = []
print('a')
'''

temp_min = 8
temp_max = 11
sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(temp_min, temp_max))
sm._A = []
print('a')

for i in range(len(patches)):
   patches[i].set_facecolor(cmap(norm(res_lst[i],min(res_lst), max(res_lst) )))


plt.colorbar(sm)
plt.ylabel("Number of galaxy datapoints", fontsize=18)
plt.xlabel("Delta Log(SFR) [M$_\odot$/ Year]", fontsize=18)

#n, bins, patches = plt.hist(all_delta_percentages, color=cmap(norm(res_lst)), bins=51)
plt.yscale('log')
plt.title("NIHAO w/o AGN")
#plt.xscale('symlog')
#plt.savefig(str("del_sfr_ELL_plots/NIHAO_hist.png"), dpi=300)
plt.show()

for i in avg_fit_dict:
    print(i, avg_fit_dict[i])



plt.scatter(sfrEndAtTime, sfrEndConc)
plt.show()
exit()
##### ORIGINAL TRY BELOW THIS

print(all_delta_percentages)

allMstar = np.array(allMstar)
print(allMstar)
all_delta_percentages = np.array(all_delta_percentages)
weighted = np.log10(allMstar + 1)  * all_delta_percentages
weighted_clean = weighted[np.logical_not(np.isnan(weighted))]
all_delta_percentages = all_delta_percentages[np.logical_not(np.isnan(all_delta_percentages))]

##### first histogram ######

n, bins, patches = plt.hist(all_delta_percentages, bins=51) #original 51
A = np.vstack((np.digitize(all_delta_percentages, bins), allConc)).T
#res = {bins[int(i)]: np.mean(A[A[:, 0] == i, 1]) for i in np.unique(A[:, 0])}
res_lst = np.array([np.mean(A[A[:, 0] == i, 1]) for i in range(len(bins))])
res_lst[np.isnan(res_lst)] = min(concentrations)
res_lst = np.array(res_lst)


#filtering!!!!!
print("n", n)
'''
print("res_lst", res_lst)
new_n = []
new_res_lst = []
for a, b in zip(n, res_lst):
    if a > 2:
        new_n.append(a)
        new_res_lst.append(b)
    else:
        new_n.append(0)
        new_res_lst.append(np.mean(concentrations))
n = new_n
res_lst = new_res_lst
'''
print('n', n)
print('bins', bins)
print("res_lst", res_lst)
print('a')
#cmap = plt.set_cmap(sm)
#sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=min(concentra    tions), vmax=max(concentrations)))
# old are 5, 12
temp_min = 8
temp_max = 11
sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(temp_min, temp_max))
sm._A = []
print('a')
#for i in range(len(patches)):
for i, (n_i, res) in enumerate(zip(n, res_lst)):
    if n_i > 2 or i == 32 or i == 33:
        patches[i].set_facecolor(cmap(norm(res_lst[i],temp_min, temp_max)))
    else:
        print(n_i)
        print(i)
        patches[i].set_color('w')
        #patches[i].set_facecolor(cmap(norm(res_lst[i],min(res_lst), max(res_lst) )))

print(g4)
plt.xlim(-15, 15)
plt.colorbar(sm)
plt.ylabel("Number of galaxy datapoints", fontsize=18)
plt.xlabel("Delta Log(SFR) [M$_\odot$/ Year]", fontsize=18)
#plt.xlim(-2, 2)

#n, bins, patches = plt.hist(all_delta_percentages, color=cmap(norm(res_lst)), bins=51)
plt.yscale('log')
#plt.xscale('symlog')
#plt.savefig(str("del_sfr_ELL_plots/FEBdeltasfrHIST_21st.png"), dpi=300)
plt.show()

print(res_lst)
print(len(res_lst))

'''
x = sfrEndAtTime
y = sfrEndConc
plt.hist2d(x, y, bins=100)
plt.xlabel('x')
plt.ylabel('y')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
plt.show()
'''


for i in avg_fit_dict:
    print(i, avg_fit_dict[i])
'''
###### second histogram, weighted this time!

n, bins, patches = plt.hist(weighted_clean, bins=51)
A = np.vstack((np.digitize(weighted_clean, bins), allConc)).T
#res = {bins[int(i)]: np.mean(A[A[:, 0] == i, 1]) for i in np.unique(A[:, 0])}
res_lst = np.array([np.mean(A[A[:, 0] == i, 1]) for i in range(len(bins))])
res_lst[np.isnan(res_lst)] = min(concentrations)
res_lst = np.array(res_lst)
print('n', n)
print('bins', bins)

print('a')
#cmap = plt.set_cmap(sm)
#sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=min(concentra    tions), vmax=max(concentrations)))
sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(vmin=min(res_lst), vmax=max(res_lst)))
sm._A = []
print('a')

for i in range(len(patches)):
   patches[i].set_facecolor(cmap(norm(res_lst[i],min(res_lst), max(res_lst) )))

plt.title("Weighted Histogram")
plt.colorbar(sm)
plt.ylabel("Number of galaxy datapoints", fontsize=18)
plt.xlabel("Delta Log(SFR) * Log(Mstar) [M${_\odot}^2$/ Year]", fontsize=18) #fix units

#n, bins, patches = plt.hist(all_delta_percentages, color=cmap(norm(res_lst)), bins=51)
plt.yscale('log')
#plt.xscale('symlog')
plt.savefig(str("deltasfrHISTweighted_21st.png"), dpi=300)
plt.show()

print(res_lst)
'''
'''
plt.hist2d(all_delta_percentages, allConc, bins=51)
plt.show()
'''
