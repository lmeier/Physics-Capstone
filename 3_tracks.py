import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import glob, os, pynbody, matplotlib
import pickle
import math
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import copy

falconBase = '/home/lem507/2018/pickles/'
#agngalDict = pickle.load(open('pickles/agnProperties.pkl', 'rb'))
#galDict =  pickle.load(open('pickles/nihaoProperties.pkl', 'rb'))
galDict =  pickle.load(open(falconBase + 'ALLjanNEW_NIHAOproperties.pkl', 'rb'))
agngalDict = pickle.load(open(falconBase + 'NIHAO_BHproperties.pkl', 'rb'))

galaxies = agngalDict.keys()

def monotonic_x(old_x, old_y):
    #check order, if looks reversed, reverse back
    if old_x[1] > old_x[-1]:
        old_x.reverse()
        old_y.reverse()

    highest_a = old_x[1]
    new_x = []
    new_y = []
    for a, b in zip(old_x, old_y):
        if a > highest_a:
            new_x.append(a)
            new_y.append(b)
            highest_a = a
    return new_x, new_y

galaxies = ['g2.79e12', 'g1.37e11', 'g5.38e11']

for g in galaxies:
    propDict = galDict[g]
    ms = propDict['ms']

    if g == 'g2.79e12':
        agnpropDict = agngalDict[g]
        agnms = agnpropDict['ms']
        agnsfr = agnpropDict['sfr']
        agnredshift= agnpropDict['z']
    #coldgas = propDict['coldgas']
    redshift = propDict['z']
    timeNew = propDict['time']
    sfr = propDict['sfr']

    SFR0 =  [29.3, 23.77, 6.80, 1.29]
    slope = [.9, .9, .9, 0.77]
    colors = ['y', 'g', 'b', 'k']
    num = 0
    fig = plt.figure(figsize=(10,10))

    expForZeroSFR = -3.4

    x = []
    y = []
    z = []
    print()
    for j in range(0, len(sfr)):
        if redshift[j] > 4:
            continue
        #check if it's within 0.05 of a whole number, label these points
        #print(redshift[j])
        z.append(redshift[j])
        x.append(ms[j])
        if sfr < 10**-3.4: #previously -2.5
            y.append(10**expForZeroSFR)
        else:
            y.append(sfr[j])
        grey = 1 - 0.9 #*random.random()

    x2 = []
    y2 = []
    z2 = []
    print()

    if g == 'g2.79e12':
        for j in range(0, len(agnsfr)):
            x2.append(agnms[j])
            z2.append(agnredshift[j])
            if sfr < 10**-3.4: #previously -2.5
                y2.append(10**expForZeroSFR)
            else:
                y2.append(agnsfr[j])

        a2, = plt.plot(x2,y2, 'o-', color='black', label='With AGN', alpha = 0.8)
        num =0
        for i in range(0, len(x)):
            if (abs(z[i] - int(z[i] + 0.1)) < 0.05) and (z[i] < 0.013 or z[i] > 0.96):
                print(z[i])
                print(num)
                plt.plot(x[i], y[i], linewidth=3,  marker='o',  color=colors[num], )
                num += 1
    #plt.plot(monotonic_x(x,y), 'o-', color=str(grey))
    #x, y =  monotonic_x(x,y)
    a1, = plt.plot(x,y, 'o-', color='red', label='Without AGN', alpha = 0.8)

    num = 0
    print(z)
    for i in range(0, len(x)):
        if (abs(z[i] - int(z[i] + 0.1)) < 0.05) and (z[i] < 0.013 or z[i] > 0.96):
            print(z[i])
            print(num)
            plt.plot(x[i], y[i], linewidth=3,  marker='o',  color=colors[num], )
            num +=1
    #x2, y2 = monotonic_x(x2,y2)
    #draw high density rectangle for red sequence
    x = [10**10.62, 10**10.88, 10**10.375, 10**10.1, 10**10.62]
    y = [10**-0.83, 10**-1.13, 10**-1.6, 10**-1.3, 10**-0.83]
    plt.plot(x,y, '-', color='r', linewidth=2, alpha =0.5)
    #draw high density triangle for blue cloud
    x = [ 10**10.05, 10**10.4, 10**8.7, 10**8.6, 10**8.1, 10**10.05] #3rd to last 10**8.25
    y = [ 10**0.5, 10**0.1, 10**-1.5, 10**-1.4, 10**-.95, 10**0.5] # 10**-1.25
    plt.plot(x,y, '-', color='r', linewidth=2, alpha = 0.5)
    #these are redshift 0.02 - 0.085
    #these are redshift 0.02 - 0.085
    for count, col in enumerate(colors[0:4]):
        fit_vals = np.asarray([10.0**4.0, 10.0**14]) #previous exponents were 6.0 and 12
        fit = SFR0[count]*((fit_vals/10**10)**slope[count])
        plt.plot(fit_vals, fit, color=col, linewidth=1.5, alpha=0.4)


    y_3_patch = mpatches.Patch(color='y', label='z = 3', alpha=0.4)
    y_2_patch = mpatches.Patch(color='g', label='z = 2', alpha=0.4)
    y_1_patch = mpatches.Patch(color='b', label='z = 1', alpha=0.4)
    y_0_patch = mpatches.Patch(color='k', label='z = 0', alpha=0.4)

    #agn_legend = mlines.Line2D([], [], color='k', marker='X',
    #markersize=15, label='AGN gals')
    if g == 'g2.79e12':
        plt.legend(handles=[a1, a2, y_3_patch, y_2_patch, y_1_patch, y_0_patch,], loc=2)
    else:
        plt.legend(handles=[a1, y_3_patch, y_2_patch, y_1_patch, y_0_patch,], loc=2)
    #ax.set_facecolor('w')
    plt.ylabel("SFR [M$_\odot$/ Year]", fontsize=18)
    plt.xlabel("M$_{\star}$ [M$_\odot$]", fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    #Renzini plus one in all directions
    plt.xlim(10**8, 10**12.5)
    plt.ylim(10**-2.5, 10**1.8)
    plt.title(g, fontsize=18)
    plt.savefig('plots/3_'+ g + 'tracks.png', dpi=300)
    plt.show()
    plt.close()


