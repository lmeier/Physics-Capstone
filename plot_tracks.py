import numpy as np
import matplotlib.pyplot as plt
import glob, os, pynbody, matplotlib
import pickle
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import copy


galDict = pickle.load(open('pickles/nihaoProperties.pkl', 'rb'))
galaxies = galDict.keys()
print(galaxies)
for g  in galaxies:

    propDict = galDict[g]

    ms = propDict['ms']
    coldgas = propDict['coldgas']
    redshift = propDict['z']
    timeNew = propDict['time']
    sfr = propDict['sfr']

    SFR0 =  [29.3, 23.77, 6.80, 1.29]
    slope = [.9, .9, .9, 0.77]
    colors = ['y', 'g', 'b', 'k']
    fig = plt.figure(figsize=(10,10))

    expForZeroSFR = -3.4

    x = []
    y = []
    print()
    for j in range(0, len(sfr)):
        x.append(ms[j])
        if sfr < 10**-3.4: #previously -2.5
            y.append(10**expForZeroSFR)
        else:
            y.append(sfr[j])
        grey = 1 - 0.9 #*random.random()
    #plt.plot(x,y, 'o-', color=str(grey))
    plt.plot(x,y, 'o-', color='r')

    #draw high density rectangle for red sequence
    x = [10**10.6, 10**10.7, 10**10.4, 10**10.3, 10**10.6]
    y = [10**-1.025, 10**-1.125, 10**-1.45, 10**-1.35, 10**-1.025]
    plt.plot(x,y, '-', color='r')

    #draw high density triangle for blue cloud
    x = [ 10**10.15, 10**10.4, 10**8.7, 10**8.6, 10**8.1, 10**10.15] #3rd to last 10**8.25
    y = [ 10**0.3, 10**0.1, 10**-1.5, 10**-1.4, 10**-.95, 10**0.3] # 10**-1.25
    plt.plot(x,y, '-', color='r')
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
    plt.legend(handles=[y_3_patch, y_2_patch, y_1_patch, y_0_patch,], loc=2)


    #ax.set_facecolor('w')
    plt.ylabel("SFR [M$_\odot$/ Year]", fontsize=18)
    plt.xlabel("M$_{\star}$ [M$_\odot$]", fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    #Renzini plus one in all directions
    plt.xlim(10**7, 10**(12.5))
    plt.ylim(10**-3.5, 10**2.0)
    plt.title(g, fontsize=18)
    plt.show()
    exit()
    plt.savefig('graphs/ms_sfr_feb/' + g +   '.png', dpi=300)
    plt.close()


