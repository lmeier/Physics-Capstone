import pynbody
import numpy as np
import pylab as plt
import pynbody.analysis.profile as profile
import glob, os, pickle
import math
import sys, ast




filename = '/scratch/database/nihao/nihao_classic/'
galaxies = ['g1.12e12', 'g2.79e12']
galaxies = ['g4.90e11', 'g8.26e11', 'g5.38e11', 'g1.37e11']
if len(sys.argv) == 1:
    gal = raw_input('What galaxy do you want to plot')
else:
    gal = sys.argv[1]

plot_bh = raw_input('Plot bh? Y for yes, anything else for no. ')
#print('here already')
#propDict = {'redshiftarray': redshiftarray, 'gasarray': gasarray, 'dm': dmarray, 'star':stararray, 'sfr': sfrarray, 'outflowarray': outflowarray, 'inflowarray': inflowarray, 'surfden1kpc': surfden1kpc, 'ssfr': ssfr, 'smass10': smass10, 'smass': smass}

#with open('pickles/6_' + gal + '.pkl', 'wb') as fp:

pD = pickle.load(open('pickles/6_' + gal + '.pkl', 'rb'))


fig, ax1 = plt.subplots()
plt.title(gal, fontsize=18)
plt.xlabel("z", fontsize=18)
ax1.plot(pD['redshiftarray'], pD['star'], 'r', linewidth=2.0, label="stars")
ax1.plot(pD['redshiftarray'], pD['dm'], 'b', linewidth=2.0,  label="dm")
ax1.plot(pD['redshiftarray'], pD['gasarray'], 'g', linewidth=2.0, label="gas")
if plot_bh == 'y':
    ax1.plot(pD['redshiftarray'], pD['bharray'], 'k', linewidth=2.0, label="bh")
ax1.grid(False)
ax1.set_ylim(10**5, 10**11)
plt.legend(loc=2, fontsize="medium", frameon=False)
plt.xlim(6, 1)
ax1.semilogy()
ax2 = ax1.twinx()
ax2.plot(pD['redshiftarray'], pD['sfr'], color="purple", label="sfr")
ax2.plot(pD['redshiftarray'], pD['inflowarray'], color="cyan", label="inflow")
ax2.plot(pD['redshiftarray'], pD['outflowarray'], color="magenta", label="outflow")
plt.legend(loc=4, fontsize="medium", frameon=False)
ax2.set_ylim(10**(-1.5), 10**4.5)
ax2.semilogy()
ax1.set_ylabel(r"log M [[$M_{\odot}$]", fontsize=18)
ax2.set_ylabel(r"log Mdot [$M_{\odot}$ / year]", fontsize=18)
plt.savefig('plots/6_' + gal + 'mt.png', dpi=300)
plt.show()


surfden1kpclog10 = []
ssfrlog10 = []

for i in pD['surfden1kpc']:
    surfden1kpclog10.append(np.log10(0.001 + i))
for i in pD['ssfr']:
    ssfrlog10.append(np.log10(0.001 + i))


fig, ax = plt.subplots()
plt.plot(surfden1kpclog10, ssfrlog10, '-o', color='k')
plt.xlabel(r'log $\sum_i$ [$M_{\odot}$/kpc$^2$]', fontsize=18)
plt.ylabel(r'log sSFR [Gyr$^-1$]', fontsize=18)
plt.xlim(7, 11)
plt.ylim(max(ssfrlog10) + 0.5, min(ssfrlog10) - 0.5)
ax.grid(True)
count = 0
last = 5
if len(sys.argv) == 3:
    list_plot_counts = ast.literal_eval(sys.argv[2])

for txt in pD['redshiftarray']:
    if len(sys.argv) == 3:
        if count in list_plot_counts:
            ax.annotate('z =' + str(round(txt, 1)), xy=((surfden1kpclog10[count]) + 0.1, ssfrlog10[count] + 0.1  ))
    else:
        if txt > 4 or txt < 0.01 or ((txt > 1.6 or txt < 0.7) and  (last - txt) > 0.5):
            ax.annotate('z =' + str(round(txt, 1)), xy=((surfden1kpclog10[count]) + 0.1, ssfrlog10[count] + 0.1  ))
            last = txt
    count += 1
plt.savefig('plots/6_' + gal + 'ssfr.png', dpi=300)
plt.show()
plt.close()


