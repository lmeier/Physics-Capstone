import pynbody
import numpy as np
import pylab as plt
import pynbody.analysis.profile as profile
import glob, os, pickle
import math
dir_out = '/home/lem507/compaction/'

filename = '/home/lem507/galaxies/g1.12e12/'

redshift = []
surfden1kpc = []
# surfdenre = [] #will use for second plots, stellar mass within Re divided by Re squared --> 1/2 sum_mass / Re^2
ssfr = []

datafiles = sorted(glob.glob(filename+'g1.12e12.0????'))
for d in datafiles:
    sim = pynbody.load(d)
    h = sim.halos()
    h1 = h[1]
    h1.physical_units()

    redshift.append(sim.properties['z'])
    ages = sim.star['age'].in_units('Gyr')
    sizes = sim.star['mass'].in_units('Msol')

    #sum_mass = sum(sizes) + 1 #add 1 to avoid divide by zero error later
    newstarmass = 0

    for age, size in zip(ages, sizes):
        if age <= 0.05: #if age is less than 50 Myr
            newstarmass += size


    cen = pynbody.analysis.halo.center(h1,mode='ssc', retcen=True)
    for j, arr in enumerate(['x', 'y', 'z']):
        sim[arr] -= cen[j]
    vcen = pynbody.analysis.halo.vel_center(h1,retcen=True,cen_size='3 kpc')
    for j, arr in enumerate(['vx','vy','vz']):
        sim[arr] = -vcen[j]
    ps10 = profile.Profile(h1.s, ndim=3, min='.001 kpc', max = '10 kpc')
    ps = profile.Profile(h1.s, ndim=3, min='.001 kpc', max='1 kpc')
    sum_mass = ps10['mass'].sum() + 1 # +1 to avoid dvide by 0, only affects early redshifts
    ssfr.append( (newstarmass / 0.05 / sum_mass))

    surfden1kpc.append(ps['mass'].sum() / math.pi)


    #compute halo angular momentum
    js = pynbody.analysis.angmom.ang_mom_vec(h1)
    trans = pynbody.analysis.angmom.calc_faceon_matrix(js)
    #rotate the whole simulations so that z goes in the direction of the halo angular mome$
    sim.transform(trans)


surfden1kpclog10 = []
ssfrlog10 = []

for i in surfden1kpc:
    surfden1kpclog10.append(np.log10(0.001 + i))
for i in ssfr:
    ssfrlog10.append(np.log10(0.001 + i))





fig, ax = plt.subplots()
plt.plot(surfden1kpclog10, ssfrlog10, '-0')
plt.xlabel(r'log $\sum_i$ [$M_{\odot}$/kpc$^2$]')
plt.ylabel(r'log sSFR [Gyr$^-1$]')
plt.xlim(7, 11)
plt.ylim(max(ssfrlog10), min(ssfrlog10))
ax.grid(True)
count = 0
for txt in redshift:
    if txt > 1.5 or (last - txt) > 0.3:
        ax.annotate('z =' + str(round(txt, 1)), xy=((surfden1kpclog10[count]), (ssfrlog10[count])))
    last = txt

plt.savefig(dir_out+'g1.12e12compaction2kpc.png')
plt.close()
