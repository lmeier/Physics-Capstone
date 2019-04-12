import pynbody
import numpy as np
import pylab as plt
import pynbody.analysis.profile as profile
import glob, os, pickle
import math

filename = '/scratch/database/nihao/nihao_classic/'
galaxies = ['g1.12e12', 'g2.79e12']
gal = galaxies[0]

redshiftarray = []
gasarray = []
dmarray = []
stararray = []
sfrarray = []
inflowarray = []
outflowarray = []

surfden1kpc = []
ssfr = []

datafiles = sorted(glob.glob(filename+gal+'/' + gal + '.0????'))
print(datafiles)
print("printed datafiles above")
for d in datafiles:
    sim = pynbody.load(d)
    if sim.properties['z'] > 5:
        continue
    h = sim.halos()
    h1 = h[1]
    h1.physical_units()
    pynbody.analysis.halo.center(h1, mode='hyb')
    pynbody.analysis.angmom.faceon(h1)

    newstarmass = 0
    stars, darkmatter, gas, inflow, outflow  = 0, 0, 0, 0, 0
    for pos, mass, age  in zip(sim.star['pos'].in_units('kpc'), sim.star['mass'].in_units('Msol'), sim.star['age'].in_units('yr')):
        if (pos[0]**2 + pos[1]**2 + pos[2]**2) < 1:
            stars += mass
        if math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2) < 1.1 and math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2) > 0.9:
            if age <= (60 * (10**6)):
                newstarmass += mass

    for pos, mass  in zip(sim.dm['pos'].in_units('kpc'), sim.dm['mass'].in_units('Msol')):
        if (pos[0]**2 + pos[1]**2 + pos[2]**2) < 1:
            darkmatter += mass
    for pos,  mass, vel  in zip(sim.gas['pos'].in_units('kpc'), sim.gas['mass'].in_units('Msol'), sim.gas['vr'].in_units("kpc yr**-1")):
        if (pos[0]**2 + pos[1]**2 + pos[2]**2) < 1:
            gas += mass
        if math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2) < 1.1 and math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2) > 0.9:
            if vel > 0:
                outflow += (mass * vel)
            else:
                inflow += (mass * vel * -1)

    gasarray.append(gas)
    stararray.append(stars)
    dmarray.append(darkmatter)
    sfrarray.append( (newstarmass / (60 * (10**6))) )
    outflowarray.append(outflow / 0.2)
    inflowarray.append(inflow / 0.2)
    redshiftarray.append(sim.properties['z'])


    smass10 = profile.Profile(h1.s, ndim=3, min='.001 kpc', max = '10 kpc')['mass'].sum()
    smass = profile.Profile(h1.s, ndim=3, min='.001 kpc', max='1 kpc')['mass'].sum()
    surfden1kpc.append(smass / math.pi)
    ssfr.append(newstarmass / 0.05 / smass)

    print sim.properties['z']









print('here already')
propDict = {'redshiftarray': redshiftarray, 'gasarray': gasarray, 'dm': dmarray, 'star':stararray, 'sfr': sfrarray, 'outflowarray': outflowarray, 'inflowarray': inflowarray, 'surfden1kpc': surfden1kpc, 'ssfr': ssfr, 'smass10': smass10, 'smass': smass}

with open('pickles/6_' + gal + '.pkl', 'wb') as fp:
    pickle.dump(propDict, fp)




fig, ax1 = plt.subplots()
plt.title(galaxy + "     1 kpc")
plt.xlabel("z")
ax1.plot(redshiftarray, stararray, 'r', linewidth=2.0, label="stars")
ax1.plot(redshiftarray, dmarray, 'k', linewidth=2.0,  label="dm")
ax1.plot(redshiftarray, gasarray, 'b', linewidth=2.0, label="gas")
ax1.grid(False)
ax1.set_ylim(10**5, 10**11)
plt.legend(loc=2, fontsize="small", frameon=False)
plt.xlim(6, 1)
ax1.semilogy()
ax2 = ax1.twinx()
ax2.plot(redshiftarray, sfrarray, color="purple", label="sfr")
ax2.plot(redshiftarray, inflowarray, color="cyan", label="inflow")
ax2.plot(redshiftarray, outflowarray, color="magenta", label="outflow")
plt.legend(loc=4, fontsize="small", frameon=False)
ax2.set_ylim(10**(-1.5), 10**4.5)
ax2.semilogy()
ax1.set_ylabel(r"log M [[$M_{\odot}$]")
ax2.set_ylabel(r"log Mdot [$M_{\odot}$ / year]")
plt.savefig('plots/' + gal + 'mt.png', dpi=300)



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

plt.savefig('plots/' + gal + 'ssfr.png', dpi=300)
plt.close()


