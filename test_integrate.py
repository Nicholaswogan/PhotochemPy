import numpy as np
from matplotlib import pyplot as plt
from PhotochemPy import PhotochemPy
import time

pt = PhotochemPy('input/species.dat', \
                 'input/reactions.rx', \
                 'input/planet.dat', \
                 'input/input_photchem.dat', \
                 'input/atmosphere.txt', \
                 'input/Sun_2.7Ga.txt')

start = time.time()
pt.integrate()
end = time.time()
print('Time to find equilibrium =',end-start,'seconds')

plot = True

if plot:
    out = pt.out_dict()
    plt.rcParams.update({'font.size': 15})
    fig,ax = plt.subplots(1,1,figsize=[9,5])
    specs = ['CH4']
    for sp in specs:
        ax.plot(out[sp],out['alt'],label=sp)

    ax.set_xscale('log')
    ax.legend()
    ax.set_ylabel('Altitude (km)')
    ax.set_xlabel('Mixing Ratio')
    # plt.savefig("Comparison.jpg",dpi = 100,bbox_inches='tight')
    plt.show()
