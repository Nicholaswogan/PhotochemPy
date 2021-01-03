import numpy as np
from matplotlib import pyplot as plt
from PhotochemPy import PhotochemPy
import time
import sys


pc = PhotochemPy('input/templates/ModernEarth/species.dat', \
              'input/templates/ModernEarth/reactions.rx', \
              'input/templates/ModernEarth/planet.dat', \
              'input/templates/ModernEarth/input_photchem.dat', \
              'input/templates/ModernEarth/atmosphere.txt', \
              'input/templates/ModernEarth/Sun_now.txt')


converged = pc.integrate(nsteps=500)


plot = False


# print(pc.out_dict()['CO2'])
if plot:
    out = pc.out_dict()
    plt.rcParams.update({'font.size': 15})
    fig,ax = plt.subplots(1,1,figsize=[9,5])
    specs = ['CH4','CO','O2','H2']
    for sp in specs:
        ax.plot(out[sp],out['alt'],label=sp)

    ax.set_xscale('log')
    ax.legend()
    ax.set_ylabel('Altitude (km)')
    ax.set_xlabel('Mixing Ratio')
    # plt.savefig("Comparison.jpg",dpi = 100,bbox_inches='tight')
    plt.show()
