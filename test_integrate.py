import numpy as np
from matplotlib import pyplot as plt
from PhotochemPy import PhotochemPy
import time

pt = PhotochemPy('input/templates/Archean+haze/species.dat', \
                 'input/templates/Archean+haze/reactions.rx', \
                 'input/templates/Archean+haze/planet.dat', \
                 'input/templates/Archean+haze/input_photchem.dat', \
                 'input/templates/Archean+haze/atmosphere.txt', \
                 'input/templates/Archean+haze/Sun_2.7Ga.txt')


pt.integrate(nsteps=1000)




# print('CH4 flux =','%.2e'%pt.surf_flux()['CH4'])


plot = False

if plot:
    out = pt.out_dict()
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
