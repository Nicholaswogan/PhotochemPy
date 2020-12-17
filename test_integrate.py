import numpy as np
from matplotlib import pyplot as plt
from PhotochemPy import PhotochemPy
import time

pt = PhotochemPy('input/templates/Hadean+HCN/species.dat', \
                 'input/templates/Hadean+HCN/reactions.rx', \
                 'input/templates/Hadean+HCN/planet.dat', \
                 'input/templates/Hadean+HCN/input_photchem.dat', \
                 'input/templates/Hadean+HCN/atmosphere.txt', \
                 'input/templates/Hadean+HCN/Sun_4.0Ga.txt')

start = time.time()
pt.integrate()
end = time.time()
print('Time to find equilibrium =',end-start,'seconds')


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
