import numpy as np
from matplotlib import pyplot as plt

# Import PhotochemPy. Do this one of the following ways!
from PhotochemPy import PhotochemPy # If installed with pip use this!

# If compiled in root directory with compile.sh, uncomment this
# (and comment out the other import command)
# import sys
# sys.path.append("..")
# from PhotochemPy import PhotochemPy

# Load input files
pc = PhotochemPy('../input/templates/Archean+haze/species.dat', \
                 '../input/templates/Archean+haze/reactions.rx', \
                 '../input/templates/Archean+haze/planet.dat', \
                 '../input/templates/Archean+haze/input_photchem.dat', \
                 '../input/templates/Archean+haze/atmosphere.txt', \
                 '../input/templates/Archean+haze/Sun_2.7Ga.txt')

# integrate to photochemical equilibirum
pc.integrate(nsteps=1000)

# Plot
plot = True
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
    plt.show()
