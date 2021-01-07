import numpy as np
from matplotlib import pyplot as plt
from PhotochemPy import PhotochemPy


# Load input files
pc = PhotochemPy('../../input/templates/Hadean+HCN/species.dat', \
                 '../../input/templates/Hadean+HCN/reactions.rx', \
                 '../../input/templates/Hadean+HCN/planet.dat', \
                 '../../input/templates/Hadean+HCN/input_photchem.dat', \
                 '../../input/templates/Hadean+HCN/atmosphere.txt', \
                 '../../input/templates/Hadean+HCN/Sun_4.0Ga.txt')

# integrate to photochemical equilibirum
pc.integrate(nsteps=1000)

# plot
input = pc.in_dict()
out = pc.out_dict()
plt.rcParams.update({'font.size': 15})
fig,ax = plt.subplots(1,1,figsize=[9,5])
species = ['H2','CO','CH4','SO2','H2S']
colors = ['C0','C1','C2','C3','C4']
for i,sp in enumerate(species):
    ax.plot(out[sp],out['alt'],colors[i]+'-',label=sp)
    ax.plot(input[sp],input['alt'],colors[i]+'--')
ax.set_xscale('log')
ax.legend()
ax.set_ylabel('Altitude (km)')
ax.set_xlabel('Mixing Ratio')
ax.set_title('Solid lines = PhotochemPy\nDashed lines = old Atmos Photochem')
plt.savefig("Hadean+HCN_validation.pdf",bbox_inches='tight')
