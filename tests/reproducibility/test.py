import numpy as np
from PhotochemPy import PhotochemPy
import os

template = '../../input/templates/Archean+haze'
star = 'Sun_2.7Ga.txt'

pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/planet.dat', \
                 template+'/input_photchem.dat', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)

# integrate to photochemical equilibirum
pc.integrate()

np.savetxt('output.txt',pc.photo.usol_out)
