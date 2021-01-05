from PhotochemPy import PhotochemPy

pc = PhotochemPy('input/templates/Archean+haze/species.dat', \
              'input/templates/Archean+haze/reactions.rx', \
              'input/templates/Archean+haze/PLANET.dat', \
              'input/templates/Archean+haze/input_photchem.dat', \
              'input/templates/Archean+haze/atmosphere.txt', \
              'input/templates/Archean+haze/Sun_2.7Ga.txt')





converged = pc.integrate()
