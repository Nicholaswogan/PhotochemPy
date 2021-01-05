from PhotochemPy import PhotochemPy

pc = PhotochemPy('input/templates/Archean+Haze/species.dat', \
              'input/templates/Archean+Haze/reactions.rx', \
              'input/templates/Archean+Haze/PLANET.dat', \
              'input/templates/Archean+Haze/input_photchem.dat', \
              'input/templates/Archean+Haze/atmosphere.txt', \
              'input/templates/Archean+Haze/Sun_2.7Ga.txt')

converged = pc.integrate()
