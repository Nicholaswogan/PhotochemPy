from PhotochemPy import PhotochemPy

template = 'Archean2Proterozoic'
sun = 'Sun_2.7Ga.txt'

pc = PhotochemPy('input/templates/'+template+'/species.dat', \
                 'input/templates/'+template+'/reactions.rx', \
                 'input/templates/'+template+'/PLANET.dat', \
                 'input/templates/'+template+'/input_photchem.dat', \
                 'input/templates/'+template+'/atmosphere.txt', \
                 'input/templates/'+template+'/'+sun)

converged = pc.integrate(method='Backward_Euler')
converged = pc.integrate(method='CVODE_BDF')
