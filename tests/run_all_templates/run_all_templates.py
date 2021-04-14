import numpy as np
from PhotochemPy import PhotochemPy

###############################
##### ModernEarth #####
###############################
template = '../../input/templates/ModernEarth'
star = 'Sun_now.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/planet.dat', \
                 template+'/input_photchem.dat', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')

###############################
##### Archean2Proterozoic #####
###############################
template = '../../input/templates/Archean2Proterozoic'
star = 'Sun_2.7Ga.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/planet.dat', \
                 template+'/input_photchem.dat', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')

########################
##### Archean+Haze #####
########################
template = '../../input/templates/Archean+haze'
star = 'Sun_2.7Ga.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/planet.dat', \
                 template+'/input_photchem.dat', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')

########################
###### Hadean+HCN ######
########################
template = '../../input/templates/Hadean+HCN'
star = 'Sun_4.0Ga.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/planet.dat', \
                 template+'/input_photchem.dat', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')

############################
##### Hadean+HCN_LLCO2 #####
############################
template = '../../input/templates/Hadean+HCN_LLCO2'
star = 'Sun_4.0Ga.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/planet.dat', \
                 template+'/input_photchem.dat', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')