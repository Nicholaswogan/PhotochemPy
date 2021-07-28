from PhotochemPy import PhotochemPy

#######################
##### ModernEarth #####
#######################
template = '../input/templates/ModernEarth'
star = 'Sun_now.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/settings.yaml', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')

###############################
##### Archean2Proterozoic #####
###############################
template = '../input/templates/Archean2Proterozoic'
star = 'Sun_2.7Ga.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/settings.yaml', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')

########################
##### Archean+Haze #####
########################
template = '../input/templates/Archean+haze'
star = 'Sun_2.7Ga.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/settings.yaml', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')

########################
###### Hadean+HCN ######
########################
template = '../input/templates/Hadean+HCN'
star = 'Sun_4.0Ga.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/settings.yaml', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')

##########################
##### MarsModern+Cl ######
##########################
template = '../input/templates/MarsModern+Cl'
star = 'Sun_now.txt'
pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/settings.yaml', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)
pc.integrate(method='Backward_Euler')
pc.integrate(method='CVODE_BDF')

####################
###### Saturn ######
####################
template = '../input/templates/Saturn'
star = 'Sun_now.txt'
pc = PhotochemPy(template+'/species.dat', \
                 '../input/templates/HadeanImpact'+'/reactions_new.rx', \
                 template+'/settings_saturn.yaml', \
                 template+'/atmosphere_saturn.txt', \
                 template+'/'+star)
# pc.integrate(method='Backward_Euler') doesn't work
pc.integrate(method='CVODE_BDF')
