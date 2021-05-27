import numpy as np
from PhotochemPy import PhotochemPy
import time
import os

template = '../../input/templates/Hadean+HCN'
star = 'Sun_4.0Ga.txt'

pc = PhotochemPy(template+'/species.dat', \
                 template+'/reactions.rx', \
                 template+'/settings.yaml', \
                 template+'/atmosphere.txt', \
                 template+'/'+star)

# integrate to photochemical equilibirum
start = time.time()
pc.integrate(nsteps=50)
end = time.time()

t = end-start

threads = os.getenv('OMP_NUM_THREADS')

fil = open('times.txt','a')
fil.write(threads+' '+'%.4e'%t+'\n')
fil.close()
