import numpy as np
from Photochem import photochem
import time

nz = 200
nq = 61
nnp = 4
nsp = 74
nr = 392
ks = 33
kj = 60
photochem.allocate_memory(nz,nq,nnp,nsp,nr,ks,kj)
# now read input files
photochem.read_species('input/species.dat')
photochem.read_reactions('input/reactions.rx')
photochem.read_planet('input/planet.dat')
photochem.read_photochem('input/input_photchem.dat')
photochem.read_atmosphere('input/atmosphere.txt')
photochem.photgrid(100.0e5)
photochem.densty()
photochem.rates()
photochem.difco()
H2O = photochem.photsatrat(22,photochem.nz)
fval = photochem.dochem(-1,22,photochem.isl,photochem.usol_init)
photochem.ltning(photochem.usol_init)
photochem.aertab()
photochem.initphoto('input/Sun_2.7Ga.txt')
photochem.flux = np.loadtxt('flux.dat')
photochem.initmie(photochem.nw,photochem.wavl,1,0)

# in the loop
start = time.time()
for i in range(1):
    prates = photochem.photo(photochem.zy,photochem.agl,photochem.io2\
                    ,photochem.ino,photochem.usol_init,photochem.kj)
end = time.time()
print('photo',(end-start)/10)

data = np.loadtxt('prates.dat')
for i in range(kj):
    print(i+1,'%.4e'%(prates[i,-1]/data[i]),'%.4e'%(prates[i,-1]),'%.4e'%(data[i]))
