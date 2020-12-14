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
photochem.read_species('input/speciesOG.dat')
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
for i in range(100):
    prates = photochem.photo(photochem.zy,photochem.agl,photochem.io2\
                    ,photochem.ino,photochem.usol_init,photochem.kj)
end = time.time()
print('photo',(end-start)/10)

start = time.time()
photochem.rainout(22,0,photochem.usol_init)
end = time.time()


start = time.time()
photochem.aercon(photochem.usol_init)
end = time.time()


start = time.time()
conver = photochem.sedmnt(photochem.frak,photochem.hcdens, photochem.ihztype,photochem.nz,photochem.np)
end = time.time()
