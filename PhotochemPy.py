import numpy as np
from Photochem import photochem

class PhotochemPy:
    def __init__(self,species_dat,reactions_dat):
        self.photo = photochem
        # In init we will allocate all memory to fortran
        # and get data (xsection, rates etc.) all loaded
        # into memory.

        self.nz = 200 # number of vertical grid points
        self.nw = 118 # number of wavelength bins

        # read species.dat to find nq, np, isl
        fil = open(species_dat,'r')
        lines = fil.readlines()
        fil.close()
        self.nq = 0
        self.isl = 0
        self.np = 0
        for line in lines:
            if line[0]=='*':
                pass
            else:
                if line.split()[1] == 'LL':
                    if line.split()[0].find('AER') != -1:
                        self.np+=1
                    self.nq+=1
                if line.split()[1] == 'SL':
                    self.isl+=1

        self.nsp = self.nq + self.isl

        fil = open(reactions_dat,'r')
        lines = fil.readlines()
        fil.close()
        self.nr = len(lines)
        # also find number of photolysis reactions and species (ks,kj)
        self.ks = 0 # number photo species
        self.kj = 0 # number photo reactions
        temp_spec = []
        for line in lines:
            if line[48:53] =='PHOTO':
                self.kj+=1
                tmp = 0
                for i in range(len(temp_spec)):
                    if line.split()[0]==temp_spec[i]:
                        tmp = 1
                if tmp == 0:
                    self.ks+=1
                    temp_spec.append(line.split()[0])


        self.photo.allocate_memory(self.nw,self.nz,self.nq,\
                                  self.np,self.nsp,self.nr,\
                                  self.ks,self.kj)
        self.photo.read_species(species_dat)
        self.photo.read_reactions(reactions_dat)

        # now we need to keep loading data and setting up the problem!!!
