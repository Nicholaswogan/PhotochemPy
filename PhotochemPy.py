import numpy as np
from Photochem import photochem
import sys

class PhotochemPy:
    def __init__(self,species_dat,reactions_rx,planet_dat,\
                      photochem_dat, atmosphere_txt, flux_txt):
        self.photo = photochem
        self.species_dat = species_dat
        self.reactions_rx = reactions_rx
        self.planet_dat = planet_dat
        self.photochem_dat = photochem_dat
        self.atmosphere_txt = atmosphere_txt
        self.flux_txt = flux_txt

        self.nz = 200 # number of vertical grid points
        # read species.dat to find nq, np, isl
        fil = open(species_dat,'r')
        lines = fil.readlines()
        fil.close()
        self.nq = 0
        self.isl = 0
        self.np = 0
        inert = 0
        self.ispec = []
        for line in lines:
            if line[0]=='*':
                pass
            else:
                if line.split()[1] == 'LL':
                    if line.split()[0].find('AER') != -1:
                        self.np+=1
                    self.nq+=1
                    # record the species
                    self.ispec.append(line.split()[0])
                if line.split()[1] == 'SL':
                    self.isl+=1
                if line.split()[1] == 'IN':
                    inert += 1

        self.nsp = self.nq + self.isl + inert

        # Read reactions.rx to find nr, kj, and ks
        fil = open(reactions_rx,'r')
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

        self.photo.allocate_memory(self.nz,self.nq,\
                                  self.np,self.nsp,self.nr,\
                                  self.ks,self.kj)
        self.photo.setup(species_dat, \
                         reactions_rx, \
                         planet_dat, \
                         photochem_dat, \
                         atmosphere_txt, \
                         flux_txt)

        self.code_run = False

    def integrate(self,nsteps=1000):
        converged = self.photo.integrate(nsteps)
        if converged == 0:
            sys.exit('Photochemical model did not converge in '\
                     +str(nsteps)+' steps')
        else:
            self.code_run = True

    def out_dict(self):
        if not self.code_run:
            sys.exit('Need to integrate before outputting a solution!')
        elif self.code_run:
            out = {}
            out['alt'] = self.photo.z/1e5
            out['den'] = self.photo.den
            out['press'] = self.photo.p
            out['T'] = self.photo.t
            for i in range(self.nq):
                out[self.ispec[i]] = self.photo.usol_out[i,:]
            return out

    def surf_flux(self):
        if not self.code_run:
            sys.exit('Need to integrate before outputing surface flux!')
        elif self.code_run:
            out = {}
            for i in range(self.nq):
                out[self.ispec[i]] = self.photo.flow[i]
            return out

    def reset(self):
        self.photo.setup(self.species_dat, \
                         self.reactions_rx, \
                         self.planet_dat, \
                         self.photochem_dat, \
                         self.atmosphere_txt, \
                         self.flux_txt)
        self.code_run = False

    def out2in(self):
        if not self.code_run:
            sys.exit('Need to integrate before setting the output as input!')
        else:
            self.photo.usol_init = self.photo.usol_out
            if self.np > 0:
                self.photo.rpar_init = self.photo.rpar
                self.photo.wfall_init = self.photo.wfall
                self.photo.aersol_init = self.photo.aersol

    def set_surfflux(self,spec,flx):
        try:
            ind = self.ispec.index(spec)
        except:
            sys.exit('species not in the model')
        if self.photo.lbound[ind] == 2 or self.photo.lbound[ind] == 3:
            self.photo.sgflux[ind] = flx
        else:
            print('mbound set to',self.photo.lbound[ind],'so the flux will not change')

    def set_surfmr(self,spec,mix):
        try:
            ind = self.ispec.index(spec)
        except:
            sys.exit('species not in the model')
        if self.photo.lbound[ind] == 1:
            self.photo.fixedmr[ind] = mix
        else:
            print('mbound set to',self.photo.lbound[ind],'so the mixing ratio will not change')

    def set_lbound(self,spec,lbound):
        try:
            ind = self.ispec.index(spec)
            self.photo.lbound[ind] = lbound
        except:
            sys.exit('species not in the model')

    def set_mbound(spec,mbound):
        try:
            ind = self.ispec.index(spec)
            self.photo.mbound[ind] = mbound
        except:
            sys.exit('species not in the model')

    # include method that makes atmosphere.txt
