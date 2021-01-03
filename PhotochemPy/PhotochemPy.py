import numpy as np
from Photochem import photochem
import sys
import os

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'
photochem.rootdir = "{:500}".format(rootdir)


class PhotochemPy:
    '''
    The PhotochemPy class.

    :ivar photo: (object) Compiled Fortran module "photochem". It has many methods and attributes.
    :ivar nq: Number of long lived species
    :ivar np: Number of particles
    :ivar isl: Number of short lived species
    :ivar ispec: List if species names
    :ivar nsp: Total number of species
    :ivar nr: Number of reactions
    :ivar ks: Number of photolysis species
    :ivar kj: Number of photolysis reactions
    :ivar species_dat: Name of the input species file.
    :ivar reactions_rx: Name of the input reactions file.
    :ivar planet_dat: Name of the input planet file.
    :ivar photochem_dat: Name of the input photochem file.
    :ivar atmosphere_txt: Name of the input atmosphere file.
    :ivar flux_txt: Name of the input solar flux file.
    :ivar code_run: If True/False then code has converged/ has not converged to equilrium.

    To import and initialize the PhotochemPy class do the following

    .. highlight:: python
    .. code-block:: python

        from PhotochemPy import PhotochemPy
        pc = PhotochemPy(species_dat, reactions_rx, planet_dat, photochem_dat, atmosphere_txt, flux_txt)

    Parameters
    ----------
    species_dat : string
        Path to input file describing the species in the photochemical model,
        and their boundary conditions.
    reactions_rx : string
        Path to input file describing the reactions in the atmosphere and
        their rates.
    planet_dat : string
        Path to input file describing the planet (surface gravity, etc.)
    photochem_dat : string
        Path to input file setting a few options for the model.
    atmosphere_txt : string
        Path to input file describing the initial atmospheric composition,
        temperature structure, eddy diffusion profile, and aersol parameters.
    flux_txt : string
        Path to input file describing the stellar flux

    '''
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
        '''
        Integrates atomsphere to photochemical equilibrium using the backward
        Euler method.

        Parameters
        ----------
        nsteps : integer, optional
            The number of steps the integrator takes to find photochemical
            equilibrium. The default value is 1000.

        Returns
        -------
        converged : bool
            If True, then the code converged to equilibrium. If False,
            the code did not converge.
        '''
        converged = self.photo.integrate(nsteps)
        if converged == 0:
            self.code_run = False
        else:
            self.code_run = True

        converged = self.code_run
        return converged

    def out_dict(self):
        '''
        Makes a dictionary of the atmosphere after integration to photochemical
        equilibrium

        Returns
        -------
        out : dict
            Dictionary containing the mixing ratio of all species in the atmosphere,
            temperature structure, total pressure, and total number density.

        Raises
        ------
        SystemExit
            When photochemical model has not been integrated to equilibrium.
        '''
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

    def in_dict(self):
        '''
        Makes a dictionary of the atmosphere before integration to photochemical
        equilibrium. This is typically the atmosphere described in the input file
        atmosphere_txt, unless the input atmosphere has been changed with the out2in
        method.

        Returns
        -------
        out : dict
            Dictionary containing the mixing ratio of all species in the input atmosphere,
            temperature structure, total pressure, and total number density.
        '''
        out = {}
        out['alt'] = self.photo.z/1e5
        out['den'] = self.photo.den
        out['press'] = self.photo.p
        out['T'] = self.photo.t
        for i in range(self.nq):
            out[self.ispec[i]] = self.photo.usol_init[i,:]
        return out

    def surf_flux(self):
        '''
        Makes dictionary of the surface fluxes of each species at photochemical
        equilibrium.

        Returns
        -------
        out : dict
            Surface flux of each species in the model in molecules/cm2/s. Positive
            flux means a flux into the atmosphere.

        Raises
        ------
        SystemExit
            When photochemical model has not been integrated to equilibrium.
        '''
        if not self.code_run:
            sys.exit('Need to integrate before outputing surface flux!')
        elif self.code_run:
            out = {}
            for i in range(self.nq):
                out[self.ispec[i]] = self.photo.flow[i]
            return out

    def reset(self):
        '''
        Resets the problem by reading in the original input files (e.g. species_dat, ...)
        '''
        self.photo.setup(self.species_dat, \
                         self.reactions_rx, \
                         self.planet_dat, \
                         self.photochem_dat, \
                         self.atmosphere_txt, \
                         self.flux_txt)
        self.code_run = False

    def out2in(self):
        '''
        Changes the initial atmosphere to the output of the past integration.

        Raises
        ------
        SystemExit
            When photochemical model has not been integrated to equilibrium.
        '''
        if not self.code_run:
            sys.exit('Need to integrate before setting the output as input!')
        else:
            self.photo.usol_init = self.photo.usol_out
            if self.np > 0:
                self.photo.rpar_init = self.photo.rpar
                self.photo.wfall_init = self.photo.wfall
                self.photo.aersol_init = self.photo.aersol

    def set_surfflux(self,spec,flx):
        '''
        Sets surface flux boundary condition for a species.
        self.photo.lbound must be 2 or 3.

        Parameters
        ----------
        spec : string
            Species to change surface flux
        flx : float
            Surface flux of species spec in molecules/cm2/s

        Raises
        ------
        SystemExit
            When species spec is not in the model
        '''
        try:
            ind = self.ispec.index(spec)
        except:
            sys.exit('species not in the model')
        if self.photo.lbound[ind] == 2 or self.photo.lbound[ind] == 3:
            self.photo.sgflux[ind] = flx
        else:
            print('mbound set to',self.photo.lbound[ind],'so the flux will not change')

    def set_mr(self,spec,mix):
        '''
        Sets surface mixing ratio boundary condition for a species.
        self.photo.lbound must be 1.

        Parameters
        ----------
        spec : string
            Species to change surface mixing ratio
        mix : float
            Surface mixing ratio of species spec.

        Raises
        ------
        SystemExit
            When species spec is not in the model
        '''
        try:
            ind = self.ispec.index(spec)
        except:
            sys.exit('species not in the model')
        if self.photo.lbound[ind] == 1:
            self.photo.fixedmr[ind] = mix
        else:
            print('mbound set to',self.photo.lbound[ind],'so the mixing ratio will not change')

    def set_lbound(self,spec,lbound):
        '''
        Sets the type of lower boundary condition

        Parameters
        ----------
        spec : string
            Species to change lbound
        lbound : integer
            New lower boundary condition type.

            - 0 = constant deposition velocity (VDEP)

            - 1 = constant mixing ratio

            - 2 = constant upward flux (SGFLUX)

            - 3 = constant vdep + vertically distributed upward flux

        Raises
        ------
        SystemExit
            When species spec is not in the model or when all(lbound != [0,1,2,3])
        '''
        if all(lbound!=np.array([0,1,2,3])):
            sys.exit('lbound must be 0, 1, 2 or 3')
        try:
            ind = self.ispec.index(spec)
            self.photo.lbound[ind] = lbound
        except:
            sys.exit('species not in the model')

    def set_mbound(self,spec,mbound):
        '''
        Sets the type of upper boundary condition. This is not used often.

        Parameters
        ----------
        spec : string
            Species to change lbound
        mbound : integer
            New lower boundary condition type.

            - 0 = CONSTANT EFFUSION VELOCITY (VEFF)  - (H and H2 set in code for molecular diffusion/diffusion limited flux)

            - 1 = constant mixing ratio - never been used so would need testing

            - 2 = CONSTANT FLUX (SMFLUX) (option for CO2/CO/0 in code)

        Raises
        ------
        SystemExit
            When species spec is not in the model or when all(mbound != [0,1,2])
        '''
        if all(lbound!=np.array([0,1,2])):
            sys.exit('lbound must be 0, 1, or 2')
        try:
            ind = self.ispec.index(spec)
            self.photo.mbound[ind] = mbound
        except:
            sys.exit('species not in the model')

    def right_hand_side(self,usol_flat):
        '''
        Returns the right-hand-side of the system of ordinary differential equations
        defining photochemistry and transport.

        Parameters
        ----------
        usol_flat : rank-1 array with bounds (self.photo.neq)
            Atmospheric composition (mixing ratios) of all species in a single array.
            Correct order of the array is given by, for example, self.photo.usol_init.flatten('F').

        Returns
        -------
        rhs : rank-1 array with bounds (self.photo.neq)
            The right hand side of the model equations (change in mixing ratio/second)
        '''
        rhs = self.photo.right_hand_side(usol_flat)
        return rhs

    def jacobian(self,usol_flat):
        '''
        Returns the jacobian of the system of ordinary differential equations
        defining photochemistry and transport.

        Parameters
        ----------
        usol_flat : rank-1 array with bounds (self.photo.neq)
            Atmospheric composition (mixing ratios) of all species in a single array.
            Correct order of the array is given by, for example, self.photo.usol_init.flatten('F').

        Returns
        -------
        jac : rank-2 array with bounds (self.photo.neq, self.photo.lda)
            The right hand side of the model equations (change in mixing ratio/second)
        '''
        jac = self.photo.jacobian(usol_flat,self.photo.lda)
        return jac

    def out2atmosphere_txt(self,filename = 'atmosphere.txt', overwrite = False):
        '''
        Writes the atomsphere at photochemical equilbirium to a file.

        Parameters
        ----------
        filename : string, optional
            Name of the output file (default is atmosphere.txt)
        overwrite: bool, optional
            If True, then a file with the name "filename" will be overwritten
            (default is False).

        Raises
        ------
        SystemExit
            When photochemical model has not been integrated to equilibrium or
            when filename is invalid.
        '''
        if not self.code_run:
            sys.exit('Need to integrate before writting a file!')
        if os.path.isfile(filename) and not overwrite:
            sys.exit(filename+' is already a file. Choose a different name, or set overwrite = True.')
        if os.path.isdir(filename):
            sys.exit(filename+' is already a directory. Choose a different name.')

        # make dict
        out = self.out_dict()
        f = {}
        f['alt'] = out['alt']
        f['press'] = out['press']
        f['temp'] = out['T']
        f['density'] = out['den']
        f['eddy'] = self.photo.edd

        for spec in self.ispec:
            f[spec] = out[spec]

        if self.photo.co2_inert==1:
            f['CO2'] = np.ones(self.nz)*self.photo.fco2

        # particle stuff
        if self.np>0:
            particles  = ['SO4AER','S8AER','HCAER','HCAER2']
            params = ['AERSOL','WFALL','RPAR']
            for j in range(len(params)):
                for i in range(self.np):
                    if params[j] == 'AERSOL':
                        f[particles[i]+'_'+params[j]] =  self.photo.aersol[:,i]
                    if params[j] == 'WFALL':
                        f[particles[i]+'_'+params[j]] =  self.photo.wfall[:,i]
                    if params[j] == 'RPAR':
                        f[particles[i]+'_'+params[j]] =  self.photo.rpar[:,i]

        # write the file
        fil = open(filename,'w')
        for key in f.keys():
            fil.write('{:25}'.format(key))
        fil.write('\n')
        for i in range(len(f['alt'])):
            for key in f.keys():
                fil.write('{:25}'.format('%.16e'%f[key][i]))
            fil.write('\n')
        fil.close()
