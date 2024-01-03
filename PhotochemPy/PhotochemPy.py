import numpy as np
from ._photochem import photochem, photochem_data, photochem_vars, photochem_wrk
import sys
import os

rootdir = os.path.dirname(os.path.realpath(__file__))+'/'
photochem_vars.rootdir = "{:500}".format(rootdir)

class PhotochemError(Exception):
    pass

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
    :ivar set_file: Name of the input settings file.
    :ivar atmosphere_txt: Name of the input atmosphere file.
    :ivar flux_txt: Name of the input solar flux file.
    :ivar code_run: If True/False then code has converged/ has not converged to equilrium.

    To import and initialize the PhotochemPy class do the following

    .. highlight:: python
    .. code-block:: python

        from PhotochemPy import PhotochemPy
        pc = PhotochemPy(species_dat, reactions_rx, set_file, atmosphere_txt, flux_txt)

    Parameters
    ----------
    species_dat : string
        Path to input file describing the species in the photochemical model,
        and their boundary conditions.
    reactions_rx : string
        Path to input file describing the reactions in the atmosphere and
        their rates.
    set_file : string
        Path to input file describing settings.
    atmosphere_txt : string
        Path to input file describing the initial atmospheric composition,
        temperature structure, eddy diffusion profile, and aersol parameters.
    flux_txt : string
        Path to input file describing the stellar flux

    '''
    def __init__(self,species_dat,reactions_rx, set_file,\
                      atmosphere_txt, flux_txt):
        self.photo = photochem
        self.data = photochem_data
        self.vars = photochem_vars
        self.wrk = photochem_wrk
        
        self.warnings = True
        
        if all(fil==None for fil in [species_dat,reactions_rx, \
                                    set_file, atmosphere_txt, flux_txt]):
            pass
        else:
            self.species_dat = species_dat
            self.reactions_rx = reactions_rx
            self.set_file = set_file
            self.atmosphere_txt = atmosphere_txt
            self.flux_txt = flux_txt

            # get species names
            fil = open(species_dat,'r')
            lines = fil.readlines()
            fil.close()
            self.ispec = [] 
            for line in lines:
                if line[0]=='*':
                    pass
                else:
                    if line.split()[1] == 'LL':
                        self.ispec.append(line.split()[0])
                    if line.split()[1] == 'SL':
                        self.ispec.append(line.split()[0])
                    if line.split()[1] == 'IN':
                        self.background_spec = line.split()[0]
                        self.ispec.append(line.split()[0])
            self.ispec.append('HV')
            self.ispec.append('M')
            err = self.photo.setup(species_dat, \
                                   reactions_rx, \
                                   set_file, \
                                   atmosphere_txt, \
                                   flux_txt)               
            if len(err.strip()) > 0:
                raise PhotochemError(err.decode("utf-8").strip())
                
            self.code_run = False

    def setup(self,species_dat,reactions_rx,set_file,\
              atmosphere_txt, flux_txt):
        '''
        In you initialize PhotochemPy with all `None` arguments, then you can run
        This to set up the atmospheres afterwords. This is necessary for some parallel
        applications (pickling errors).

        Parameters
        ----------
        species_dat : string
            Path to input file describing the species in the photochemical model,
            and their boundary conditions.
        reactions_rx : string
            Path to input file describing the reactions in the atmosphere and
            their rates.
        set_file : string
            Path to input file describing the settings.
        atmosphere_txt : string
            Path to input file describing the initial atmospheric composition,
            temperature structure, eddy diffusion profile, and aersol parameters.
        flux_txt : string
            Path to input file describing the stellar flux
        '''

        self.species_dat = species_dat
        self.reactions_rx = reactions_rx
        self.set_file = set_file
        self.atmosphere_txt = atmosphere_txt
        self.flux_txt = flux_txt

        # get species names
        fil = open(species_dat,'r')
        lines = fil.readlines()
        fil.close()
        self.ispec = [] 
        for line in lines:
            if line[0]=='*':
                pass
            else:
                if line.split()[1] == 'LL':
                    self.ispec.append(line.split()[0])
                if line.split()[1] == 'SL':
                    self.ispec.append(line.split()[0])
                if line.split()[1] == 'IN':
                    self.background_spec = line.split()[0]
                    self.ispec.append(line.split()[0])
        self.ispec.append('HV')
        self.ispec.append('M')
        err = self.photo.setup(species_dat, \
                               reactions_rx, \
                               set_file, \
                               atmosphere_txt, \
                               flux_txt)               
        if len(err.strip()) > 0:
            raise PhotochemError(err.decode("utf-8").strip())
            
        self.code_run = False
        
    def test_for_reproducibility(self):
        u0 = self.vars.usol_init.flatten(order='F').copy()
        u1 = self.vars.usol_init.flatten(order='F').copy()*2.0
        self.right_hand_side(0,u0)
        rhs1 = self.right_hand_side(0,u1)
        rhs2 = self.right_hand_side(0,u1)
        should_be_true = np.all(np.isclose(rhs1,rhs2,rtol=1.0e-8,atol=1.0e-30))
        if not should_be_true:
            raise PhotochemError("There is a problem with the right-hand-side. "+\
                                 "Two calls with the same inputs gave different results.")        

    def integrate(self,nsteps=1000,method='Backward_Euler',rtol = 1e-3, atol = 1e-27, fast_and_loose = True):
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
        if method == "CVODE_BDF":
            self.vars.max_cvode_steps = nsteps
            converged, err = self.photo.cvode_equilibrium(rtol,atol,fast_and_loose)
            if len(err.strip()) > 0:
                raise PhotochemError(err.decode("utf-8").strip())
        elif method == "Backward_Euler":
            converged, err = self.photo.integrate(nsteps)
            if len(err.strip()) > 0:
                raise PhotochemError(err.decode("utf-8").strip())

        if not converged:
            self.code_run = False
        else:
            self.code_run = True

            # check redox conservation
            if np.abs(self.vars.redox_factor) > 1e-3 and self.warnings:
                print('Warning, redox conservation is not very good.')
                print('redox factor =','%.2e'%self.vars.redox_factor)
                
            # check for mixing ratios greater than 1
            if np.max(self.vars.usol_out) > 1 and self.warnings:
                print('Warning, some mixing ratios are greater than 1.')
                
        return self.code_run

    def evolve(self,t0,usol_start,t_eval,rtol = 1.0e-3, atol= 1e-27, nsteps = 1000000, \
               fast_and_loose = True, outfile = None, overwrite = False, amount2save = 1):
        """Evolves the atmosphere with the CVODE BDF integrator from Sundials.

        Parameters
        ----------
        t0 : float
            Starting time (s)
        usol_start : Array{float,2}
            Initial conditions. nq by nz array of atmospheric mixing ratios.
        t_eval : Vector{float}
            Times to evaluate the solution (s)
        rtol : float
            Relative tolerance. Probably don't go higher than 1e-3.
        atol : float
            Absolute tolerance. About 1e-25 works well for rtol=1e-3.
            For low rtol (~1e-5) then use rtol=~1e-30.
        fast_and_loose : bool
            If 1, then will use a fast approximation to the jacobian.
            If 0, then CVODE will compute a more accurate jacobian (slowly).
        outfile : string
            If a file path is given, the the solution will be appended to the file "outfile"
            throughout the simulation. If this is used, then None is returned

        Returns
        -------
        solution : Array{float,3}
            Array of dimension [len(t_eval),nq,nz] containing mixing ratios of the atmosphere at
            each time.
        """
        if usol_start.shape != (self.data.nq, self.data.nz):
            raise PhotochemError('usol_start is the wrong shape')
            
        self.vars.max_cvode_steps = nsteps

        # in this case num_sol = len(t_eval)
        if outfile == None:
            num_sol = len(t_eval)
            solution, success, err = self.photo.cvode(t0,usol_start,t_eval,rtol,atol,fast_and_loose)
            if len(err.strip()) > 0:
                raise PhotochemError(err.decode("utf-8").strip())
            return solution
        else:
            if os.path.isfile(outfile) and not overwrite:
                raise PhotochemError(outfile,' is already a file.')
            success, err = self.photo.cvode_save(t0,usol_start,t_eval,rtol,atol,fast_and_loose,outfile,amount2save)
            if len(err.strip()) > 0:
                raise PhotochemError(err.decode("utf-8").strip())
            if not success:
                raise PhotochemError('CVODE returned an error.')
            return None

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
            raise PhotochemError('Need to integrate before outputting a solution!')
        elif self.code_run:
            out = {}
            out['alt'] = self.data.z/1e5
            out['den'] = self.vars.den
            out['press'] = self.vars.p
            out['T'] = self.vars.t
            for i in range(self.data.nq):
                out[self.ispec[i]] = self.vars.usol_out[i,:]
            for i in range(self.data.nq,self.data.nq+self.data.isl):
                out[self.ispec[i]] = self.wrk.d[i]/self.vars.den
            out[self.ispec[-3]] = self.wrk.d[-3]/self.vars.den
            out[self.ispec[-2]] = self.wrk.d[-2]/self.vars.den
            out[self.ispec[-1]] = self.wrk.d[-1]/self.vars.den
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
        out['alt'] = self.data.z/1e5
        out['den'] = self.vars.den
        out['press'] = self.vars.p
        out['T'] = self.vars.t
        for i in range(self.data.nq):
            out[self.ispec[i]] = self.vars.usol_init[i,:]
        for i in range(self.data.nq,self.data.nq+self.data.isl):
            out[self.ispec[i]] = self.wrk.d[i]/self.vars.den
        out[self.ispec[-1]] = self.wrk.d[-3]/self.vars.den
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
            raise PhotochemError('Need to integrate before outputing surface flux!')
        elif self.code_run:
            out = {}
            for i in range(self.data.nq):
                out[self.ispec[i]] = self.vars.flow[i]
            return out

    def reset(self):
        '''
        Resets the problem by reading in the original input files (e.g. species_dat, ...)
        '''
        err = self.photo.setup(self.species_dat, \
                             self.reactions_rx, \
                             self.set_file, \
                             self.atmosphere_txt, \
                             self.flux_txt)          
        if len(err.strip()) > 0:
            raise PhotochemError(err.decode("utf-8").strip())
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
            raise PhotochemError('Need to integrate before setting the output as input!')
        else:
            self.vars.usol_init = self.vars.usol_out
            if self.data.np > 0:
                self.vars.rpar_init = self.wrk.rpar
                self.vars.wfall_init = self.wrk.wfall
                self.vars.aersol_init = self.wrk.aersol

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
            if ind+1 > self.data.nq:
                raise PhotochemError('Only long-lived species can have boundary conditions')
        except ValueError:
            raise PhotochemError('species not in the model')
        if self.vars.lbound[ind] == 2:
            self.vars.sgflux[ind] = flx
        elif self.vars.lbound[ind] == 3:
            self.vars.distflux[ind] = flx
        else:
            print('mbound set to',self.vars.lbound[ind],'so the flux will not change')

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
            if ind+1 > self.data.nq:
                raise PhotochemError('Only long-lived species can have boundary conditions')
        except ValueError:
            raise PhotochemError('species not in the model')
        if self.vars.lbound[ind] == 1:
            self.vars.fixedmr[ind] = mix
        else:
            print('lbound set to',self.vars.lbound[ind],'so the mixing ratio will not change')

    def set_vdep(self,spec,vdep):
        '''
        Sets surface deposition velocity boundary condition for a species.
        self.photo.lbound must be 0 or 3.

        Parameters
        ----------
        spec : string
            Species to change surface mixing ratio
        vdep : float
            Surface deposition velocity of species spec (cm/s)

        Raises
        ------
        SystemExit
            When species spec is not in the model
        '''
        try:
            ind = self.ispec.index(spec)
            if ind+1 > self.data.nq:
                raise PhotochemError('Only long-lived species can have boundary conditions')
        except ValueError:
            raise PhotochemError('species not in the model')
        if self.vars.lbound[ind] == 0 or self.vars.lbound[ind] == 3:
            self.vars.vdep[ind] = vdep
            self.vars.vdep0[ind] = vdep
        else:
            print('lbound set to',self.vars.lbound[ind],'so the vdep will not change')

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
        if all(lbound != np.array([0,1,2,3])):
            raise PhotochemError('lbound must be 0, 1, 2 or 3')
        try:
            ind = self.ispec.index(spec)
            if ind+1 > self.data.nq:
                raise PhotochemError('Only long-lived species can have boundary conditions')
            self.vars.lbound[ind] = lbound
        except ValueError:
            raise PhotochemError('species not in the model')

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
            raise PhotochemError('lbound must be 0, 1, or 2')
        try:
            ind = self.ispec.index(spec)
            if ind+1 > self.data.nq:
                raise PhotochemError('Only long-lived species can have boundary conditions')
            self.vars.mbound[ind] = mbound
        except ValueError:
            raise PhotochemError('species not in the model')

    def right_hand_side(self,tn,usol_flat):
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
        rhs, err = self.photo.right_hand_side(tn, usol_flat)
        if len(err.strip()) > 0:
            raise PhotochemError(err.decode("utf-8").strip())
        return rhs

    def jacobian(self,tn,usol_flat):
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
        ldaa = self.data.nq+self.data.nq+1
        jac, err = self.photo.jacobian(tn,usol_flat,ldaa)
        if len(err.strip()) > 0:
            raise PhotochemError(err.decode("utf-8").strip())
        return jac
        
    def stm2photochem(self, stm, sol_stm, nz = 200, P_top = 2.1427e-08, \
                      rpar_layer = None, smallest = 1e-30, zero_out = []):
        
        sol_dry = stm.dry_end_atmos(sol_stm)

        stm_names = set(stm.gas.species_names)
        pc_names = set(self.ispec[0:self.data.nq])
        inter = list(pc_names.intersection(stm_names))
        diff = list(pc_names.difference(stm_names))
        
        usol_layer = np.ones(self.data.nq)*smallest
        for sp in inter:
            ind = self.ispec.index(sp)
            usol_layer[ind] = np.clip(sol_dry[sp],smallest,np.inf)
            
        for sp in zero_out:
            ind = self.ispec.index(sp)
            usol_layer[ind] = smallest
            
        rpar_default = [1.e-5, 1.e-5, 1.e-7, 1.e-7]
        if rpar_layer == None:
            rpar_layer = np.array([rpar_default[i] for i in range(self.data.np)])
        else:
            if len(rpar_layer) != self.data.np:
                raise PhotochemError("Input rpar_layer must have the sane length as self.data.np")
        
        P_surf = sol_dry['Psurf']
        
        err = self.photo.steam2photochem(nz,P_surf,P_top,usol_layer,rpar_layer)
        if len(err.strip()) > 0:
            raise PhotochemError(err.decode("utf-8").strip())

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
            raise PhotochemError('Need to integrate before writting a file!')
        if os.path.isfile(filename) and not overwrite:
            raise PhotochemError(filename+' is already a file. Choose a different name, or set overwrite = True.')
        if os.path.isdir(filename):
            raise PhotochemError(filename+' is a directory. Choose a different name.')

        # make dict
        out = self.out_dict()
        f = {}
        f['alt'] = out['alt']
        f['press'] = out['press']
        f['temp'] = out['T']
        f['density'] = out['den']
        f['eddy'] = self.vars.edd

        for i in range(self.data.nq):
            f[self.ispec[i]] = out[self.ispec[i]]

        # particle stuff
        if self.data.np > 0:
            particles  = ['SO4AER','S8AER','HCAER','HCAER2']
            params = ['AERSOL','WFALL','RPAR']
            for j in range(len(params)):
                for i in range(self.data.np):
                    if params[j] == 'AERSOL':
                        f[particles[i]+'_'+params[j]] =  self.wrk.aersol[:,i]
                    if params[j] == 'WFALL':
                        f[particles[i]+'_'+params[j]] =  self.wrk.wfall[:,i]
                    if params[j] == 'RPAR':
                        f[particles[i]+'_'+params[j]] =  self.wrk.rpar[:,i]

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
        
    def production_and_loss(self, spec, tn, usol_in):
        
        try:
            i = self.ispec.index(spec)
        except ValueError:
            raise PhotochemError('species not in the model')
        
        usol_in = np.asfortranarray(usol_in)
        # compute right hand side to update all data for usol_in
        self.right_hand_side(tn, usol_in.flatten(order='F'))
        
        # loss
        ll = self.data.numl[i]
        nz = self.data.nz
        nsp = self.data.nsp
        loss = np.empty((ll,nz))
        loss_react = ['' for a in range(ll)]
        kk = 0
        jj = 0
        for ii in range(ll):
            k = self.data.iloss[0,i,ii]-1
            n = self.data.jchem[0,k]-1 # reactant 1
            m = self.data.jchem[1,k]-1 # reactant 2
            loss[jj] = self.wrk.a[k]*self.wrk.d[n,:]*self.wrk.d[m,:]
            l1 = [self.data.jchem[2,k]-1,self.data.jchem[3,k]-1,self.data.jchem[4,k]-1]
            loss_react[jj] = self.ispec[n]+" + "+self.ispec[m]+" => "
            for l in l1:
                if l != -1:
                    loss_react[jj] += self.ispec[l]+" + "     
            loss_react[jj] = loss_react[jj][:-3]
            jj += 1
            
        pp = self.data.nump[i]
        production = np.empty((pp,nz))
        prod_react = ['' for a in range(pp)]
        kk = 0
        jj = 0
        for ii in range(pp):
            k = self.data.iprod[i,ii]-1
            n = self.data.jchem[0,k]-1 # reactant 1
            m = self.data.jchem[1,k]-1 # reactant 2
            production[jj] = self.wrk.a[k]*self.wrk.d[n,:]*self.wrk.d[m,:]

            p1 = [self.data.jchem[2,k]-1,self.data.jchem[3,k]-1,self.data.jchem[4,k]-1]
            prod_react[jj] = self.ispec[n]+" + "+self.ispec[m]+" => "
            for p in p1:
                if p != -1:
                    prod_react[jj] += self.ispec[p]+" + "     
            prod_react[jj] = prod_react[jj][:-3]
            jj += 1

        integ_loss = np.sum(loss,axis=1)*(self.data.dz[0])
        integ_prod = np.sum(production,axis=1)*(self.data.dz[0])
        
        # sort
        ind = np.argsort(integ_prod)
        prod_react = np.array(prod_react)[ind].tolist()[::-1]
        production = production[ind][::-1]
        integ_prod = integ_prod[ind][::-1]

        ind = np.argsort(integ_loss)
        loss_react = np.array(loss_react)[ind].tolist()[::-1]
        loss = loss[ind][::-1]
        integ_loss = integ_loss[ind][::-1]
        
        pl = {}
        pl['loss'] = loss
        pl['loss_react'] = loss_react
        pl['loss_integ'] = integ_loss
        pl['production'] = production
        pl['prod_react'] = prod_react
        pl['prod_integ'] = integ_prod
        
        return pl
        
    def koxy(self):
        flux = self.surf_flux()
        sol = self.out_dict()

        redoxstate = self.data.redoxstate[:self.data.nq]

        ind = np.where(redoxstate<0)[0]
        Fred = 0.0
        for i in ind:
            Fred += flux[self.ispec[i]]*-redoxstate[i]
            esc = self.vars.veff[i]*self.vars.den[-1]*sol[self.ispec[i]][-1]
            Fred -= esc*-redoxstate[i]

        ind = np.where(redoxstate>0)[0]
        Foxi = 0.0
        for i in ind:
            Foxi += flux[self.ispec[i]]*redoxstate[i]
            esc = self.vars.veff[i]*self.vars.den[-1]*sol[self.ispec[i]][-1]
            Foxi -= esc*redoxstate[i]
            
        for i,sp in enumerate(self.ispec[:self.data.nq]):
            if sp == 'CO' or sp == 'O':
                pass
            else:
                if self.vars.mbound[i] != 0:
                    raise PhotochemError('Upper boundary must be effusion velocity'//
                                        ' to compute koxy (CO and O are an exception)')
            
        return Foxi/Fred
        
        
