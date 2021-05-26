import pickle
import numpy as np
from scipy.io import FortranFile, FortranEOFError
import os

# now outdated and I will remove eventually
def save2pickle(filename,outname,delete = False):
    fil = open(filename,'r')
    lines = fil.readlines()
    fil.close()

    specs = lines[0].split()

    inds = []
    for i,line in enumerate(lines):
        if line[0]=='t':
            inds.append(i-1)
    inds.append(len(lines))

    nz = inds[1]-inds[0]-2
    nt = len(inds)-1

    times = np.zeros(nt)
    solution = np.zeros([nt,len(specs),nz])

    for i in range(len(inds)-1):
        t = float(lines[inds[i]+1].split()[2])
        times[i] = t
        n = int(lines[inds[i]+1].split()[5])
        for j,line in enumerate(lines[inds[i]+2:inds[i+1]]):
            solution[i,:,j] = [float(a) for a in line.split()]

    out_dict = {}
    out_dict['time'] = times
    out_dict['solution'] = solution
    out_dict['species'] = specs

    fil = open(outname,'wb')
    pickle.dump(out_dict,fil)
    fil.close()

    if delete == True:
        os.remove(filename)

def read_evolve_output(filename):
    '''reads fortran binary file exported by cvode.f90/cvode_save
    '''
    f = FortranFile(filename, 'r')
    nq = f.read_record(np.int32)[0]
    npp = f.read_record(np.int32)[0]
    nw = f.read_record(np.int32)[0]
    kj = f.read_record(np.int32)[0]
    ks = f.read_record(np.int32)[0]
    nz = f.read_record(np.int32)[0]
    tmp = f.read_record(np.dtype("S8"))
    background_gas = np.char.strip(tmp.astype('U'))[0]
    tmp = f.read_record(np.dtype("S8"))
    ispec = np.char.strip(tmp.astype('U')).tolist()
    photospec = f.read_record(np.int32)
    
    # new
    nsp = f.read_record(np.int32)[0]
    nr = f.read_record(np.int32)[0]
    nmax = f.read_record(np.int32)[0]
    jchem = f.read_record(np.int32).reshape((5,nr),order='F')
    iprod = f.read_record(np.int32).reshape((nsp,nmax),order='F')
    iloss = f.read_record(np.int32).reshape((2,nsp,nmax),order='F')
    nump = f.read_record(np.int32)
    numl = f.read_record(np.int32)
    A = f.read_record(np.double).reshape((nr,nz),order='F')
    # new
    
    z = f.read_record(np.double)
    T = f.read_record(np.double)
    edd = f.read_record(np.double)
    rpar = f.read_record(np.double).reshape((nz,npp),order='F')
    wavl = f.read_record(np.double)
    flux = f.read_record(np.double)

    nt = f.read_record(np.int32)[0]
    amount2save = f.read_record(np.int32)[0]

    photospecies = [ispec[j-1] for j in photospec]
    
    sol = {}
    sol['background_gas'] = background_gas
    sol['ispec'] = ispec
    sol['photospecies'] = photospecies
    sol['alt'] = z
    sol['T'] = T
    sol['edd'] = edd
    sol['rpar'] = rpar
    sol['wavl'] = wavl
    sol['flux'] = flux
    
    sol['jchem'] = jchem
    sol['iprod'] = iprod
    sol['iloss'] = iloss
    sol['nump'] = nump
    sol['numl'] = numl
    sol['A'] = A
    
    for sp in ispec:
        sol[sp] = np.empty((nt,nz))
    sol['time'] = np.empty((nt))
    sol['den'] = np.empty((nt,nz))
    if amount2save == 1:
        for sp in ispec[:nq]:
            sol[sp+'_chemprod'] = np.empty((nt,nz))
            sol[sp+'_chemloss'] = np.empty((nt,nz))
        sol['P'] = np.empty((nt,nz))
        sol['surf_radiance'] = np.empty((nt,nw))
        for sp in photospecies:
            sol[sp+"_prates"] = np.empty((nt,nz))
    
    for j in range(nt):
        try:
            tmp = f.read_record(np.int32)
            t = f.read_record(np.double)[0]
            D = f.read_record(np.double).reshape((nsp+2,nz),order='F')
            den = f.read_record(np.double)
            if amount2save == 1:
                P = f.read_record(np.double)
                surf_radiance = f.read_record(np.double)
                photorates = f.read_record(np.double).reshape((ks,nz),order='F')
                yp = f.read_record(np.double).reshape((nq,nz),order='F')
                loss = f.read_record(np.double).reshape((nq,nz),order='F')

            # save
            for i,sp in enumerate(ispec):
                sol[sp][j] = D[i]/den
            sol['time'][j] = t
            sol['den'][j] = den
            if amount2save == 1:
                for i,sp in enumerate(ispec[:nq]):
                    sol[sp+"_chemprod"][j] = yp[i]
                    sol[sp+"_chemloss"][j] = loss[i]
                sol['P'][j] = P
                sol['surf_radiance'][j] = surf_radiance
                for i,sp in enumerate(photospecies):
                    sol[sp+'_prates'][j] = photorates[i,:]  
        except FortranEOFError:
            break

    if j != nt-1:        
        # need to delete j:
        for i,sp in enumerate(ispec):
            sol[sp] = np.delete(sol[sp],slice(j,nt),axis=0)
        sol['time'] = np.delete(sol['time'],slice(j,nt),axis=0)
        sol['den'] = np.delete(sol['den'],slice(j,nt),axis=0)
        if amount2save == 1:
            for i,sp in enumerate(ispec[:nq]):
                sol[sp+"_chemprod"] = np.delete(sol[sp+"_chemprod"],slice(j,nt),axis=0)
                sol[sp+"_chemloss"] = np.delete(sol[sp+"_chemloss"],slice(j,nt),axis=0)
            sol['P'] = np.delete(sol['P'],slice(j,nt),axis=0)
            sol['surf_radiance'] = np.delete(sol['surf_radiance'],slice(j,nt),axis=0)
            for i,sp in enumerate(photospecies):
                sol[sp+'_prates'] = np.delete(sol[sp+'_prates'],slice(j,nt),axis=0)
    
    f.close()
    
    return sol

def prod_and_loss(sol, spec, j):
    '''Compute production and loss reactions (vs z) for a given species.
    Must input sol from read_evolve_output
    '''
    i = sol['ispec'].index(spec)
    # loss
    ll = sol['numl'][i]
    nz = len(sol['T'])
    nsp = sol['iprod'].shape[0]
    loss = np.empty((ll,nz))
    loss_react = ['' for a in range(ll)]
    kk = 0
    jj = 0
    for ii in range(ll):
        k = sol['iloss'][0,i,ii]-1
        n = sol['jchem'][0,k]-1 # reactant 1
        m = sol['jchem'][1,k]-1 # reactant 2
        # try:
        # if sol['ispec'][m] == "HV":
        #     loss[jj] = sol['A'][k]*sol[sol['ispec'][n]][j]*sol['den'][j]
        #     l1 = [sol['jchem'][2,k]-1,sol['jchem'][3,k]-1,sol['jchem'][4,k]-1]
        #     loss_react[jj] = sol['ispec'][n]+" + "+'hv'+" => "
        #     for l in l1:
        #         if l != -1:
        #             loss_react[jj] += sol['ispec'][l]+" + "     
        #     loss_react[jj] = loss_react[jj][:-3]
        # else:
        loss[jj] = sol['A'][k]*sol[sol['ispec'][n]][j]*sol[sol['ispec'][m]][j]*sol['den'][j]**2
        l1 = [sol['jchem'][2,k]-1,sol['jchem'][3,k]-1,sol['jchem'][4,k]-1]
        loss_react[jj] = sol['ispec'][n]+" + "+sol['ispec'][m]+" => "
        for l in l1:
            if l != -1:
                loss_react[jj] += sol['ispec'][l]+" + "     
        loss_react[jj] = loss_react[jj][:-3]
        jj += 1
        # except:
        #     kk += 1
        #     pass
    kk = 0
    loss = np.delete(loss,slice(ll-kk,ll),axis=0)
    loss_react = np.delete(loss_react,slice(ll-kk,ll)).tolist()

    # production
    pp = sol['nump'][i]
    nz = len(sol['T'])
    production = np.empty((pp,nz))
    prod_react = ['' for a in range(pp)]
    kk = 0
    jj = 0
    for ii in range(pp):
        k = sol['iprod'][i,ii]-1
        n = sol['jchem'][0,k]-1 # reactant 1
        m = sol['jchem'][1,k]-1 # reactant 2
        # try:
        production[jj] = sol['A'][k]*sol[sol['ispec'][n]][j]*sol[sol['ispec'][m]][j]*sol['den'][j]**2

        p1 = [sol['jchem'][2,k]-1,sol['jchem'][3,k]-1,sol['jchem'][4,k]-1]
        prod_react[jj] = sol['ispec'][n]+" + "+sol['ispec'][m]+" => "
        for p in p1:
            if p != -1:
                prod_react[jj] += sol['ispec'][p]+" + "     
        prod_react[jj] = prod_react[jj][:-3]
        jj += 1
        # except:
        # kk += 1
        # pass
    kk = 0
    production = np.delete(production,slice(pp-kk,pp),axis=0)
    prod_react = np.delete(prod_react,slice(pp-kk,pp)).tolist()

    integ_loss = np.sum(loss,axis=1)*(sol['alt'][1]-sol['alt'][0])
    integ_prod = np.sum(production,axis=1)*(sol['alt'][1]-sol['alt'][0])

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
            

        
def stm2atmosphere_txt(pc, stm, sol_stm, filename = 'impact_atmosphere.txt', overwrite=False, smallest = 1e-30):
    
    if os.path.isfile(filename) and not overwrite:
        raise Exception(filename+' is already a file. Choose a different name, or set overwrite = True.')
    if os.path.isdir(filename):
        raise Exception(filename+' is a directory. Choose a different name.')
    
    sol_dry = stm.dry_end_atmos(sol_stm)
    
    stm_names = set(stm.gas.species_names)
    pc_names = set(pc.ispec[0:pc.data.nq])
    inter = list(pc_names.intersection(stm_names))
    diff = list(pc_names.difference(stm_names))

    nz = pc.data.nz
    f = {}
    out = pc.in_dict()
    f['alt'] = out['alt']
    f['temp'] = out['T']
    f['eddy'] = pc.vars.edd

    for sp in inter:
        f[sp] = np.clip(sol_dry[sp],smallest,np.inf)*np.ones(nz)
    for sp in diff:
        f[sp] = np.ones(nz)*smallest

    # particle stuff
    if pc.data.np > 0:
        particles  = ['SO4AER','S8AER','HCAER','HCAER2']
        params = ['AERSOL','WFALL','RPAR']
        for j in range(len(params)):
            for i in range(pc.data.np):
                if params[j] == 'AERSOL':
                    f[particles[i]+'_'+params[j]] =  pc.wrk.aersol[:,i]
                if params[j] == 'WFALL':
                    f[particles[i]+'_'+params[j]] =  pc.wrk.wfall[:,i]
                if params[j] == 'RPAR':
                    f[particles[i]+'_'+params[j]] =  pc.wrk.rpar[:,i]

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
    
    print('Surface Pressure (bars) =', '%.4f'%sol_dry['Psurf'])
