import pickle
import numpy as np
import os

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
