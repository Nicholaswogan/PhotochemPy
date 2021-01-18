import numpy as np
from multiprocessing import Pool
from PhotochemPy import PhotochemPy
import sys
import time
import os

pc = PhotochemPy(None,None,None,None,None,None)

def wrapper(inpt):
    move_ind,init_ind,solutions,aersol_props,params,param_space,nsteps = inpt
    pc.photo.usol_init = solutions[init_ind]
    pc.photo.rpar_init = aersol_props[init_ind][0]
    pc.photo.wfall_init = aersol_props[init_ind][1]
    pc.photo.aersol_init = aersol_props[init_ind][2]
    for i,key in enumerate(params.keys()):
        if key == 'CO2':
            pc.photo.fco2 = 10**param_space[i][move_ind]
        else:
            pc.set_mr(key,10**param_space[i][move_ind])
    converged = pc.integrate(nsteps = nsteps)
    # converged = True
    if not converged:
        return [converged,np.nan,np.nan,np.nan,np.nan]
    if converged:
        aersol_props_converged = [pc.photo.rpar,pc.photo.wfall,pc.photo.aersol]
        return [converged,pc.photo.usol_out,aersol_props_converged,pc.out_dict(),pc.surf_flux()]

def sweep(files, params, max_processes = None, verbose=True, nsteps_list=[1000], nsteps_init = 1000,random_seed=True,seeds=None):

    species_dat,reactions_rx,planet,input_phot,atmosphere,sun = files
    pc.setup(species_dat,reactions_rx,planet,input_phot,atmosphere,sun)

    # If no max_processes given, then use os.cpu_count()
    if type(max_processes)!=int:
        sys.exit('max_processes must be an integer')
    if max_processes == None:
        max_processes = os.cpu_count()

    # check that all params are mixing ratios
    for key in params.keys():
        if key == 'CO2':
            pass
        else:
            ind = pc.ispec.index(key)
            if pc.photo.lbound[ind] != 1:
                sys.exit('All species must have fixed mixing ratios (lbound = 1)\n'+
                          key+' has lbound = '+str(pc.photo.lbound[ind]))

    # Make a grid for the parameters
    param_space = np.meshgrid(*(v for _, v in params.items()),indexing='ij')
    progress = np.zeros(param_space[0].shape,dtype=int)
    solutions = np.array(np.zeros(param_space[0].shape,dtype=int).tolist(),dtype=object)
    solutions_dict = np.array(np.zeros(param_space[0].shape,dtype=int).tolist(),dtype=object)*np.nan
    fluxes = np.array(np.zeros(param_space[0].shape,dtype=int).tolist(),dtype=object)*np.nan
    aersol_props = np.array(np.zeros(param_space[0].shape,dtype=int).tolist(),dtype=object)

    # find the "closest" part of the grid to input atmosphere
    closest = []
    for key in params.keys():
        if key == 'CO2':
            mr = np.log10(pc.photo.fco2)
        else:
            ind = pc.ispec.index(key) # index
            mr = np.log10(pc.photo.fixedmr[ind])
        closest.append(np.abs(params[key]-mr).argmin())
    closest = tuple(closest)
    start_vals = [params1[closest] for params1 in param_space]
    # now input the closest part of the grid
    for i,key in enumerate(params.keys()):
        if key == 'CO2':
            pc.photo.fco2 = 10**start_vals[i]
        else:
            pc.set_mr(key,10**start_vals[i])

    # find equlibrium at first grid point
    print('Trying to move to closest point on the grid...')
    converged = pc.integrate(nsteps = nsteps_init)
    if not converged:
        sys.exit("Didn't converge to first grid point.\n"\
        +"Try increasing nsteps_init, or starting with an atmosphere.txt closer to the grid.")
    print('Successful.')
    pc.out2in()
    pc.photo.verbose = False

    # save first grid point
    progress[closest] = 1
    solutions[closest] = pc.photo.usol_out
    solutions_dict[closest] = pc.out_dict()
    fluxes[closest] = pc.surf_flux()
    aersol_props[closest] = [pc.photo.rpar,pc.photo.wfall,pc.photo.aersol]


    if random_seed:
        if seeds==None:
            seeds = max_processes
        if seeds>len(np.where(progress==-0)[0]):
            seeds = len(np.where(progress==0)[0])
        print('Randomly seeding space with',seeds,'seeds')
        indexes = []
        for j in range(seeds):
            duplicates = True
            while duplicates:
                index  = []
                for i in range(len(param_space[0].shape)):
                    ind = np.random.randint(param_space[0].shape[i])
                    index.append(ind)
                index = tuple(index)
                duplicates = any([index==ind for ind in indexes])
                if len(indexes)==0:
                    duplicates=False
            indexes.append(index)

        inpts = [[indexes[i],closest,solutions,aersol_props,params,param_space,nsteps_init] for i in range(seeds)]
        with Pool(max_processes) as p:
            outs = p.map(wrapper,inpts)

        for i,out in enumerate(outs):
            if not out[0]: # if no convergence
                progress[indexes[i]] = 0
            elif out[0]: # if convergence
                progress[indexes[i]] = 1
                solutions[indexes[i]] = out[1]
                aersol_props[indexes[i]] = out[2]
                solutions_dict[indexes[i]] = out[3]
                fluxes[indexes[i]] = out[4]

        print('Done seeding.')
    # now start the algorithm
    track_progress = [progress.copy()]
    Epoch = 0
    start_loop = time.time()
    for kk,nsteps in enumerate(nsteps_list):
        progress[progress==2] = 0
        if verbose:
            print()
            print('Solving for grid with nsteps =',nsteps)
            print(' ',end='')
            [print('=',end='') for i in range(38)]
            print()
        while 1:
            # find the possible moves
            possible_moves, closest_ones = find_neighboring_zeros(progress)
            inpts = [[possible_moves[i],closest_ones[i],solutions,aersol_props,params,param_space,nsteps] for i in range(len(possible_moves))]
            processes = np.min([len(possible_moves),max_processes])

            if verbose:
                print('| ITERATION','{:<25}'.format(str(Epoch)),' |\n'\
                      '{:<26}'.format('| Converged ='),'{:<10}'.format(str(len(np.where(progress==1)[0]))),' |\n'\
                      '{:<26}'.format('| Failed to Converged ='),'{:<10}'.format(str(len(np.where(progress==2)[0]))),' |\n'\
                      '{:<26}'.format('| Unattempted ='),'{:<10}'.format(str(len(np.where(progress==0)[0]))),' |\n'\
                      '{:<26}'.format('| Possible Moves ='),'{:<10}'.format(str(len(possible_moves))),' |')
            start = time.time()
            # do all possible moves in parallel
            with Pool(processes) as p:
                outs = p.map(wrapper,inpts)
            end = time.time()

            if verbose:
                print('{:<26}'.format('| Time to try all moves ='),'{:<10}'.format(('%.0f'%(end-start))+' s'),' |')
                print('|','{:<36}'.format(''),'|')

            # now we collect our wits!
            for i,out in enumerate(outs):
                if not out[0]: # if no convergence
                    progress[possible_moves[i]] = 2
                    solutions[possible_moves[i]] = out[1]
                    aersol_props[possible_moves[i]] = out[2]
                    solutions_dict[possible_moves[i]] = out[3]
                    fluxes[possible_moves[i]] = out[4]
                elif out[0]: # if convergence
                    progress[possible_moves[i]] = 1
                    solutions[possible_moves[i]] = out[1]
                    aersol_props[possible_moves[i]] = out[2]
                    solutions_dict[possible_moves[i]] = out[3]
                    fluxes[possible_moves[i]] = out[4]

            Epoch+=1
            possible_moves, _ = find_neighboring_zeros(progress)
            track_progress.append(progress.copy())

            if len(possible_moves)==0:
                if verbose:
                    print('|','{:<36}'.format('NO MORE POSSIBLE MOVES'),'|')
                    print(' ',end='')
                    [print('=',end='') for i in range(38)]
                    print()
                break
    end_loop = time.time()
    if verbose:
        print()
        print('Exiting the program')
        print('Time of all Iterations =','%.0f'%(end_loop-start_loop),'seconds')
    return [solutions_dict,fluxes,track_progress]


def find_neighboring_zeros(array):
    # first find the ones
    ones = np.array(np.where(array==1)).T
    # now find the zeros
    possible_moves = []
    closest_ones = []
    for ind in ones: # consider every one
        for dim in range(array.ndim): # loop over dimensions
            move = np.zeros(array.ndim,dtype=int)
            for i in [-1,1]: # look in either direction
                # the potential move
                move[dim] = i
                check = move+ind
                if array.shape[dim]-1 < check[dim] or 0 > check[dim]: # if move is beyond array bounds, then pass
                    break
                # else, check the value
                val = array[tuple(check)]
                if val == 0: # if val = 0 then save its location
                    # save its location
                    possible_moves.append(check)
                    # save location of where you checked from
                    closest_ones.append(ind)
    # now remove duplicates
    temp = [tuple(x) for x in possible_moves]
    uniques = [temp.index(x) for x in set(temp)]
    possible_moves = [tuple(possible_moves[u]) for u in uniques]
    closest_ones = [tuple(closest_ones[u]) for u in uniques]
    return possible_moves,closest_ones
