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
