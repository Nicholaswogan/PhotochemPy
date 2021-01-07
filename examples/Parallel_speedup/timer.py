import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess

max_threads = int(os.cpu_count()/2)

cmd1 = 'rm times.txt'
subprocess.run(cmd1.split())

cmd2 = 'python test.py'

for i in range(1,max_threads+1):
    os.environ['OMP_NUM_THREADS'] = str(i)
    subprocess.run(cmd2.split())

threads,times = np.loadtxt('times.txt').T


plt.rcParams.update({'font.size': 15})
fig,ax = plt.subplots(1,1,figsize=[9,5])

ax.plot(threads,(times)/(times[0]),'C0o-')

ax.set_ylim(0,ax.get_ylim()[1])
ax.grid()
ax.set_ylabel('Time to run / Time to run on 1 thread')
ax.set_xlabel('Number of Threads Used')

ax.set_title('Parallel processing speedup in PhotochemPy',pad=15)

plt.savefig("Parallel_speedup.pdf",bbox_inches='tight')
