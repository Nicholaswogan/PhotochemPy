import numpy as np
import subprocess

cmd = 'python test.py'
num_runs = 2

out = []
for i in range(num_runs):
    subprocess.run(cmd.split())
    out.append(np.loadtxt('output.txt'))

for i in range(1,num_runs):
    print(np.all(out[0]==out[i]))
    assert(np.all(out[0]==out[i]))

cmd1 = 'rm output.txt'
subprocess.run(cmd1.split())
