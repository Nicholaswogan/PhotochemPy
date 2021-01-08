import subprocess
import os

cmd = 'sh valgrind_test.sh'
subprocess.run(cmd.split())

# Hadean template
fil = open('Hadean.txt','r')
linesH = fil.readlines()
fil.close()

i = 0
for line in linesH:
    message = 'ERROR SUMMARY: '
    ind = line.find('ERROR SUMMARY: ')
    if ind>0:
        errors = int(line[len(message)+ind])
        i = 1

if i!=1:
    sys.exit('did not parse valgrind log properly')

if errors != 0:
    print('Errors in Hadean template!')
else:
    print('No memory errors in Hadean template')
assert(errors==0)


# Archean template
fil = open('Archean.txt','r')
linesH = fil.readlines()
fil.close()

i = 0
for line in linesH:
    message = 'ERROR SUMMARY: '
    ind = line.find('ERROR SUMMARY: ')
    if ind>0:
        errors = int(line[len(message)+ind])
        i = 1

if i!=1:
    sys.exit('did not parse valgrind log properly')

if errors != 0:
    print('Errors in Archean template!')
else:
    print('No memory errors in Archean template')
assert(errors==0)

# Modern template
fil = open('Modern.txt','r')
linesH = fil.readlines()
fil.close()

i = 0
for line in linesH:
    message = 'ERROR SUMMARY: '
    ind = line.find('ERROR SUMMARY: ')
    if ind>0:
        errors = int(line[len(message)+ind])
        i = 1

if i!=1:
    sys.exit('did not parse valgrind log properly')

if errors != 0:
    print('Errors in Modern template!')
else:
    print('No memory errors in Modern template')
assert(errors==0)


cmd = 'rm Hadean.txt Archean.txt Modern.txt'
subprocess.run(cmd.split())
