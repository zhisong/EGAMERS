#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f = open("fqc.out")
line1 = f.readline()
ns = list(map(int, line1.split()))
neigen = ns[0]
nr = ns[1]
nregam = ns[2]

r = np.zeros((nr), dtype=float)
omglocal = np.zeros((nr), dtype=float)
omgegamlocal = np.zeros((nregam), dtype=float)
regamlocal = np.zeros((nregam), dtype=float)
gammaegamlocal = np.zeros((nregam), dtype=float)
omgglobal = np.zeros((neigen), dtype=float)
gammaglobal = np.zeros((neigen), dtype=float)

plt.figure(1)
for i1 in range(neigen):
    line = f.readline()
    number_r = list(map(float, line.split()))
    omgglobal[i1] = (number_r[0])
    gammaglobal[i1] = number_r[1]

for i2 in range(nr):
    line = f.readline()
    number_r = list(map(float, line.split()))
    r[i2] = (number_r[0])
    omglocal[i2] = (number_r[1])

for i3 in range(nregam):
    line = f.readline()
    number_r = list(map(float, line.split()))
    regamlocal[i3] = number_r[0]
    omgegamlocal[i3] = number_r[1]
    gammaegamlocal[i3] = number_r[2]

ll,=plt.plot(r, omglocal)
llegam,=plt.plot(regamlocal, omgegamlocal)
for i1 in range(neigen):
    gg, =plt.plot([0, 1], [omgglobal[i1], omgglobal[i1]])
#    leg = plt.legend([gg], ['Global Mode ' + str(i1+1)])
#    ax = plt.gca().add_artist(leg)

plt.legend([ll,llegam, gg], ['GAM continuum','EGAM continuum', 'Global EGAM'])
plt.ylabel(r'$Re(\Omega)$')
plt.xlabel('r')
plt.title(r'$\gamma/\omega$ = ' + str(gammaglobal[0]/omgglobal[0]))

plt.figure(2)
llegam,=plt.plot(regamlocal, gammaegamlocal)
for i1 in range(neigen):
    gg, =plt.plot([0, 1], [gammaglobal[i1], gammaglobal[i1]])
plt.legend([llegam, gg], ['EGAM continuum', 'Global EGAM'])
plt.ylabel(r'$\gamma=Im(\Omega)$')
plt.xlabel('r')
plt.title(r'$\gamma/\omega$ = ' + str(gammaglobal[0]/omgglobal[0]))
plt.show()