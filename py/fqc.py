#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f = open("fqc.out")
line1 = f.readline()
ns = map(int, line1.split())
neigen = ns[0]
nr = ns[1]

r = np.zeros((nr), dtype=float)
omglocal = np.zeros((nr), dtype=float)
omgglobal = np.zeros((neigen), dtype=float)
gammaglobal = np.zeros((neigen), dtype=float)

plt.figure(1)
for i1 in range(neigen):
    line = f.readline()
    number_r = map(float, line.split())
    omgglobal[i1] = (number_r[0])
    gammaglobal[i1] = number_r[1]

for i2 in range(nr):
    line = f.readline()
    number_r = map(float, line.split())
    r[i2] = (number_r[0])
    omglocal[i2] = (number_r[1])


ll,=plt.plot(r, omglocal)
for i1 in range(neigen):
    gg, =plt.plot([0, 1], [omgglobal[i1], omgglobal[i1]])
#    leg = plt.legend([gg], ['Global Mode ' + str(i1+1)])
#    ax = plt.gca().add_artist(leg)

fs = 14

plt.legend([ll], ['Thermal GAM'], fontsize=fs)
plt.ylabel(r'$\omega$ (rad/s)', fontsize=fs)
plt.xlabel('r/a',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
# plt.title(r'$\gamma/\omega$ = ' + str(gammaglobal[0]/omgglobal[0]))
plt.tight_layout()
plt.show()
