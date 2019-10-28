#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f = open("field.out")
line1 = f.readline()
ns = map(int, line1.split())
neigen = ns[0]
nr = ns[1]

r = np.zeros((nr), dtype=float)
err = np.zeros((nr), dtype=float)
eri = np.zeros((nr), dtype=float)

for i1 in range(neigen):
    for i2 in range(nr):
        line = f.readline()
        number_r = map(float, line.split())
        r[i2] = number_r[0]
        err[i2] = (number_r[1])
        eri[i2] = (number_r[2])
    era = np.sqrt(np.square(err) + np.square(eri))
    plt.figure(i1)
    rep,=plt.plot(r,err)
    imp,=plt.plot(r,eri)
    abp,=plt.plot(r,era)
    plt.ylabel('Er')
    plt.xlabel('r')
    plt.legend([rep, imp, abp],['Real (Mode ' + str(i1+1) + ')','Imag (Mode ' + str(i1+1) + ')', 'Abs (Mode ' + str(i1+1) + ')'])

    fs = 20
    fs2 = 24
    plt.figure()
    plt.plot(r,era)
    plt.ylabel(r'$E_r$ (normalized)',fontsize=fs2)
    plt.xlabel('r/a',fontsize=fs2)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    # plt.savefig('radial_structure_50.png')
    plt.tight_layout()

    plt.show()
