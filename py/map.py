#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

f = open("map.out")
line1 = f.readline()
headerline = map(int, line1.split())
itrap = headerline[0]
icop = headerline[1]
ictp = headerline[2]

if itrap == 1 :
    line2 = f.readline()
    trapline = map(int, line2.split())
    ipphin = trapline[0]
    neen = trapline[1]
    
    pphigridn = np.zeros((ipphin, neen), dtype=float)
    een       = np.zeros((ipphin, neen), dtype=float)
    fqcn      = np.zeros((ipphin, neen), dtype=float)

    for i1 in range(ipphin):
        for i2 in range(neen):
            line = f.readline()
            data1 = map(float, line.split())
            pphigridn[i1, i2] = data1[0]
            een[i1,i2] = data1[1]
            fqcn[i1,i2] = data1[2]
    
    plt.figure(1)
    CS = plt.contourf(pphigridn.reshape(ipphin, neen), een.reshape(ipphin, neen), fqcn.reshape(ipphin, neen),cmap=plt.cm.bone)
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel(r'$\omega_b$')
    plt.title('Trapped particle frequency map')
    plt.xlabel(r'$P_\varphi / e \Psi$')
    plt.ylabel('E (keV)')
    plt.show()

