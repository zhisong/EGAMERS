#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

f = open("map.out")
line1 = f.readline()
headerline = map(int, line1.split())
itrap = headerline[0]
icop = headerline[1]
ictp = headerline[2]

if (ictp == 1):
    line2 = f.readline()
    ctpline = map(int, line2.split())
    ipphin = ctpline[0]
    neen = ctpline[1]
    
    pphigridn = np.zeros((ipphin, neen), dtype=float)
    een       = np.zeros((ipphin, neen), dtype=float)
    fqcn      = np.zeros((ipphin, neen), dtype=float)
    dfdee     = np.zeros((ipphin, neen), dtype=float)
    gam       = np.zeros((ipphin, neen), dtype=float)

    for i1 in range(ipphin):
        for i2 in range(neen):
            line = f.readline()
            data1 = map(float, line.split())
            pphigridn[i1, i2] = data1[0]
            een[i1,i2] = data1[1]
            fqcn[i1,i2] = data1[2]        
            dfdee[i1,i2] = data1[3]
            # gam[i1,i2] = data1[2] - 213266.51149959548
    
    fs = '14'

    plt.figure(figsize=(7,5))

    CS = plt.contourf(pphigridn.reshape(ipphin, neen), een.reshape(ipphin, neen), fqcn.reshape(ipphin, neen),cmap=plt.cm.bone)
    cbar = plt.colorbar(CS, format=OOMFormatter(3, mathText=True))
    cbar.ax.set_ylabel(r'$\omega_b$ (rad/s) ', fontsize=fs)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.update_ticks()
    plt.xlabel(r'$P_\varphi / e \Psi$', fontsize=fs)
    plt.ylabel('E (keV)', fontsize=fs)

    plt.savefig('ctp_fqc_map.png', bbox_inches='tight', dpi=300)


    plt.figure(figsize=(7,5))

    CS = plt.contourf(pphigridn.reshape(ipphin, neen), een.reshape(ipphin, neen), dfdee.reshape(ipphin, neen),cmap=plt.cm.bone)
    cbar = plt.colorbar(CS,format=OOMFormatter(-4, mathText=True))
    cbar.ax.set_ylabel(r'$\frac{\partial f_0}{\partial E}$', fontsize='18')
    cbar.formatter.set_powerlimits((0, 0))
    cbar.update_ticks()
    plt.xlabel(r'$P_\varphi / e \Psi$',fontsize=fs)
    plt.ylabel('E (keV)', fontsize=fs)

    plt.savefig('ctp_grid_weight.png', bbox_inches='tight',dpi=300)

    # plt.figure(3)
    # CS = plt.contourf(pphigridn.reshape(ipphin, neen), een.reshape(ipphin, neen), gam.reshape(ipphin, neen),cmap=plt.cm.bone)
    # cbar = plt.colorbar(CS)
    # cbar.ax.set_ylabel(r'$\omega_b - \omega$')
    # plt.title('Ctp gam map')
    # plt.xlabel(r'$P_\varphi / e \Psi$')
    # plt.ylabel('E (keV)')

# plt.show()