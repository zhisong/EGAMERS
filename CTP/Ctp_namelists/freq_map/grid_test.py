#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

f = open("map.out")
line1 = f.readline()
headerline = map(int, line1.split())
itrap = headerline[0]
icop = headerline[1]
ictp = headerline[2]

line2 = f.readline()
trapline = map(int, line2.split())
npphin = trapline[0]
neen = trapline[1]

pphigridn = []
een       = []

for line in f:
    data = map(float, line.split())
    pphigridn.append(data[0])
    een.append(data[1])

plt.figure(2)
plt.scatter(pphigridn, een)
plt.title('Trapped grid')
plt.xlabel(r'$P_\varphi / e \Psi$')
plt.ylabel('E (keV)')

plt.show()