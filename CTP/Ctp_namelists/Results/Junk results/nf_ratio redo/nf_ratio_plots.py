import matplotlib.pyplot as plt
import numpy as np

nfs = []
omegas = []
gammas = []


with open('nf_ratio_results') as f:

    line1 = f.readline()
    header = line1.split()
    for line in f:
        data = [float(x) for x in line.split()]
        nfs.append(data[0])
        omegas.append(data[1])
        gammas.append(data[2])

growth_percent = [gammas[i]/omegas[i] for i in range(len(omegas))]

plt.subplot(2,1,1)
plt.gca().set_title('Omega')
plt.scatter(nfs,omegas)
# plt.xlim(0,0.2)
# plt.subplot(3,1,2)
# plt.gca().set_title('Gamma')
# plt.scatter(nfs,gammas)
# plt.xlim(0,0.2)
plt.subplot(2,1,2)
plt.gca().set_title('gamma/omega')
plt.scatter(nfs,growth_percent)

plt.show()


