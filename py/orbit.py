import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

f = open("orbit.out")
line1 = f.readline()
omgline = map(float, line1.split())
omega = omgline[0]
line2 = f.readline()
ns = map(int, line2.split())
n = ns[0]

r = np.array((0.,)*n)
theta = np.array((0.,)*n)

for i1 in range(n):
    line = f.readline()
    data = map(float, line.split())
    r[i1] = data[0]
    theta[i1] = data[1]

x = r * np.cos(theta)
y = r * np.sin(theta)

thetam = np.arange(0., 2*pi, 2*pi/n)
xm = np.cos(thetam)
ym = np.sin(thetam)

plt.figure(1)
plt.plot(x ,y ,xm, ym)
plt.axis('equal')
plt.title(r'$\omega_b$ = ' + str(omega))
plt.show()
