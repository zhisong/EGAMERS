import numpy as np
import matplotlib.pyplot as plt

x = [1.7, 2.0, 1.4, 2.2, 1.2]
y = [4.8, 3.2, 6.0, 7.0, 7.0]

deg = 4

p = np.polyfit(x,y,deg)

def f(x):
	#val = p[0]*x**5 + p[1]*x**4 + p[2]*x**3 + p[3]*x**2 + p[4]*x + p[5]
	val = p[0]
	for i in range(0,deg):
		val = val*x + p[i+1]
	#return val
	
	return val

vecfunc = np.vectorize(f)
xs = np.linspace(1.0,2.2,1000)
T = vecfunc(xs)

plt.figure()
plt.plot (xs, T, 'bo', xs, T, 'k')

print(p)
print(p/1000)

plt.show()
