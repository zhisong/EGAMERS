from array import array
import matplotlib.pyplot as plt

f = open("fqc.out")
line1 = f.readline()
ns = map(int, line1.split())
neigen = ns[0]
nr = ns[1]

r = array('d')
omglocal = array('d')
omgglobal = array('d')

plt.figure(1)
for i1 in range(neigen):
    line = f.readline()
    number_r = map(float, line.split())
    omgglobal.append(number_r[0])

for i2 in range(nr):
    line = f.readline()
    number_r = map(float, line.split())
    r.append(number_r[0])
    omglocal.append(number_r[1])


ll,=plt.plot(r, omglocal)
for i1 in range(neigen):
    gg, =plt.plot([0, 1], [omgglobal[i1], omgglobal[i1]])
#    leg = plt.legend([gg], ['Global Mode ' + str(i1+1)])
#    ax = plt.gca().add_artist(leg)

plt.legend([ll], ['Thermal GAM'])
plt.ylabel(r'$Re(\Omega)$')
plt.xlabel('r')
plt.show()
