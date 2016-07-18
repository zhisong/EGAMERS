from array import array
import matplotlib.pyplot as plt

f = open("field.out")
line1 = f.readline()
ns = map(int, line1.split())
neigen = ns[0]
nr = ns[1]

r = array('d')
err = array('d')
eri = array('d')

for i1 in range(neigen):
    for i2 in range(nr):
        line = f.readline()
        number_r = map(float, line.split())
        r.append(number_r[0])
        err.append(number_r[1])
        eri.append(number_r[2])
    plt.figure(i1)
    rep,=plt.plot(r,err)
    imp,=plt.plot(r,eri)
    plt.ylabel('Er')
    plt.xlabel('r')
    plt.legend([rep, imp],['Real (Mode ' + str(i1+1) + ')','Imag (Mode ' + str(i1+1) + ')'])
    plt.show()
