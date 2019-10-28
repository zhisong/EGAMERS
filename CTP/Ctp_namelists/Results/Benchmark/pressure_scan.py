import matplotlib.pyplot as plt
import numpy as np

ys = []
omegas = []
gammas = []

with open('results1_8.txt') as f:
    Z = f.readline()
    gam = float(f.readline())*2*np.pi
    line2 = f.readline()
    header = line2.split()
    for line in f:
        data = [float(x) for x in line.split()]
        ys.append(data[0]/2./2.11) # y values too big by 2*
        omegas.append(data[1])
        gammas.append(data[2])

def setlabel(ax, label, loc=2, borderpad=0.6, **kwargs):
    legend = ax.get_legend()
    if legend:
        ax.add_artist(legend)
    line, = ax.plot(np.NaN,np.NaN,color='none',label=label)
    label_legend = ax.legend(handles=[line],loc=loc,handlelength=0,handleheight=0,handletextpad=0,borderaxespad=0,borderpad=borderpad,frameon=False,**kwargs)
    label_legend.remove()
    ax.add_artist(label_legend)
    line.remove()

growth_percent = [gammas[i]/omegas[i] for i in range(len(omegas))]
omegas_percent = [o/gam for o in omegas]

fs = '14'

fig1, (ax1,ax2) = plt.subplots(2)

fig1.set_size_inches(4.5,4.5)
ax1.scatter(ys,omegas_percent)
ax1.set_xlabel('Y', fontsize=fs)
ax1.set_ylabel(r'$\omega/\omega_{GAM}$', fontsize=fs)
ax1.yaxis.set_ticks_position('both')
ax1.set_ylim([0,1.2])
ax1.set_yticks([i*0.3 for i in range(5)])
ax1.set_xlim([0,8])
# ax1.set_xticks([i for i in range(7)])
ax1.tick_params(labelsize=fs)
setlabel(ax1,'(a)',fontsize=fs)

ax2.scatter(ys,growth_percent)
ax2.set_ylabel(r'$\gamma/\omega$', fontsize=fs)
ax2.set_xlabel('Y', fontsize=fs)
# ax2.set_ylim([0.0,0.12])
# ax2.set_yticks([i*0.04 for i in range(4)])
ax2.set_xlim([0,8])
# ax2.set_xticks([i*0.04 for i in range(4)])
ax2.tick_params(labelsize=fs)

ax2.yaxis.set_ticks_position('both')
setlabel(ax2,'(b)',fontsize=fs)


fig1.tight_layout()
# fig1.savefig('bench1_8_zoom.png', dpi=200)

plt.show()