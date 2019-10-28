import matplotlib.pyplot as plt
import numpy as np

lambdas = []
omegas = []
gammas = []
gams = []

with open('lambda0_results.txt') as f:
    Y = f.readline()
    line2 = f.readline()
    header = line2.split()
    for line in f:
        data = [float(x) for x in line.split()]
        lambdas.append(data[0])
        omegas.append(data[1])
        gammas.append(data[2])
        gams.append(data[3])

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
omegas_percent = [omegas[i]/gams[i]/2/np.pi for i in range(len(omegas))]

fs = '12'

fig1, (ax1,ax2) = plt.subplots(2)

fig1.set_size_inches(4.5,4.5)
ax1.scatter(lambdas,omegas_percent)
ax1.set_xlabel(r'$\Lambda_{0}$', fontsize=fs)
ax1.set_ylabel(r'$\omega/\omega_{GAM}$', fontsize=fs)
ax1.yaxis.set_ticks_position('both')
ax1.set_ylim([0.3,0.6])
# ax1.set_xlim([0,1.05])
setlabel(ax1,'(a)',fontsize=fs)

ax2.scatter(lambdas,growth_percent)
ax2.set_ylabel(r'$\gamma/\omega$', fontsize=fs)
ax2.set_xlabel(r'$\Lambda_{0}$', fontsize=fs)
ax2.set_ylim([0.0,1.0])
# ax2.set_xlim([0,1.05])

ax2.yaxis.set_ticks_position('both')
setlabel(ax2,'(b)',fontsize=fs)


fig1.tight_layout()
fig1.savefig('pitch_Y_'+Y.rstrip("\n")+'.png', dpi=300)

# plt.show()