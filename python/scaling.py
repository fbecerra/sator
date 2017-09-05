
from all import *

figheight = 10.
screen_ratio = 1.
fontsize = 9.

mp.rcParams['font.size'] = fontsize
mp.rc('legend',**{'fontsize': fontsize})

figsize = get_figsize(screen_ratio, figheight)

fig = plt.figure(figsize = figsize)

x1 = np.array([16., 32., 64., 128., 256., 512.])
y1 = np.array([53619., 31333., 16029., 9072., 5701., 4882.], dtype = 'd')

x2 = np.array([1, 2, 4, 8, 16])
y2 = np.array([1940, 1022, 576, 396, 480], dtype = 'd')

ideal1 = x1 / np.min(x1)
ideal2 = x2 / np.min(x2)

y1 = 1 / y1
y1 /= np.min(y1)

y2 = 1 / y2
y2 /= np.min(y2)

ax1 = fig.add_subplot(111)

ax1.plot(x1, ideal1, 'red', label = r'${\rm Ideal}$')
ax1.plot(x1, y1, 'black', label = r'${\rm Simulation}$')

ax1.set_xlim(12, 600)
ax1.set_ylim(0.8, 40)

ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.set_xlabel(r'${\rm Number\,of\,cores}$')
ax1.set_ylabel(r'${\rm Scaling}$')

ax1.set_title(r'${\rm MPI\,Scaling\,Test}$')

ax1.legend(loc = 2)

#ax2 = fig.add_subplot(122)

#ax2.plot(x2, ideal2, 'red')
#ax2.plot(x2, y2, 'black')

#ax2.set_xlim(0.8, 20)
#ax2.set_ylim(0.8, 20)

#ax2.set_xscale('log')
#ax2.set_yscale('log')

#ax2.set_xlabel(r'${\rm Number\,of\,threads}$')
#ax2.set_ylabel(r'${\rm Scaling}$')

#ax2.set_title(r'${\rm Tree\,walk\,-\,OpenMP}$')

fig.tight_layout()

plt.savefig('fig_scaling.pdf')
