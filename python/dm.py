
from all import *

comm.Barrier()

t1 = time.clock()

s = readsnap('mh2x1', -1, verbose = 1)

s.readblocks([1, 2, 3, 4, 5], ['pos', 'vel', 'id', 'mass', 'u', 'rho', 'vol', 'chem', 'gamma'])

comm.Barrier()

t2 = time.clock()

if rank == 0: print 'reading took', t2 - t1, 'secs'


#s.mass = s.mass[:s.npart_local[0]]

comm.Barrier()

t1 = time.clock()


bins = 1000

xval = s.x
yval = s.y

xmin = comm.allreduce(xval.min(), op = MPI.MIN)
xmax = comm.allreduce(xval.max(), op = MPI.MAX)

ymin = comm.allreduce(yval.min(), op = MPI.MIN)
ymax = comm.allreduce(yval.max(), op = MPI.MAX)

weights = s.mass
#weights = np.ones(s.mass.size)

hist, xbins, ybins = np.histogram2d(xval, yval, bins, [[xmin, xmax], [ymin, ymax]], weights = weights)

hist = comm.allreduce(hist, op = MPI.SUM)

idx = hist != 0

hist[idx] = np.log10(hist[idx])

hist = np.ma.masked_where(hist == 0, hist)

cmap = mp.cm.jet

cmap.set_bad('w')

figwidth = 20
figheight = 20
fontsize = 12
labelfac = 1.2

labelsize = labelfac * fontsize

mp.rcParams['font.size'] = fontsize

figsize = [figwidth * cm_to_inch, figheight * cm_to_inch]

fig = plt.figure(figsize = figsize)

ax = fig.add_subplot(111)

#ax.set_xlabel(r'${\rm log}\,n_{\rm H}\,[{\rm cm}^{-3}]$', fontsize = labelsize)
#ax.set_ylabel(r'${\rm log}\,T\,[{\rm K}]$', fontsize = labelsize)

im = ax.imshow(np.rot90(hist), cmap = cmap, extent = [xmin, xmax, ymin, ymax], aspect = 'auto', interpolation = 'none')

fig.tight_layout()

plt.savefig('out.eps')

comm.Barrier()

t2 = time.clock()

if rank == 0: print 'rest took', t2 - t1, 'secs'

#x = np.arange(10)
#y = x**2

#mp.plot(x, y)
#mp.xlabel('$\gamma$')
#mp.savefig('out.eps')




# mpi_print('x', rank, s.x.min(), s.x.max(), s.x.size)
# mpi_print('y', rank, s.y.min(), s.y.max(), s.y.size)
# mpi_print('z', rank, s.z.min(), s.z.max(), s.z.size)
# mpi_print('velx', rank, s.velx.min(), s.velx.max(), s.velx.size)
# mpi_print('vely', rank, s.vely.min(), s.vely.max(), s.vely.size)
# mpi_print('velz', rank, s.velz.min(), s.velz.max(), s.velz.size)
# mpi_print('id', rank, s.id.min(), s.id.max(), s.id.size)
# mpi_print('mass', rank, s.mass.min(), s.mass.max(), s.mass.size)

# if s.npartall[0] > 0:

#     mpi_print('u', rank, s.u.min(), s.u.max(), s.u.size)
#     mpi_print('nh', rank, s.nh.min(), s.nh.max(), s.nh.size)
#     mpi_print('vol', rank, s.vol.min(), s.vol.max(), s.vol.size)
#     mpi_print('temp', rank, s.temp.min(), s.temp.max(), s.temp.size)

#     mpi_print('gamma', rank, s.gamma.min(), s.gamma.max(), s.gamma.size)
#     mpi_print('mu', rank, s.mu.min(), s.mu.max(), s.mu.size)

#     mpi_print('h2i', rank, s.h2i.min(), s.h2i.max(), s.h2i.size)
#     mpi_print('hii', rank, s.hii.min(), s.hii.max(), s.hii.size)
#     mpi_print('dii', rank, s.dii.min(), s.dii.max(), s.dii.size)
#     mpi_print('hdi', rank, s.hdi.min(), s.hdi.max(), s.hdi.size)
#     mpi_print('heii', rank, s.heii.min(), s.heii.max(), s.heii.size)
#     mpi_print('heiii', rank, s.heiii.min(), s.heiii.max(), s.heiii.size)
