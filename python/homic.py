
from all import *
from com import *

nside = 128
boxsize = 1
nh = 1e-2
temp = 1e3
mfac = 1e-1

file = '/ptmp/mpa/tgreif/iso1/iso1_ic'

npart = nside**3
rhocloud = nh * protonmass / hydrogen_massfrac / unit_mass * unit_length**3
mcloud = 4 * pi * nh * protonmass / hydrogen_massfrac * (boxsize / 2 * unit_length)**3 / unit_mass
mu = 4 / (1 + 3 * hydrogen_massfrac)
utherm = temp / (gamma_adb - 1) / mu * boltzmann / protonmass * unit_mass / unit_energy

maxrad = np.sqrt(3) / 2 * boxsize
radfac = 1e-2


id = np.arange(npart, dtype = 'i')

x = (id / nside**2 + 0.5) * boxsize / nside - boxsize / 2
y = (id / nside % nside + 0.5) * boxsize / nside - boxsize / 2
z = (id % nside + 0.5) * boxsize / nside - boxsize / 2

id += 1

r = np.sqrt(x**2 + y**2 + z**2)

theta = np.zeros(npart, dtype = 'd')

idx = r != 0

theta[idx] = np.arccos(z[idx] / r[idx])

phi = np.zeros(npart, dtype = 'd')

idx = x > 0

phi[idx] = np.arctan(y[idx] / x[idx])

idx = x == 0

phi[idx] = np.sign(y[idx]) * pi / 2

idx = (x < 0) & (y >= 0)

phi[idx] = np.arctan(y[idx] / x[idx]) + pi

idx = (x < 0) & (y < 0)

phi[idx] = np.arctan(y[idx] / x[idx]) - pi

rcloud = boxsize / 4

idx = r < rcloud

r[idx] = (r[idx] / rcloud)**2

x[idx] = r[idx] * np.sin(theta[idx]) * np.cos(phi[idx])
y[idx] = r[idx] * np.sin(theta[idx]) * np.sin(phi[idx])
z[idx] = r[idx] * np.cos(theta[idx])

x += boxsize / 2
y += boxsize / 2
z += boxsize / 2

# r = np.zeros(npart, dtype = 'd')

# for i in np.arange(npart):

#     rtmp = 0.

#     while rtmp < radfac * maxrad:

#         rtmp = np.random.random_sample()**(2) * maxrad

#     r[i] = rtmp

# theta = np.arccos(2 * np.random.random_sample(npart) - 1)
# phi = 2 * pi * np.random.random_sample(npart)

# x = boxsize / 2 + r * np.sin(theta) * np.cos(phi)
# y = boxsize / 2 + r * np.sin(theta) * np.sin(phi)
# z = boxsize / 2 + r * np.cos(theta)

# idx = (x >= 0) & (x < boxsize) & (y >=0) & (y < boxsize) & (z >= 0) & (z < boxsize)

# npart = np.sum(idx)

# x = x[idx]
# y = y[idx]
# z = z[idx]

pos = np.zeros([npart, 3], dtype = 'd')

pos[:, ::3] = x.reshape(npart, 1)
pos[:, 1::3] = y.reshape(npart, 1)
pos[:, 2::3] = z.reshape(npart, 1)

# rad = np.sqrt((x - boxsize / 2)**2 + (y - boxsize / 2)**2 + (z - boxsize / 2)**2)

# idx = rad >= boxsize / 2

vel = np.zeros([npart, 3], dtype = 'd')



mass = np.zeros(npart, dtype = 'd')
mass.fill(mcloud / npart)

#mass[idx] = mfac * mcloud / npart

u = np.zeros(npart, dtype = 'd')
u.fill(utherm)

f = open(file, mode = 'wb')

nfiles = 1
hubbleparam = 1.0
blksize = 256

f.write(struct.pack('i', blksize))
f.write(struct.pack('i', npart))
np.zeros(23, dtype = 'i').tofile(f)
f.write(struct.pack('i', npart))
np.zeros(6, dtype = 'i').tofile(f)
f.write(struct.pack('i', nfiles))
f.write(struct.pack('d', boxsize))
np.zeros(4, dtype = 'i').tofile(f)
f.write(struct.pack('d', hubbleparam))
np.zeros(24, dtype = 'i').tofile(f)
f.write(struct.pack('i', blksize))

blksize = 3 * 8 * npart

f.write(struct.pack('i', blksize))
pos.tofile(f)
f.write(struct.pack('i', blksize))

f.write(struct.pack('i', blksize))
vel.tofile(f)
f.write(struct.pack('i', blksize))

blksize = 1 * 4 * npart

f.write(struct.pack('i', blksize))
id.tofile(f)
f.write(struct.pack('i', blksize))

blksize = 1 * 8 * npart

f.write(struct.pack('i', blksize))
mass.tofile(f)
f.write(struct.pack('i', blksize))

f.write(struct.pack('i', blksize))
u.tofile(f)
f.write(struct.pack('i', blksize))

f.close()
