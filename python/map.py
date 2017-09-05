
from all import *

snapfile = '/ptmp/mpa/tgreif/h4r1/h4r1_109'

bins = 300
width = 1e-6

val_min_nh = 0#10**11
val_max_nh = 1e16

val_min_temp = 0#10**2.8
val_max_temp = 3e3

figwidth = 28
figheight = 28
fontsize = 12
labelfac = 1.2

labelsize = labelfac * fontsize

mp.rcParams['font.size'] = fontsize

figsize = [figwidth * cm_to_inch, figheight * cm_to_inch]

f = open(snapfile, mode = 'rb')

buf1 = struct.unpack('i', f.read(4))[0]

npart = np.array(struct.unpack(repr(nsnaptypes) + 'i', f.read(24)))
masstable = np.array(struct.unpack(repr(nsnaptypes) + 'd', f.read(48)))
time = struct.unpack('d', f.read(8))[0]
redshift = struct.unpack('d', f.read(8))[0]
flag_sfr = struct.unpack('i', f.read(4))[0]
flag_feedback = struct.unpack('i', f.read(4))[0]
npartall = np.array(struct.unpack(repr(nsnaptypes) + 'i', f.read(24)))
flag_cooling = struct.unpack('i', f.read(4))[0]
numfiles = struct.unpack('i', f.read(4))[0]
boxsize = struct.unpack('d', f.read(8))[0]
omega_m = struct.unpack('d', f.read(8))[0]
omega_lambda = struct.unpack('d', f.read(8))[0]
hubbleparam = struct.unpack('d', f.read(8))[0]

ngas = npart[0]

print npart
print masstable
print time
print redshift
print numfiles
print boxsize
print omega_m
print omega_lambda
print hubbleparam

f.seek(96, 1)

buf2 = struct.unpack('i', f.read(4))[0]

buf1 = struct.unpack('i', f.read(4))[0]

pos = np.fromfile(f, dtype = 'd', count = 3 * ngas)
pos *= 1 / (1 + redshift) / hubbleparam

buf2 = struct.unpack('i', f.read(4))[0]

buf1 = struct.unpack('i', f.read(4))[0]

vel = np.fromfile(f, dtype = 'd', count = 3 * ngas)

buf2 = struct.unpack('i', f.read(4))[0]

buf1 = struct.unpack('i', f.read(4))[0]

id = np.fromfile(f, dtype = 'i', count = ngas)

buf2 = struct.unpack('i', f.read(4))[0]

buf1 = struct.unpack('i', f.read(4))[0]

mass = np.fromfile(f, dtype = 'd', count = ngas)
mass *= unit_mass / solar_mass / hubbleparam

buf2 = struct.unpack('i', f.read(4))[0]

buf1 = struct.unpack('i', f.read(4))[0]

u = np.fromfile(f, dtype = 'd', count = ngas)
u *= unit_energy / unit_mass

buf2 = struct.unpack('i', f.read(4))[0]

buf1 = struct.unpack('i', f.read(4))[0]

rho = np.fromfile(f, dtype = 'd', count = ngas)
rho *= hubbleparam**2 * (1 + redshift)**3 * unit_mass / unit_length**3
nh = rho * hydrogen_massfrac / protonmass

buf2 = struct.unpack('i', f.read(4))[0]

buf1 = struct.unpack('i', f.read(4))[0]

vol = np.fromfile(f, dtype = 'd', count = ngas)
vol *= 1 / (1 + redshift)**3 / hubbleparam**3

hsml = (3. / 4 / np.pi)**(1. / 3) * vol**(1. / 3)

buf2 = struct.unpack('i', f.read(4))[0]

temp = u * 2 / 3 * protonmass / boltzmann


x = pos[: : 3]
y = pos[1 : : 3]
z = pos[2 : : 3]

idx = nh == max(nh)

cx = x[idx]
cy = y[idx]
cz = z[idx]

x -= cx
y -= cy
z -= cz

idx = np.logical_and(np.logical_and(np.fabs(x) < width / 2., np.fabs(y) < width / 2.), np.fabs(z) < width / 2.)

x = x[idx]
y = y[idx]
z = z[idx]

mass = mass[idx]
nh = nh[idx]
hsml = hsml[idx]
temp = temp[idx]

ngas = np.sum(idx)

print ngas


imin = np.minimum(np.maximum(((x + width / 2 - hsml) / width * bins).astype('i'), 0), bins - 1)
imax = np.minimum(np.maximum(((x + width / 2 + hsml) / width * bins).astype('i'), 0), bins - 1)

jmin = np.minimum(np.maximum(((y + width / 2 - hsml) / width * bins).astype('i'), 0), bins - 1)
jmax = np.minimum(np.maximum(((y + width / 2 + hsml) / width * bins).astype('i'), 0), bins - 1)

for i in np.arange(2):

   if i == 0:

      var = nh

      val_min = val_min_nh
      val_max = val_max_nh

      cmap = mp.cm.jet

      var_ending = 'nh'

   if i == 1:

      var = temp

      val_min = val_min_temp
      val_max = val_max_temp

      cmap = mp.cm.gist_heat

      var_ending = 'temp'

   map = np.zeros((bins, bins))
   sum = np.zeros((bins, bins))

   map_part = nh * nh * mass * var
   sum_part = nh * nh * mass

   for n in np.arange(ngas):

      map[imin[n] : imax[n] + 1, jmin[n]: jmax[n] + 1] += map_part[n]

      sum[imin[n] : imax[n] + 1, jmin[n]: jmax[n] + 1] += sum_part[n]

   img = np.zeros((bins, bins))

   idx = sum != 0

   img[idx] = np.log10(map[idx] / sum[idx])

   if val_min == 0:

      vmin = np.amin(img[img != 0])

   else:

      img[np.logical_and(img < np.log10(val_min), img != 0)] = np.log10(val_min)

      vmin = np.log10(val_min)

   if val_max == 0:

      vmax = np.amax(img[img != 0])

   else:

      img[np.logical_and(img > np.log10(val_max), img != 0)] = np.log10(val_max)

      vmax = np.log10(val_max)

   print vmin, vmax

   img = np.ma.masked_where(img == 0, img)

   cmap.set_bad('w')

   fig = plt.figure(figsize = figsize)

   ax = fig.add_subplot(111)

   #ax.set_xlabel(r'${\rm log}\,n_{\rm H}\,[{\rm cm}^{-3}]$', fontsize = labelsize)
   #ax.set_ylabel(r'${\rm log}\,T\,[{\rm K}]$', fontsize = labelsize)

   im = ax.imshow(np.rot90(img), cmap = cmap, vmin = vmin, vmax = vmax, extent = [0, width, 0, width], aspect = 'auto')#, interpolation = 'none')

   #fig.colorbar(im)

   fig.tight_layout()

   outfile_eps = var_ending + '.eps'

   plt.savefig(outfile_eps)

print 'PYTHON: Done!'
