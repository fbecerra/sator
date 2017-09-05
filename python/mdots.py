#!/opt/apps/python/epd/7.3.2/bin/python

from all import *

filesystem = 1
ending = 'radial'

figheight = 20.
screen_ratio = 1.

fontsize = 7.

mp.rcParams['font.size'] = fontsize

base, snapnum_start, snapnum_end, numsnaps, numsnaps_rank, rank_start = get_args(sys.argv)

time = np.array([])
mdot10 = np.array([])
mdot100 = np.array([])
mdotsink = np.zeros(138 - snapnum_start + 1)

for m in np.arange(numsnaps_rank):

   loc_start = rank_start + m

   snapnum = snapnum_start + loc_start

   path, snapnum, snapfile = setup_snapfile(filesystem, snapnum, base, ending)

   f = open(snapfile, mode = 'rb')

   num_blocks = struct.unpack('i', f.read(4))[0]
   max_string_len = struct.unpack('i', f.read(4))[0]
   string_len = str(max_string_len) + 's'
   bins = struct.unpack('i', f.read(4))[0]
   length_unit = struct.unpack('i', f.read(4))[0]
   comoving_units = struct.unpack('i', f.read(4))[0]
   xmin = struct.unpack('d', f.read(8))[0]
   xmax = struct.unpack('d', f.read(8))[0]
   xval = np.fromfile(f, dtype = 'd', count = bins)

   num_columns, num_rows = get_columns_and_rows(num_blocks, screen_ratio)

#   xlabel = get_label(length_unit, 3)

   num_blocks = 2

   for i in np.arange(num_blocks):

      yblock = f.read(max_string_len).replace("\x00", "")

      ylabel = get_label(yblock)

      yval = np.fromfile(f, dtype = 'd', count = bins)

   f.close()

   mdot10 = np.append(mdot10, yval[42])
   mdot100 = np.append(mdot100, yval[57])
   time = np.append(time, m)
  
   print 'PYTHON: Snap', loc_start, 'done!'

mdotsink = np.append(mdotsink, [7.42046039, 7.55955781, 6.86611228, 7.3292258,  .3793539, 6.57096285,
  5.87793936,  6.24945087,  5.43382213,  5.21302302,  4.50692787,  3.0195085,
  4.57795601,  4.47964136,  4.41160949,  4.28617551,  4.21918158,  3.67805776,
  4.24440728,  2.03702162,  3.98624239,  4.02798025,  4.14043194,  3.52941832])

figsize = get_figsize(screen_ratio, figheight)

fig = plt.figure(figsize = figsize)

ax = fig.add_subplot(1, 1, 1)

ax.plot(time, mdot10, 'black', label='Mdot10')
ax.plot(time, mdot100, 'blue', label='Mdot100')
ax.plot(time, mdotsink, 'red', label='MdotSink')

#ax.set_xlim([xmin, xmax])

ax.legend(loc='upper left')

ax.set_xlabel('Snapshot number', fontsize=12)

ax.set_ylabel(ylabel, fontsize=12)

fig.tight_layout()

write_image(path, base, snapnum, numsnaps, 'out')

plt.close()
