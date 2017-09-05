#!/opt/apps/python/epd/7.3.2/bin/python

from all import *

filesystem = 1
ending = 'radial'

figheight = 20.
screen_ratio = 1.

fontsize = 7.

mp.rcParams['font.size'] = fontsize

base, snapnum_start, snapnum_end, numsnaps, numsnaps_rank, rank_start = get_args(sys.argv)

#print 'hola1'
#print 'Size: ', size, ' and Rank:', rank, '\n'
#every_snapshot = 15
bases = ["ah1w3f0", "ah1w3r0", "ah1w3"]
snapnums =[23, 22, 157]

for m in np.arange(3):

   loc_start = rank_start + m

#   snapnum = snapnum_start + loc_start*every_snapshot
   snapnum = snapnums[m]
   base = bases[m]

#   print 'snapnum ', snapnum, ', loc_start', loc_start, '\n'
#   print 'rank_start', rank_start

   path, snapnum, snapfile = setup_snapfile(filesystem, snapnum, base, ending)

   f = open(snapfile, mode = 'rb')
#   print 'snapfile: ', snapfile

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

   figsize = get_figsize(screen_ratio, figheight)

   if loc_start == 0:
      fig = plt.figure(figsize = figsize)

   xlabel = get_label(length_unit, 3)

   for i in np.arange(num_blocks):

      yblock = f.read(max_string_len).replace("\x00", "")

      ylabel = get_label(yblock)

      yval = np.fromfile(f, dtype = 'd', count = bins)

      ax = fig.add_subplot(num_rows, num_columns, i + 1)

      l1 = ax.plot(xval, yval, label=snapnum) #, 'black')

#      ax.set_xlim([xmin, xmax])

      ax.set_xlabel(xlabel, fontsize = 12)
      ax.set_ylabel(ylabel, fontsize = 12)
     
#      if m==numsnaps_rank/every_snapshot and i==num_blocks-1:
#         ax.legend(bbox_to_anchor=(1.05,0), loc='lower left', borderaxespad=0.)

   f.close()

   fig.tight_layout()

   print 'PYTHON: Snap', loc_start, 'done!'

#print 'hola2'

write_image(path, base, snapnum, numsnaps, 'out')
