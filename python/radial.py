#!/opt/apps/python/epd/7.3.2/bin/python

from all import *

filesystem = 1
ending = 'radial'

figheight = 20.
screen_ratio = 1.

fontsize = 7.

mp.rcParams['font.size'] = fontsize

base, snapnum_start, snapnum_end, numsnaps, numsnaps_rank, rank_start = get_args(sys.argv)

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

   figsize = get_figsize(screen_ratio, figheight)

   fig = plt.figure(figsize = figsize)

   xlabel = get_label(length_unit, 3)

   for i in np.arange(num_blocks):

      yblock = f.read(max_string_len).replace("\x00", "")

      ylabel = get_label(yblock)

      yval = np.fromfile(f, dtype = 'd', count = bins)

      ax = fig.add_subplot(num_rows, num_columns, i + 1)

      ax.plot(xval, yval, 'black')

      if yblock == 'q':
         ax.plot([xmin, xmax], [1, 1], 'grey', linestyle = '--')

      ax.set_xlim([xmin, xmax])

      ax.set_xlabel(xlabel, fontsize=12)
      ax.set_ylabel(ylabel, fontsize=12)
     
   f.close()

   fig.tight_layout()

   write_image(path, base, snapnum, numsnaps, 'out')

   plt.close()

   print 'PYTHON: Snap', loc_start, 'done!'

