
from all import *

filesystem = 1
ending = 'mbe'

figheight = 10.
screen_ratio = 2.

fontsize = 7.

xmin = 0
xmax = 100

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
   dist = struct.unpack('d', f.read(8))[0]
   xval_min = struct.unpack('d', f.read(8))[0]
   xval_max = struct.unpack('d', f.read(8))[0]
   xval = np.fromfile(f, dtype = 'd', count = bins)

#   dist = 10**dist
#   xval_min = 10**xval_min
#   xval_max = 10**xval_max
#   xval = 10**xval

   num_columns, num_rows = get_columns_and_rows(num_blocks, screen_ratio)

   figsize = get_figsize(screen_ratio, figheight)

   fig = plt.figure(figsize = figsize)

   xlabel = get_label(length_unit, 3)
   print dist
   for i in np.arange(num_blocks):

      yblock = f.read(max_string_len).replace("\x00", "")

      ylabel = get_label(yblock)

      yval = np.fromfile(f, dtype = 'd', count = bins)
      print yval

      ax = fig.add_subplot(num_rows, num_columns, i + 1)

      ax.plot(xval, yval, 'black')

      ax.plot([xval_min, xval_max], [0, 0], 'grey', linestyle = '--')

      if i == 0:
         ax.plot([dist, dist], [max(yval) - 0.1, max(yval)], 'black')

      ax.set_xlim([xval_min, xval_max])

      #if i == 1:
         #ax.set_ylim([xmin, xmax])

      ax.set_xlabel(xlabel)
      ax.set_ylabel(ylabel)

   f.close()

   fig.tight_layout()

   write_image(path, base, snapnum, numsnaps, ending)

   print 'PYTHON: Snap', loc_start, 'done!'
