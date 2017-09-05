
from all import *

filesystem = 1
ending = 'pspace'

figheight = 20.
screen_ratio = 1.

fontsize = 7

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
   xblock = f.read(max_string_len).replace("\x00", "")
   xval_min = struct.unpack('d', f.read(8))[0]
   xval_max = struct.unpack('d', f.read(8))[0]
   xval = np.fromfile(f, dtype = 'd', count = bins)

   num_columns, num_rows = get_columns_and_rows(num_blocks, screen_ratio)

   figsize = get_figsize(screen_ratio, figheight)

   fig = plt.figure(figsize = figsize)

   for i in np.arange(num_blocks):

      yblock = f.read(max_string_len).replace("\x00", "")

      xlabel = get_label(xblock)
      ylabel = get_label(yblock)

      cmap = mp.cm.jet

      cmap.set_bad('w')

      yval_min = struct.unpack('d', f.read(8))[0]
      yval_max = struct.unpack('d', f.read(8))[0]

      yval = np.fromfile(f, dtype = 'd', count = bins)

      yval = np.ma.masked_where(yval == min_real_number, yval)

      val = np.fromfile(f, dtype = 'd', count = bins**2)

      val = np.reshape(val, (bins, bins))
      val = np.rot90(val)
      val = np.ma.masked_where(val == min_real_number, val)

      ax = fig.add_subplot(num_rows, num_columns, i + 1)

      im = ax.imshow(val, cmap = cmap, extent = [xval_min, xval_max, yval_min, yval_max], aspect = 'auto')

      ax.plot(xval, yval, 'black')

      ax.set_xlabel(xlabel, fontsize=12)
      ax.set_ylabel(ylabel, fontsize=12)

   f.close()

   fig.tight_layout()

   write_image(path, base, snapnum, numsnaps, 'out')

   print 'PYTHON: Snap', loc_start, 'done!'
