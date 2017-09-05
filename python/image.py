
from all import *
from mpl_toolkits.axes_grid1 import AxesGrid

filesystem = 1
ending = 'img'

figheight = 15.
screen_ratio = 1.

num_sub_snaps = 1

fontsize = 8

ax = 0

mp.rcParams['font.size'] = fontsize

figsize = get_figsize(screen_ratio, figheight)

voro_flag, base, snapnum_start, snapnum_end, num_snaps, num_snaps_rank, rank_start = get_args(sys.argv, 1)

for i in np.arange(num_snaps_rank):

   loc_start = rank_start + i

   snapnum = snapnum_start + loc_start

   for j in np.arange(num_sub_snaps):

      path, snap_string, snapfile = setup_snapfile(filesystem, snapnum, base, ending, num_sub_snaps)

      if voro_flag:

         snapfile = path + '/' + base + '/' + 'density_proj_' + snapnum

         f = open(snapfile, mode = 'rb')

         f.seek(4)

         num_blocks = 1
         bins = struct.unpack('i', f.read(4))[0]
         num_sinks = 0
         width = 1.

      else:

         f = open(snapfile, mode = 'rb')

         num_blocks = struct.unpack('i', f.read(4))[0]
         max_string_len = struct.unpack('i', f.read(4))[0]
         string_len = str(max_string_len) + 's'
         xbins = struct.unpack('i', f.read(4))[0]
         ybins = struct.unpack('i', f.read(4))[0]
         num_sinks = struct.unpack('i', f.read(4))[0]
         width = struct.unpack('d', f.read(8))[0]
         height = struct.unpack('d', f.read(8))[0]
         time = struct.unpack('d', f.read(8))[0]
         length_unit = struct.unpack('i', f.read(4))[0]
         comoving_units = struct.unpack('i', f.read(4))[0]

         num_sinks = 1
         if num_sinks:

            sinkpos = np.fromfile(f, dtype = 'd', count = 2 * num_sinks)

            sinkx = sinkpos[0 : : 2]
            sinky = sinkpos[1 : : 2]
            print sinkx, sinky

      num_columns, num_rows = get_columns_and_rows(num_blocks, screen_ratio)

      fig = plt.figure(figsize = figsize)

      grid = AxesGrid(fig, 111, nrows_ncols = [num_rows, num_columns], ngrids = num_blocks, axes_pad = 0.4, share_all = True, cbar_location = 'top', cbar_mode = 'each', cbar_pad = '0%')

      for k in np.arange(num_blocks):

         if voro_flag:
            block = 'temp'
         else:
            block = f.read(max_string_len).replace("\x00", "")

         print block
         label = get_label(block)

         xlabel = get_label(length_unit, 0)
         ylabel = get_label(length_unit, 1)

         cmap = get_color_map(block)

         cmap.set_bad('w')

         val_min_dim = []
         val_max_dim = []
         val_dim = []

         for dim in np.arange(3):

            val_min_dim.append(struct.unpack('d', f.read(8))[0])
            val_max_dim.append(struct.unpack('d', f.read(8))[0])
            val_dim.append(np.fromfile(f, dtype = 'd', count = xbins * ybins))

         val_min = val_min_dim[ax]
         val_max = val_max_dim[ax]
         val = val_dim[ax]

         if voro_flag:
            val[val != 0] = np.log10(val[val != 0])

         val = np.reshape(val, (xbins, ybins))
         val = np.rot90(val)
         val = np.ma.masked_where(val == 0, val)

         im = grid[k].imshow(val, cmap = cmap, extent = [0, width, 0, height]) # , vmin = val_min, vmax = val_max)

         #if num_sinks:

            #grid[k].plot(sinkx, sinky, 'wo')
            #grid[k].set_xlim([xval_min, xval_max])
            #grid[k].set_ylim([yval_min, yval_max])

         grid.cbar_axes[k].colorbar(im)
         grid.cbar_axes[k].axis["right"].toggle(ticks = True, ticklabels = True, label = True)
         grid.cbar_axes[k].set_ylabel(label, va = 'top', fontsize = fontsize + 1)

         grid[k].set_xlim(0, 10)
         grid[k].set_ylim(0, 10)
         grid[k].set_xlabel(xlabel, fontsize = fontsize + 1)
         grid[k].set_ylabel(ylabel, fontsize = fontsize + 1)

      f.close()

      write_image(path, base, snap_string, num_snaps_rank, ending, num_sub_snaps, j)

      plt.close()

      print 'PYTHON: #', snapnum_start + i, ',', j, 'of', num_snaps_rank * num_sub_snaps, 'done!'
