
from all import *

filesystem = 0
ending = 'mov'

fig_width = 12.
fontsize = 12.

args = len(sys.argv)

if args < 6:
   print 'Not enough arguments specified! Exiting...'
   comm.Abort()

base = sys.argv[1]
snapnum_start = int(sys.argv[2])
snapnum_end = int(sys.argv[3])

if snapnum_end < snapnum_start: snapnum_end = snapnum_start

num_snaps = snapnum_end - snapnum_start + 1

num_sub_snaps = int(sys.argv[4])
target_block = sys.argv[5]

if target_block == 'nh':
   vmin_init = float(sys.argv[6])
   vmax_init = float(sys.argv[7])
   vmin_target = float(sys.argv[8])
   vmax_target = float(sys.argv[9])

if base == 'mh1x2':
   flag_cosm = 1
else:
   flag_cosm = 0

hubble_param = 0.7
omega_m = 0.27
redshift_cosm = 23.8925390917
dt_cosm = 277.775778778722326

ax = 0

mp.rc('axes', edgecolor = 'w')
mp.rc('axes', labelcolor = 'w')
mp.rc('xtick', c = 'w')
mp.rc('ytick', c = 'w')
mp.rc('font', size = fontsize)

if snapnum_end < snapnum_start:
   snapnum_end = snapnum_start

num_snaps = snapnum_end - snapnum_start + 1

for i in np.arange(num_snaps):

   snapnum = snapnum_start + i

   for j in np.arange(num_sub_snaps):

      path, snap_string, snapfile = setup_snapfile(filesystem, snapnum, base, ending, num_sub_snaps, j)

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

      for k in np.arange(num_blocks):

         block = f.read(max_string_len).replace("\x00", "")

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

         if block == target_block:

            label = get_label(block)

            xlabel = r'$L\,[pc]$'

            cmap = get_color_map(block)

            cmap.set_bad('w')

            if block == 'nh':
               vmin = vmin_init + j * (vmin_target - vmin_init) / max((num_sub_snaps - 1), 1)
               vmax = vmax_init + j * (vmax_target - vmax_init) / max((num_sub_snaps - 1), 1)

            if block == 'temp':
               val[val != 0] = 10**val[val != 0]
               vmin = 10.
               vmax = 2000.
               label = r'$T\,[{\rm K}]$'

            val = np.reshape(val, (xbins, ybins))
            val = np.rot90(val)
            val = np.ma.masked_where(val == 0, val)

            fig = plt.figure(figsize = [fig_width, fig_width * height / width])

            axes = fig.add_axes([0, 0, 1, 1])

            axes.get_xaxis().set_visible(False)
            axes.get_yaxis().set_visible(False)

            im = axes.imshow(val, cmap = cmap, extent = [0, width, 0, height], vmin = vmin, vmax = vmax)

            logo = Image.open('logo.jpg').rotate(180).transpose(0)

            axes = fig.add_axes([0.023, 0.004, 0.12, 0.12])

            axes.get_xaxis().set_visible(False)
            axes.get_yaxis().set_visible(False)

            axes.imshow(logo)

            scale_width = 0.25

            axes = fig.add_axes([0.72, 0.07, scale_width, 0.0])

            axes.get_yaxis().set_visible(False)

            axes.xaxis.tick_bottom()

            axes.set_xlim([0, scale_width * width])

            axes.ticklabel_format(style = 'sci', scilimits = (0, 0), axis = 'x', fontsize = fontsize)

            for ticklabel in axes.get_xticklabels():
               ticklabel.set_fontsize(fontsize - 2)

            axes.set_xlabel(xlabel, fontsize = fontsize + 1)

            caxes = fig.add_axes([0.93, 0.64, 0.02, 0.3])

            cbar = fig.colorbar(im, cax = caxes)

            for ticklabel in caxes.get_yticklabels():
               ticklabel.set_fontsize(fontsize - 2)

            caxes.xaxis.set_label_position('top')

            caxes.set_xlabel(label, labelpad = 12, fontsize = fontsize + 1)

            fig.text(0.031, 0.11, r'${\rm Thomas\,H.\,Greif}$', color = 'w', fontsize = fontsize)

            fig.text(0.34, 0.96, r'${\rm Collapse\,of\,primordial\,star-forming\,cloud}$', color = 'w', fontsize = fontsize + 4)

            if flag_cosm:
               redshift = 1. / time - 1.
            else:
               redshift = redshift_cosm

            string = r'$z\,=\,$' + r'$%.2f$' % redshift

            fig.text(0.02, 0.96, string, color = 'w', fontsize = fontsize + 1)

            if flag_cosm:
               dt = 2. / 3. / hubble / hubble_param / omega_m * time**1.5 / sec_per_year / 1e6
            else:
               dt = dt_cosm + time * unit_time / sec_per_year / 1e6

            string = r'$t_{\rm H}\,=\,$' + r'$%.5f$' % dt + r'$\,{\rm Myr}$'

            fig.text(0.02, 0.92, string, color = 'w', fontsize = fontsize + 1)

      f.close()

      write_image(path, base, snap_string, num_snaps, ending, num_sub_snaps, j)

      print 'PYTHON: #', snapnum_start + i, ',', j, 'of', num_snaps * num_sub_snaps, 'done!'
