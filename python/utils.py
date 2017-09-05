
from include import *


def block_exit():

    print 'Unknown block! Exiting...'

    comm.Abort()


def get_args(argv, img_flag = 0, base = '', snapnum_start = 0, snapnum_end = 0):

    args = len(argv)

    if img_flag:
        req_args = 4
    else:
        req_args = 3

    if args < req_args:
        print 'Not enough arguments specified! Exiting...'
        comm.Abort()

    if img_flag:
        voro_flag = int(argv[1])

    base = argv[req_args - 2]

    snapnum_start = snapnum_end = int(argv[req_args - 1])

    if args > req_args:
        snapnum_end = int(argv[req_args])

    if snapnum_end < snapnum_start: snapnum_end = snapnum_start

    num_snaps = snapnum_end - snapnum_start + 1

    num_snaps_rank = num_snaps / size

    rank_start = rank * num_snaps_rank

    rest = num_snaps - size * num_snaps_rank

    if rank < rest:
        num_snaps_rank += 1
        rank_start += rank
    else:
        rank_start += rest

    if img_flag:
        return voro_flag, base, snapnum_start, snapnum_end, num_snaps, num_snaps_rank, rank_start
    else:
        return base, snapnum_start, snapnum_end, num_snaps, num_snaps_rank, rank_start


def get_figsize(screen_ratio, figheight):

    figsize = [screen_ratio * figheight * cm_to_inch, figheight * cm_to_inch]

    return figsize


def com(s):

    rank_max = -1

    cx = -1.
    cy = -1.
    cz = -1.

    loc_nh_max = s.nh.max()

    glob_nh_max = comm.allreduce(loc_nh_max, op = MPI.MAX)

    if loc_nh_max == glob_nh_max:

        rank_max = rank

        idx = s.nh == loc_nh_max

        cx = s.x[idx]
        cy = s.y[idx]
        cz = s.z[idx]

        #print s.id[idx], s.nh[idx]

    rank_max = np.array(comm.allgather(rank_max))

    rank_max = rank_max[rank_max > -1][0]

    cx = comm.bcast(cx, rank_max)
    cy = comm.bcast(cy, rank_max)
    cz = comm.bcast(cz, rank_max)

    return cx, cy, cz


def register_own_cmaps():

    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 0.0, 0.0),
                       (0.5, 1.0, 1.0),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),
                      (0.5, 1.0, 1.0),
                      (1.0, 0.0, 0.0))}

    cmap = mp.colors.LinearSegmentedColormap('my_colormap', cdict, N = 256)

    mp.cm.register_cmap(name = 'own1', cmap = cmap)


def set_img_min_max(val, val_min, val_max):

    if val_min == 0:

        vmin = np.amin(val[val != 0])

    else:

        val[np.logical_and(val < np.log10(val_min), val != 0)] = np.log10(val_min)

        vmin = np.log10(val_min)

    if val_max == 0:

        vmax = np.amax(val[val != 0])

    else:

        val[np.logical_and(val > np.log10(val_max), val != 0)] = np.log10(val_max)

        vmax = np.log10(val_max)

    return val, vmin, vmax


def setup_snapfile(filesystem, snapnum, base, ending = '', num_sub_snaps = 1, subnum = 0):

    systype = os.environ['SYSTYPE']

    if systype == 'odyssey-opteron':
        path = '/n/hernquistfs2/fbecerra'
    elif systype == 'odin':
        path = '/ptmp/mpa/tgreif'
    elif systype == 'stampede':
#        path = '/scratch/02563/fbecerra'
#        path = '/scratch/02563/fbecerra/images'
        path = '/scratch/00025/tgreif'
    elif systype == 'lonestar':
        path = '/scratch/00025/tgreif'
    else:
        print 'System type not recognized! Exiting...'
        comm.Abort()

    if snapnum == -1:
        snap_string = 'ic'
    else:
        snap_string = repr(snapnum).zfill(3)

    if ending != '':
        ending = '.' + ending

    if num_sub_snaps == 1:
        sub_string = ''
    else:
        sub_string = '_' + repr(subnum).zfill(3)

    snapfile = path + '/' + base + '/' + base + '_' + snap_string + sub_string + ending
#    snapfile = path + '/' + base + '_' + snap_string + sub_string + ending

    return path, snap_string, snapfile


def check_snapnum_range(snapnum):

    if snapnum < 0 or snapnum >= 10000:

        if rank == 0: print 'snap_num =', snapnum, 'not allowed! Exiting...'

        comm.Abort()

    else:

        return snapnum


def get_color_map(block):

    if block == 'nh':
        return mp.cm.jet
    elif block == 'temp':
        return mp.cm.gist_heat 
    elif block == 'gravacc':
        return mp.cm.jet
    elif block == 'gradp':
        return mp.cm.jet
    elif block == 'abhm':
        return mp.cm.gist_stern
    elif block == 'abh2':
        return mp.cm.gist_stern
    elif block == 'abhii':
        return mp.cm.get_cmap('own1')
    elif block == 'abdii':
        return mp.cm.gist_stern
    elif block == 'abhd':
        return mp.cm.gist_stern
    elif block == 'abheii':
        return mp.cm.gist_stern
    elif block == 'abheiii':
        return mp.cm.gist_stern
    elif block == 'gamma':
        return mp.cm.jet
    elif block == 'pdvrate':
        return mp.cm.gist_stern
    elif block == 'h2rate':
        return mp.cm.gist_stern
    elif block == 'cierate':
        return mp.cm.gist_stern
    elif block == 'chemrate':
        return mp.cm.gist_stern
    elif block == 'escfrac':
        return mp.cm.gist_stern
    elif block == 'allowref':
        return mp.cm.gist_stern
    elif block == 'divvel':
        return mp.cm.gist_stern
    elif block == 'cool':
        return mp.cm.gist_stern
    elif block == 'collapse':
        return mp.cm.get_cmap('own1')
    elif block == 'entro':
        return mp.cm.gist_stern
    elif block == 'sigma_mass':
        return mp.cm.gist_stern
    elif block == 'csnd':
        return mp.cm.gist_stern
    elif block == 'omega_rot':
        return mp.cm.gist_stern
    elif block == 'q':
        return mp.cm.gist_stern
    elif block == 'sidm_density':
        return mp.cm.jet
    else:
        block_exit()


def get_columns_and_rows(num_blocks, screen_ratio):

    num_columns_real = float(screen_ratio) * np.sqrt(float(num_blocks))
    num_columns = int(num_columns_real)

    if num_columns_real > num_columns:
        num_columns += 1

    if num_columns > num_blocks:
        num_columns = num_blocks

    num_rows_real = float(num_blocks) / float(num_columns)
    num_rows = int(num_rows_real)

    if num_rows_real > num_rows:
        num_rows += 1

    return num_columns, num_rows


def get_label(block, dir = 0):

    if block == 0:
        if dir == 0:
            label = r'$x\,[kpc]$'
        elif dir == 1:
            label = r'$y\,[kpc]$'
        else:
            label = r'$r\,[kpc]$'
    elif block == 1:
        if dir == 0:
            label = r'$x\,[pc]$'
        elif dir == 1:
            label = r'$y\,[pc]$'
        else:
            label = r'$r\,[pc]$'
    elif block == 2:
        if dir == 0:
            label = r'$x\,[au]$'
        elif dir == 1:
            label = r'$y\,[au]$'
        else:
            label = r'$r\,[au]$'
    elif block == 'enc':
        label = r'${\rm log}\,M_{\rm enc}\,[{\rm M}_\odot]$'
    elif block == 'mbe_ratio':
        label = r'${\rm log}\left(M_{\rm enc}/M_{\rm BE}\right)$'
    elif block == 'acc_ratio':
        label = r'${\rm log}\left(t_{\rm ff}/t_{\rm acc}\right)$'
    elif block == 'vrad':
#        label = r'$v_{\rm rad}\,[{\rm km}\,{\rm s}^{-1}]$'
#    elif block == 'accrate':
        label = r'$\dot{M}\,[{\rm M}_\odot / {\rm yr}]$' 
    elif block == 'frac_rad':
        label = r'$v_{\rm rad}/c_{\rm s}$'
    elif block == 'vrot':
        label = r'$v_{\rm rot}\,[{\rm km}\,{\rm s}^{-1}]$'
    elif block == 'frac_rot':
        label = r'$v_{\rm rot}/c_{\rm s}$'
    elif block == 'vturb':
        label = r'$v_{\rm turb}\,[{\rm km}\,{\rm s}^{-1}]$'
    elif block == 'frac_turb':
        label = r'$v_{\rm turb}/\left|v_{\rm rad}\right|$'
    elif block == 'kep_ratio':
        label = r'$v_{\rm rot}/v_{\rm kep}$'
    elif block == 'mach':
        label = r'$\mathcal{M}$'
    elif block == 'nh':
        label = r'${\rm log}\,n_{\rm H}\,[{\rm cm}^{-3}]$'
    elif block == 'dnh':
        label = r'${\rm log}\,\sigma_\delta$'
    elif block == 'temp':
        label = r'${\rm log}\,T\,[{\rm K}]$'
    elif block == 'gravacc':
        label = r'$a_{\rm grav}$'
    elif block == 'gradp':
        label = r'$\nabla P$'
    elif block == 'abhm':
        label = r'${\rm log}\,{\rm y}_{{\rm H}^-}$'
    elif block == 'abh2':
        label = r'${\rm log}\,{\rm y}_{{\rm H}_2}$'
    elif block == 'abhii':
        label = r'${\rm log}\,{\rm y}_{\rm HII}$'
    elif block == 'abdii':
        label = r'${\rm log}\,{\rm y}_{\rm DII}$'
    elif block == 'abhd':
        label = r'${\rm log}\,{\rm y}_{\rm HD}$'
    elif block == 'abheii':
        label = r'${\rm log}\,{\rm y}_{\rm HeII}$'
    elif block == 'abheiii':
        label = r'${\rm log}\,{\rm y}_{\rm HeIII}$'
    elif block == 'gamma':
        label = r'$\gamma$'
    elif block == 'geff':
        label = r'$\gamma_{\rm eff}$'
    elif block == 'escfrac':
        label = r'${\rm log}\,f_{\rm esc}$'
    elif block == 'pdvrate':
        label = r'${\rm log}\left(t_{\rm pdv}/t_{\rm ff}\right)$'
    elif block == 'h2rate':
        label = r'${\rm log}\left(t_{\rm H_2}/t_{\rm ff}\right)$'
    elif block == 'cierate':
        label = r'${\rm log}\left(t_{\rm CIE}/t_{\rm ff}\right)$'
    elif block == 'chemrate':
        label = r'${\rm log}\left(t_{\rm 3b}/t_{\rm ff}\right)$'
    elif block == 'allowref':
        label = r'${\rm AllowRef}$'
    elif block == 'divvel':
        label = r'$\nabla V$'
    elif block == 'cool':
        label = r'${\rm log}\left(t_{\rm cool}/t_{\rm ff}\right)$'
    elif block == 'collapse':
        label = r'${\rm log}\left(t_{\rm ff}/t_{\rm c_s}\right)$'
    elif block == 'entro':
        label = r'$S$'
    elif block == 'sigma_mass':
        label = r'${\rm log}\,\Sigma\,[{\rm g}\,{\rm cm}^{-2}]$'
    elif block == 'csnd':
        label = r'${\rm c}_{\rm s}\,[{\rm km}\,{\rm s}^{-1}]$'
    elif block == 'omega_rot':
        label = r'${\rm log}\,\Omega\,[{\rm yr}^{-1}]$'
    elif block == 'q':
        label = r'${\rm Q}$'
    elif block == 'sidm_density':
        label = r'${\rm log}\,n_{\rm H}\,[{\rm cm}^{-3}]$'
    else:
        block_exit()

    return label


def get_min_max_val(block):

    if block == 'temp':
        val_min = 2.1
        val_max = 2.9
    else:
        val_min = 0
        val_max = 0

    return val_min, val_max


def get_plot_block():

    f = open('py_plot_blocks.txt', mode = 'r')

    columns = []

    for line in f:

        read_block, block_flag = line.split()

        block_flag = int(block_flag)

        columns.append([read_block, block_flag])

    plot_block = dict(columns)

    f.close()

    return plot_block


def mpi_print(*vals):

    for i in range(size):

        if i == rank:

            for val in vals:

                print val,

            print

        comm.Barrier()


def write_image(path, base, snap_string, num_snaps, ending, num_sub_snaps = 1, subnum = 0):

    if num_sub_snaps == 1:
        sub_string = ''
    else:
        sub_string = '_' + repr(subnum).zfill(3)

#    string = path + '/' + base + '/' + base + '_' + snap_string + sub_string + '_' + ending
    string = path + '/' + base + '_' + snap_string + sub_string + '_' + ending

    if num_snaps == 1 and num_sub_snaps == 1:

        #outfile_eps = ending + '.eps'
        outfile_pdf = 'data/' + ending + '.pdf'
        outfile_png = 'data/' + ending + '.png'
        outfile_jpg = 'data/' + ending + '.jpg'

    else:

        #outfile_eps = ending + '.eps'
        outfile_pdf = string + '.pdf'
        outfile_png = string + '.png'
        outfile_jpg = string + '.jpg'

    #plt.savefig(outfile_eps)
    plt.savefig(outfile_pdf)
    plt.savefig(outfile_png)
#    Image.open(outfile_png).save(outfile_jpg, quality = 100)

#    rm_png = 'rm -f ' + outfile_png

    #os.system(rm_png)
