
from all import *

imgmap_src = os.environ['HOME'] + '/arepo/sator/swig/imgmap.py'
imgmap = imp.load_source('imgmap', imgmap_src)

filesystem = 1
ending = 'img'

width = 25
length_unit = 'pc'
comoving_units = 0
bins = 300

var_flag = 0

val_min = 0
val_max = 0

figwidth = 28
figheight = 28
fontsize = 12
labelfac = 1.2

comm.Barrier()

t1 = time.clock()

base, snapnum_start, snapnum_end, numsnaps, numsnaps_rank, rank_start = get_args(sys.argv)

for m in np.arange(numsnaps_rank):

   loc_start = rank_start + m

   snapnum = snapnum_start + loc_start

   if numsnaps > 1: 
      verbose = 0
   else:
      verbose = 1

#    if numsnaps > 1: 

#       verbose = 0

#       s = readsnap(base, snapnum_end, length_unit = length_unit, comoving_units = comoving_units, verbose = verbose)

#       s.readblocks([0], ['pos', 'rho'])

#       cx, cy, cz = com(s)

#       print snapnum

#    else: 

#       verbose = 1

   s = readsnap(filesystem, base, snapnum, length_unit = length_unit, comoving_units = comoving_units, verbose = verbose)

   if s.npartall[0] == 0:
      read_snaptypes = [1]
   else:
      read_snaptypes = [0]

   #read_snaptypes = [0, 1]

   read_blocks = ['pos', 'id', 'mass', 'u', 'rho', 'vol']#, 'chem', 'gamma']

   s.readblocks(read_snaptypes, read_blocks)

#   print min(s.x), max(s.x)
#   print min(s.y), max(s.y)
#   print min(s.z), max(s.z)
#   print min(s.vx), max(s.vx)
#   print min(s.vy), max(s.vy)
#   print min(s.vz), max(s.vz)
#   print min(s.mass), max(s.mass)
#   print min(s.id), max(s.id)
   print max(s.nh)
   print s.id[s.nh == max(s.nh)]
#   print min(s.temp), max(s.temp)
#   print min(s.hsml), max(s.hsml)

#    if numsnaps == 1:

#       if s.npartall[0] == 0:

#          cx = cy = cz = s.boxsize / 2.

#       else:

   cx, cy, cz = com(s)

   #cx = 344.93085172 / s.hubbleparam
   #cy = 357.03698736 / s.hubbleparam
   #cz = 347.82869235 / s.hubbleparam

   #cx = 0.5 * s.boxsize
   #cy = 0.5 * s.boxsize
   #cz = 0.5 * s.boxsize

#    idf = 1005603351

#    idxh = s.id == idf

#    cx = s.x[idxh]
#    cy = s.y[idxh]
#    cz = s.z[idxh]

#    print s.id[idxh], cx, cy, cz

   s.x -= cx
   s.y -= cy
   s.z -= cz

#    print s.nh[s.id == idf]
#    print s.temp[s.id == idf]
#    print s.mass[s.id == idf]

#    rad = np.sqrt(s.x**2 + s.y**2 + s.z**2)

#    idx = rad < 0.0001

#    print s.x[idx]


#    idx = rad > 0.1

#    nhmax = max(s.nh[idx])

#    print nhmax

#    idx = s.nh == nhmax

#    print s.id[idx]


#    if(snapnum == 48):

#       xx = s.x * s.hubbleparam
#       yy = s.y * s.hubbleparam
#       zz = s.z * s.hubbleparam   

#       newrad = np.sqrt(xx**2 + yy**2 + zz**2)



#       idxbla = np.logical_and(newrad > 0.005, newrad < 0.01)

#       id_list = s.id[idxbla]
#       temp_list = s.temp[idxbla]

#       print id_list
#       print temp_list

#    else:

#       fac = 0.
#       maxfac = 0.

#       c = np.in1d(s.id, id_list)

#       id_list2 = s.id[c]
#       temp_list2 = s.temp[c]

#       if id_list.size != id_list2.size:

#          print 'something wrong 2!'

#       for i in np.arange(id_list.size):

#          for j in np.arange(id_list.size):

#             if id_list[i] == id_list2[j]:

#                if temp_list2[j] >= temp_list[i]:

#                   fac = temp_list2[j] / temp_list[i]

#                else:

#                   fac = temp_list[i] / temp_list2[j]


#                if fac > maxfac:

#                   maxfac = fac

#                   maxi = i
#                   maxj = j

#       print maxfac, id_list[maxi], temp_list[maxi], id_list2[maxj], temp_list2[maxj]


   idx = np.logical_and(np.logical_and(np.fabs(s.x) < width / 2., np.fabs(s.y) < width / 2.), np.fabs(s.z) < width / 2.)

   s.x = s.x[idx]
   s.y = s.y[idx]
   s.z = s.z[idx]

   s.mass = s.mass[idx]

   if s.npartall[0] == 0:
      s.nh = 1.
      var = s.mass
      s.hsml = 50 * 0.017 / s.hubbleparam
   else:
      s.nh = s.nh[idx]
      s.hsml = s.hsml[idx]

   if var_flag == 0:
      var = s.nh
      cmap = mp.cm.jet
   elif var_flag == 1:
      var = s.temp
      cmap = mp.cm.gist_heat
   else:
      print 'Unknown variable! Exiting...'
      comm.Abort()

   if var[0] != s.nh[0]:
      var = var[idx]

   imin = np.maximum(np.minimum(((s.x + width / 2. - s.hsml) / width * bins).astype('int32'), bins - 1), 0)
   imax = np.maximum(np.minimum(((s.x + width / 2. + s.hsml) / width * bins).astype('int32'), bins - 1), 0)

   jmin = np.maximum(np.minimum(((s.y + width / 2. - s.hsml) / width * bins).astype('int32'), bins - 1), 0)
   jmax = np.maximum(np.minimum(((s.y + width / 2. + s.hsml) / width * bins).astype('int32'), bins - 1), 0)

   img_part = s.nh * s.nh * s.mass * var

   sum_part = s.nh * s.nh * s.mass

   img = np.zeros([bins, bins])
   sum = np.zeros([bins, bins])

   imgmap.range(imin, imax, jmin, jmax, img_part, sum_part, img, sum)

#       if(snapnum == 2):

#          img_save = img

#       else:

#          new = img / img_save

#          print np.max(new)

      #img = comm.reduce(img, op = MPI.SUM)
      #sum = comm.reduce(sum, op = MPI.SUM)

      #if rank == 0:

   idx_valid = img != 0

   img[idx_valid] = np.log10(img[idx_valid] / sum[idx_valid])

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

   img = np.ma.masked_where(img == 0, img)

   cmap.set_bad('w')

   labelsize = labelfac * fontsize

   mp.rcParams['font.size'] = fontsize

   figsize = [figwidth * cm_to_inch, figheight * cm_to_inch]

   fig = plt.figure(figsize = figsize)

   ax = fig.add_subplot(111)

   #ax.set_xlabel(r'${\rm log}\,n_{\rm H}\,[{\rm cm}^{-3}]$', fontsize = labelsize)
   #ax.set_ylabel(r'${\rm log}\,T\,[{\rm K}]$', fontsize = labelsize)

   im = ax.imshow(np.rot90(img), cmap = cmap, vmin = vmin, vmax = vmax, extent = [0, width, 0, width], aspect = 'auto')#, interpolation = 'none')

   #fig.colorbar(im)

   fig.tight_layout()

   write_image(s.path, base, loc_start, snapnum, numsnaps, ending)

   if rank == 0 and numsnaps > 1:

      print 'Iter ', m, ' done!'

comm.Barrier()

t2 = time.clock()

if rank == 0: print 'took', t2 - t1, 'secs'
