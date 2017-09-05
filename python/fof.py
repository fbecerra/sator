
from include import *
from utils import *

class readfof:

    def __init__(self, filesystem = 0, base = '/', fofnum = 0, flag_double = 1, length_unit = 'kpc', comoving_units = 0):

        self.base = base
        self.flag_double = flag_double

        self.path, self.fofnum, self.snapfile = setup_snapfile(filesystem, fofnum, base)

        self.foffile = self.fofbase = self.path + '/' + self.base + '/fof_tab_' + self.fofnum

        if not os.path.isfile(self.foffile):

            if self.fofnum != 'ic':
                self.fofbase = self.path + '/' + self.base + '/' + 'groups' + '_' + self.fofnum + '/fof_tab_' + self.fofnum

            self.foffile = self.fofbase + '.0'

            if not os.path.isfile(self.foffile):

                print 'Could not find fof file! Exiting...'

                comm.Abort()

        self.length_unit = length_unit
        self.comoving_units = comoving_units

        self.readheader()

    def readheader(self):

        f = open(self.foffile, mode = 'rb')

        self.nbytes = struct.unpack('i', f.read(4))[0]

        self.ngroups = struct.unpack('i', f.read(4))[0]
        self.nsubgroups = struct.unpack('i', f.read(4))[0]
        self.nids = struct.unpack('i', f.read(4))[0]
        self.totngroups = struct.unpack('i', f.read(4))[0]
        self.totnsubgroups = struct.unpack('i', f.read(4))[0]

        dummy = struct.unpack('i', f.read(4))[0]
        self.totnids = struct.unpack('i', f.read(4))[0]
        dummy2 = struct.unpack('i', f.read(4))[0]
        self.numfiles  = struct.unpack('i', f.read(4))[0]
        dummy3 = struct.unpack('i', f.read(4))[0]

        self.time = struct.unpack('d', f.read(8))[0]
        self.redshift = struct.unpack('d', f.read(8))[0]
        self.hubbleparam = struct.unpack('d', f.read(8))[0]
        self.boxsize = struct.unpack('d', f.read(8))[0]
        self.omega_m = struct.unpack('d', f.read(8))[0]
        self.omega_lambda = struct.unpack('d', f.read(8))[0]
        self.flag_doubleprecision = struct.unpack('i', f.read(4))[0]

#        print self.foffile, self.nbytes, self.ngroups, self.nsubgroups, self.nids, self.totngroups, self.totnsubgroups, self.totnids, dummy, dummy2, self.numfiles, dummy3
#        print self.time, self.redshift, self.hubbleparam, self.boxsize, self.omega_m, self.omega_lambda, self.flag_doubleprecision

        f.close()

    def readgroups(self):

        if self.flag_doubleprecision == 0:
            float_type = 'f'
        else:
            float_type = 'd'

        self.len = np.zeros(self.totngroups, 'i')
        self.mass = np.zeros(self.totngroups, float_type)
        self.pos = np.zeros(3 * self.totngroups, float_type)
        self.cm = np.zeros(3 * self.totngroups, float_type)
        self.vel = np.zeros(3 * self.totngroups, float_type)
        self.lentype = np.zeros(6 * self.totngroups, 'i')
        self.masstype = np.zeros(6 * self.totngroups, float_type)
        self.ids = np.zeros(self.totnids, 'i')

        count = 0
        countids = 0

        for i in range(self.numfiles):

            if self.numfiles > 1:
                self.foffile = self.fofbase + '.' + repr(i)
            else:
                self.foffile = self.fofbase

            self.readheader()

            f = open(self.foffile, mode = 'rb')

            blksize = struct.unpack('i', f.read(4))[0]

            self.ngroups = struct.unpack('i', f.read(4))[0]

            f.seek(blksize + 12)

            if self.ngroups != 0:

               self.len[count: count + self.ngroups] = np.fromfile(f, dtype = 'i', count = self.ngroups)
               print self.len, np.sum(self.len)
   
               f.seek(8, 1)
   
               self.mass[count: count + self.ngroups] = np.fromfile(f, dtype = float_type, count = self.ngroups)
   
               f.seek(8, 1)
   
               self.pos[3 * count: 3 * count + 3 * self.ngroups] = np.fromfile(f, dtype = float_type, count = 3 * self.ngroups)
   
               f.seek(8, 1)
   
               self.cm[3 * count: 3 * count + 3 * self.ngroups] = np.fromfile(f, dtype = float_type, count = 3 * self.ngroups)
   
               f.seek(8, 1)
   
               self.vel[3 * count: 3 * count + 3 * self.ngroups] = np.fromfile(f, dtype = float_type, count = 3 * self.ngroups)
   
               f.seek(8, 1)
   
               self.lentype[6 * count: 6 * count + 6 * self.ngroups] = np.fromfile(f, dtype = 'i', count = 6 * self.ngroups)
   
               f.seek(8, 1)
   
               self.masstype[6 * count: 6 * count + 6 * self.ngroups] = np.fromfile(f, dtype = float_type, count = 6 * self.ngroups)
   
               f.seek(8, 1)

#               print self.foffile, self.len, self.mass, self.pos, self.cm, self.vel, self.lentype, self.masstype, self.nids, countids

            if self.nids != 0:

               self.ids[countids: countids + self.nids] = np.fromfile(f, dtype = 'i', count = self.nids)

            f.close()

            count += self.ngroups
            countids += self.nids
