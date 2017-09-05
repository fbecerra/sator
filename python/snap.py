
from include import *
from utils import *

class readsnap:

    def __init__(self, filesystem = 0, base = '/', snapnum = 0, flag_dist = 0, flag_double = 1, length_unit = 'kpc', comoving_units = 0, verbose = 1):

        self.base = base
        self.flag_dist = flag_dist
        self.flag_double = flag_double

        self.path, self.snapnum, self.snapfile = setup_snapfile(filesystem, snapnum, base)

        self.snapbase = self.snapfile

        if not os.path.isfile(self.snapfile):

            if self.snapnum != 'ic':
                self.snapbase = self.path + '/' + self.base + '/' + 'snapdir' + '_' + self.snapnum + '/' + self.base + '_' + self.snapnum

            self.snapfile = self.snapbase + '.0'

            if not os.path.isfile(self.snapfile):
                print 'Could not find snapshot file! Exiting...'
                comm.Abort()

        if self.flag_dist == 0:

            size_snap = 1
            rank_snap = 0

        else:

            size_snap = size
            rank_snap = rank

        self.length_unit = length_unit
        self.comoving_units = comoving_units

        self.readheader(verbose)

    def readheader(self, verbose):

        f = open(self.snapfile, mode = 'rb')

        f.seek(28)

        self.masstable = np.array(struct.unpack(repr(nsnaptypes) + 'd', f.read(48)))
        self.time = struct.unpack('d', f.read(8))[0]
        self.redshift = struct.unpack('d', f.read(8))[0]
        self.flag_sfr = struct.unpack('i', f.read(4))[0]
        self.flag_feedback = struct.unpack('i', f.read(4))[0]
        self.npartall = np.array(struct.unpack(repr(nsnaptypes) + 'i', f.read(24)))
        self.flag_cooling = struct.unpack('i', f.read(4))[0]
        self.numfiles = struct.unpack('i', f.read(4))[0]
        self.boxsize = struct.unpack('d', f.read(8))[0]
        self.omega_m = struct.unpack('d', f.read(8))[0]
        self.omega_lambda = struct.unpack('d', f.read(8))[0]
        self.hubbleparam = struct.unpack('d', f.read(8))[0]

        f.close()

        self.npart = np.zeros((self.numfiles, nsnaptypes), dtype = 'i')

        for i in range(self.numfiles):

            if self.numfiles > 1:

                self.snapfile = self.snapbase + '.' + repr(i)

            f = open(self.snapfile, mode = 'rb')

            f.seek(4)

            self.npart[i, :] = np.array(struct.unpack(repr(nsnaptypes) + 'i', f.read(24)))

            f.close()

        if self.length_unit == 'kpc':
            fac = 1
        elif self.length_unit == 'pc':
            fac = 1e3
        elif self.length_unit == 'au':
            fac = unit_length / astronomical_unit
        else:
            print 'Length unit not implemented! Exiting...'
            comm.Abort()

        fac /= self.hubbleparam

        self.boxsize *= fac

#        if rank == 0 and verbose == 1:
#            print 'mass = ', self.masstable
#            print 'time = ', self.time
#            print 'redshift = ', self.redshift
#            print 'flag_sfr = ', self.flag_sfr
#            print 'flag_feedback = ', self.flag_feedback
#            print 'npartall = ', self.npartall
#            print 'flag_cooling = ', self.flag_cooling
#            print 'numfiles = ', self.numfiles
#            print 'boxsize = ', self.boxsize
#            print 'omega_m = ', self.omega_m
#            print 'omega_lambda = ', self.omega_lambda
#            print 'hubbleparam = ', self.hubbleparam

        self.npart_des = np.zeros((size_snap, nsnaptypes), dtype = 'i')

        self.npart_des[range(size_snap), :] = self.npartall / size_snap

        for j in range(nsnaptypes):

            self.npart_des[range(size_snap) < self.npartall[j] % size_snap, j] += 1

        self.npart_local = self.npart_des[rank_snap]

        self.snap_list = [[[] for j in range(nsnaptypes)] for i in range(size_snap)]
        self.list_start = [[[] for j in range(nsnaptypes)] for i in range(size_snap)]
        self.list_end = [[[] for j in range(nsnaptypes)] for i in range(size_snap)]

        for j in range(nsnaptypes):

            partcount = 0

            rank_idx = 0

            for i in range(self.numfiles):

                if self.npart[i, j] > 0:

                    partcount += self.npart[i, j]

                    self.snap_list[rank_idx][j].append(i)
                    self.list_start[rank_idx][j].append(0)

                    while partcount > self.npart_des[rank_idx, j]:

                        partcount -= self.npart_des[rank_idx, j]

                        self.list_end[rank_idx][j].append(self.npart[i, j] - partcount - 1)

                        rank_idx += 1

                        self.snap_list[rank_idx][j].append(i)

                        self.list_start[rank_idx][j].append(self.npart[i, j] - partcount)

                    self.list_end[rank_idx][j].append(self.npart[i, j] - 1)

        count_npart = np.zeros(nsnaptypes, dtype = 'i')
        count_des = np.zeros((size_snap, nsnaptypes), dtype = 'i')

        for i in range(size_snap):

            for j in range(nsnaptypes):

                if self.snap_list[i][j]:

                    for k in range(len(self.snap_list[i][j])):

                        count_npart[j] += self.list_end[i][j][k] - self.list_start[i][j][k] + 1
                        count_des[i, j] += self.list_end[i][j][k] - self.list_start[i][j][k] + 1

        if np.any(count_npart != self.npartall) or np.any(count_des != self.npart_des):

            print 'Inconsistent particle numbers! Exiting...'

            comm.Abort()

    def get_block_info(self, type_in, snap_idx, target_block, list_start):

        self.offset = 264

        for block in snapblocks:

            self.offset += 4

            # Number of dimensions

            if block == 'pos' or block == 'vel' or block == 'gravacc' or block == 'gradp':
                self.ndims = 3
            elif block == 'chem':
                self.ndims = 6
            elif block == 'rates':
                self.ndims = 4
            else:
                self.ndims = 1

            # Number of elements

            if block == target_block:
                type = type_in
            else:
                type = nsnaptypes

            if block == 'pos' or block == 'vel' or block == 'id' or block == 'gravacc':
                nelements = self.npart[snap_idx, 0: type].sum()
            elif block == 'mass':
                nelements = self.npart[snap_idx, self.masstable[0: type] == 0].sum()
            else:
                if block == target_block:
                    nelements = 0
                else:
                    nelements = self.npart[snap_idx, 0]

            if block == target_block:
                nelements += list_start

            # Type and bytes

            if block == 'id' or block == 'allowref':
                self.dtype = 'i'
                bytes_per_element = 4
            else:
                if self.flag_double == 1:
                    self.dtype = 'd'
                    bytes_per_element = 8
                else:
                    self.dtype = 'f'
                    bytes_per_element = 4

            # Offset in bytes

            self.offset += self.ndims * nelements * bytes_per_element

            if block != target_block:
                self.offset += 4

            # Number of elements

            if block == target_block:

                return

    def set_block_value(self, block, data):

        if block == 'pos':

            if self.length_unit == 'kpc':
                fac = 1
            elif self.length_unit == 'pc':
                fac = 1e3
            elif self.length_unit == 'au':
                fac = unit_length / astronomical_unit
            else:
                print 'Length unit not implemented! Exiting...'
                comm.Abort()

            if self.comoving_units == 0:
                fac /= (1 + self.redshift)

            fac /= self.hubbleparam

            self.x = fac * data[::3]
            self.y = fac * data[1::3]
            self.z = fac * data[2::3]

        elif block == 'vel':

            fac = np.sqrt(1 + self.redshift)**-1

            self.vx = fac * data[::3]
            self.vy = fac * data[1::3]
            self.vz = fac * data[2::3]

        elif block == 'id':

            self.id = data

        elif block == 'mass':

            fac = unit_mass / solar_mass / self.hubbleparam

            self.mass = fac * data

        elif block == 'u':

            fac = unit_energy / unit_mass

            self.u = fac * data

            fac = self.mu * (self.gamma - 1) * protonmass / boltzmann

            self.temp = fac * self.u

        elif block == 'rho':

            fac = self.hubbleparam**2 * (1 + self.redshift)**3 * unit_mass / unit_length**3

            self.rho = fac * data

            fac = hydrogen_massfrac / protonmass

            self.nh = fac * self.rho

        elif block == 'vol':

            if self.length_unit == 'kpc':
                fac = 1
            elif self.length_unit == 'pc':
                fac = 1e9
            elif self.length_unit == 'au':
                fac = (unit_length / astronomical_unit)**3

            if self.comoving_units == 0:
                fac /= (1 + self.redshift)**3

            fac /= self.hubbleparam**3

            self.vol = fac * data

            fac = (3. / 4 / np.pi)**(1. / 3)

            self.hsml = fac * self.vol**(1. / 3)

        elif block == 'delaunay':

            self.delaunay = data

        elif block == 'gravacc':

            self.gravaccx = data[::3]
            self.gravaccy = data[1::3]
            self.gravaccz = data[2::3]

        elif block == 'gradp':

            self.gradpx = data[::3]
            self.gradpy = data[1::3]
            self.gradpz = data[2::3]

        elif block == 'chem':

            self.h2i = data[::6]
            self.hii = data[1::6]
            self.dii = data[2::6]
            self.hdi = data[3::6]
            self.heii = data[4::6]
            self.heiii = data[5::6]

            self.abe = self.hii + self.dii + self.heii + 2 * self.heiii

            idx = self.abe > 1 + abde + 2 * abhe
            
            self.abe[idx] = 1 + abde + 2 * abhe

            self.mu = (1 + 4 * abhe) / (1 + abhe - self.h2i + self.abe)

        elif block == 'gamma':

            self.gamma = data

        elif block == 'rates':

            self.pdvrate = data[::4]
            self.h2rate = data[1::4]
            self.cierate = data[2::4]
            self.chemrate = data[3::4]

        elif block == 'allowref':

            self.allowref = data

    def readblocks(self, types = snaptypes, blocks = snapblocks):

        self.gamma = np.array(gamma_adb)
        self.mu = np.array(mu_prim)

        if 'u' in blocks and ('gamma' in blocks or 'chem' in blocks):
            del blocks[blocks.index('u')]
            blocks.append('u')

        data_table = [[] for i in range(nsnapblocks)]
        data_offset = [0 for i in range(nsnapblocks)]

        for type in types:

            count = 0

            for i in self.snap_list[rank_snap][type]:

                for block in blocks:

                    if type == 0 or block == 'pos' or block == 'vel' or block == 'id' or block == 'mass':

                        self.get_block_info(type, i, block, self.list_start[rank_snap][type][count])

                        data_idx = snapblocks.index(block)

                        if len(data_table[data_idx]) == 0:

                            nelements = 0

                            for tmp_type in types:
                                if tmp_type == 0 or block == 'pos' or block == 'vel' or block == 'id' or block == 'mass':
                                    nelements += self.ndims * self.npart_local[tmp_type]

                            data_table[data_idx] = np.zeros(nelements, self.dtype)

                        data = data_table[data_idx]

                        nelements = self.ndims * (self.list_end[rank_snap][type][count] - self.list_start[rank_snap][type][count] + 1)

                        if block == 'mass' and self.masstable[type] != 0:

                            data[data_offset[data_idx]: data_offset[data_idx] + nelements] = self.masstable[type]

                        else:

                            if self.numfiles > 1:
                                self.snapfile = self.snapbase + '.' + repr(i)
                            else:
                                self.snapfile = self.snapbase

                            f = open(self.snapfile, mode = 'rb')

                            f.seek(self.offset)

                            data[data_offset[data_idx]: data_offset[data_idx] + nelements] = np.fromfile(f, dtype = self.dtype, count = nelements)

                            f.close()

                        data_offset[data_idx] += nelements

                count += 1

        for block in blocks:

            data_idx = snapblocks.index(block)

            if len(data_table[data_idx]) != 0:

                data = data_table[data_idx]

                self.set_block_value(block, data)
