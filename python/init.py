
from all import *

s = readsnap('h4', 221, verbose = 1)

s.readblocks([0, 1, 2, 3, 4, 5], ['pos', 'vel', 'id', 'mass', 'u', 'rho', 'vol', 'chem', 'gamma'])

s.mass = s.mass[:s.npart_local[0]]
