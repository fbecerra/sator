
from all import *

args = len(sys.argv)

if args < 3:
    print 'Not enough arguments specified! Exiting...'
    comm.Abort()

base = sys.argv[1]
snapnum = int(sys.argv[2])

fof = readfof(0, base, snapnum)

fof.readgroups()

tgroup = 0

print fof.mass[tgroup] * unit_mass / solar_mass / fof.hubbleparam
print fof.pos[3 * tgroup: 3 * tgroup + 3]
