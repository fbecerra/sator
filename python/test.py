from all import *

rank            = 0
size            = 1
filesystem      = 1         # 0:local, 1:DM
verbose         = 1
base            = 'nahw1rs1'
initial_file    = 5
final_file      = 19
threshold_mass  = 1e8       # in Msun

time = np.array([])
mass = np.array([])

for snap in range(initial_file, final_file+1):
	g = readsnap(filesystem, base, snap)
	g.readblocks([0, 1, 2, 3, 4, 5], ['pos', 'vel', 'id', 'mass'])

	idx = np.where(g.id == 1025028687)[0]
        if len(idx) > 0:
		length = len(g.mass)
		time = np.append(time, g.time * unit_time / sec_per_year)
	        mass = np.append(mass, g.mass[idx])
	#print snap, idx, g.x[idx], g.y[idx], g.z[idx], g.mass[idx]
        print g.time * unit_time / sec_per_year, g.mass[idx]

print time, mass

