from all import *

rank            = 0
size            = 1
filesystem      = 1         # 0:local, 1:DM
verbose         = 1
base            = 'ah1dm'
initial_file    = 1
final_file      = 61
threshold_mass  = 1e8       # in Msun

for num_fof in range(initial_file, final_file+1):
	g = readfof(filesystem, base, num_fof)
	g.readheader()
	g.readgroups()
        if len(g.mass) <= 0:
           continue
        if g.mass.max()/g.hubbleparam > threshold_mass*solar_mass/unit_mass:
		print "Halo found at redshift z=", g.redshift, " and time t=", g.time*unit_time/sec_per_year, "yr in snapshot", num_fof
		break

if (num_fof == final_file) and (g.mass.max()/g.hubbleparam < threshold_mass*solar_mass/unit_mass):
	print "Halo not found :("
	sys.exit()

##### Halo properties #####
halo_mass  = g.mass.max()
mmg        = g.mass == halo_mass
g.pos      = np.array(g.pos).reshape((len(g.pos)/3.,-1))
halo_pos   = g.pos[mmg,:]/g.hubbleparam/(1+g.redshift)		#[kpc] and comoving coordinates
g.vel      = np.array(g.vel).reshape((len(g.vel)/3.,-1))
halo_vel   = g.vel[mmg,:]*(1+g.redshift)		#[km s-1] and comoving coordinates

##### Output #####
print "\t***Halo properties***"
print "Position\t\t:", g.pos[mmg,:][0], " in code units and ", halo_pos[0], " in kpc"
print "Velocity\t\t:", g.vel[mmg,:][0], " in code units and ", halo_vel[0], " in km*s-1"
print "Mass\t\t\t:", halo_mass, " in code units and ", halo_mass*unit_mass/solar_mass/g.hubbleparam, " in Msun"

GroupOffset = np.zeros(g.totngroups, int)
for i in range (1,g.totngroups):
	GroupOffset[i] = GroupOffset[i-1] + g.len[i-1]
offset = GroupOffset[mmg][0]
left   = g.len[mmg][0]


##### Reading snapshot #####
s = readsnap(filesystem, base, num_fof, comoving_units=0, verbose = verbose)
s.readheader(verbose)
s.readblocks([0, 1, 2, 3, 4, 5], ['pos', 'vel', 'id', 'mass'])
		# pos in [kpc], vel in [km s-1], mass in [g]. All of them in comoving coordinates c_c=0


##### Angular momentum #####
idxs      = np.argsort(s.id)
sortedsid = s.id[idxs]
idxg      = np.searchsorted(sortedsid, g.ids[offset:offset+left])
id_index  = idxs[idxg]
angmom_x  = np.sum(s.mass[id_index]*((s.y[id_index]-halo_pos[0][1])*(s.vz[id_index]-halo_vel[0][2])-(s.z[id_index]-halo_pos[0][2])*(s.vy[id_index]-halo_vel[0][1])))
angmom_y  = np.sum(s.mass[id_index]*((s.z[id_index]-halo_pos[0][2])*(s.vx[id_index]-halo_vel[0][0])-(s.x[id_index]-halo_pos[0][0])*(s.vz[id_index]-halo_vel[0][2])))
angmom_z  = np.sum(s.mass[id_index]*((s.x[id_index]-halo_pos[0][0])*(s.vy[id_index]-halo_vel[0][1])-(s.y[id_index]-halo_pos[0][1])*(s.vx[id_index]-halo_vel[0][0])))
J         = np.sqrt(angmom_x**2+angmom_y**2+angmom_z**2)*solar_mass*unit_length*unit_velocity   # [h-2 erg]


##### Halo properties, see Barkana & Loeb (2001) #####
omega_k   = 0
omega_mz  = s.omega_m*(1+s.redshift)**3/(s.omega_m*(1+s.redshift)**3+s.omega_lambda+omega_k*(1+s.redshift)**2)
d         = omega_mz - 1
delta_c   = 18*np.pi**2+82*d-39*d**2
r_vir     = 0.784*(halo_mass*unit_mass/solar_mass/10**8)**(1./3.)*((s.omega_m/omega_mz)*delta_c/(18*np.pi**2))**(-1./3.)*((1+s.redshift)/10)**(-1)
v_c       = 23.4*(halo_mass*unit_mass/solar_mass/10**8)**(1./3.)*((s.omega_m/omega_mz)*delta_c/(18*np.pi**2))**(1./6.)*((1+s.redshift)/10)**(1./2.)
T_vir     = 1.98e4*(mu_prim/0.6)*(halo_mass*unit_mass/solar_mass/10**8)**(2./3.)*((s.omega_m/omega_mz)*delta_c/(18*np.pi**2))**(1./3.)*((1+s.redshift)/10)
Eb        = 5.45e53*(halo_mass*unit_mass/solar_mass/10**8)**(5./3.)*((s.omega_m/omega_mz)*delta_c/(18*np.pi**2))**(1./3.)*((1+s.redshift)/10)
spinparam = J*s.hubbleparam**2*np.sqrt(Eb)/gravity/(halo_mass*unit_mass)**(5./2.)

##### Output #####
print "\t***More halo properties calculated!***"
print "Virial radius\t\t: ", r_vir, "[h-1 kpc], ", r_vir/s.hubbleparam, "[kpc]"
print "Circular velocity\t: ", v_c, "[km s-1]"
print "Virial Temperature\t: ", T_vir, "[K]"
print "Binding energy\t\t: ", Eb, "[h-1 erg], ", Eb/s.hubbleparam, "[erg]"
print "Spin parameter\t\t: ", spinparam

