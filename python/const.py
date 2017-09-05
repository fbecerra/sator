
pi = 3.1415927
gamma_adb = 5. / 3
hydrogen_massfrac = 0.76
abhe = (1 - hydrogen_massfrac) / 4 / hydrogen_massfrac
abde = 2.6e-5
mu_prim = 1 / (hydrogen_massfrac * (1 + abhe))
protonmass = 1.6726e-24
boltzmann = 1.3806e-16
electron_volt = 1.60219e-12
gravity = 6.672e-8
solar_mass = 1.989e33
hubble = 3.2407789e-18
sec_per_year = 3.155e7
astronomical_unit = 1.49598e13
solar_radius = 6.955e10

unit_length = 3.085678e21
unit_mass = 1.989e43
unit_velocity = 1e5
unit_time = unit_length / unit_velocity
unit_energy = unit_mass * unit_length**2 / unit_time**2

cm_to_inch = 0.393700787

nsnaptypes = 6
snapblocks = ['pos', 'vel', 'id', 'mass', 'u', 'rho', 'vol', 'delaunay', 'gravacc', 'gradp', 'chem', 'gamma', 'rates', 'allowref']

snaptypes = range(nsnaptypes)
nsnapblocks = len(snapblocks)
