
from all import *

def get_str(field, count):
  if field == 'cool':
    if count == 0:
      label = 'HydroHeat'
    if count == 1:
      label = 'H2 rovib'
    if count == 2:
      label = 'H2 CIE'
    if count == 3:
      label = 'HI Ly-a'
    if count == 4:
      label = 'H- cont.'
    if count == 5:
      label = 'IC'
  elif field == 'ab':
    if count == 0:
      label = 'H-'
    if count == 1:
      label = 'H2'
    if count == 2:
      label = 'H+'
  return label

num_cooling_omukai = 3 

figheight = 24
#figheight = 12
fontsize = 12

#nplots = 7
nplots = 3

abs_min = 1e-306

infile = 'data/tgchem.dat'

mp.rcParams['font.size'] = fontsize
mp.rc('legend',**{'fontsize': 0.7 * fontsize})

#figsize = [nplots * figheight / 4 * cm_to_inch, figheight * cm_to_inch]
figsize = [ figheight  * cm_to_inch, figheight * cm_to_inch]

fig = plt.figure(figsize = figsize)

lines = ["-", "--", "-.", ":"]
nlines = len(lines)
linecycler = cycle(lines)

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorcycler = cycle(colors)

f = open(infile, mode = 'rb')

num_iter = struct.unpack('i', f.read(4))[0]

num_abundances = struct.unpack('i', f.read(4))[0]
num_rates = struct.unpack('i', f.read(4))[0]
num_tot_rates = num_abundances * 2 * num_rates
num_chemical = num_abundances * num_rates
num_cooling = 1 + struct.unpack('i', f.read(4))[0]
num_entries = 3 + num_abundances + num_tot_rates + num_chemical + num_cooling

nh_min = np.log10(struct.unpack('d', f.read(8))[0])
nh_max = np.log10(struct.unpack('d', f.read(8))[0])

abstol = np.log10(np.fromfile(f, dtype = 'd', count = num_abundances))
min_abstol = max_real_number

for j in np.arange(num_abundances):
    min_abstol = min(abstol[j], min_abstol)

data = np.fromfile(f, dtype = 'd', count = num_iter * num_entries)

data[data < 0] = -data[data < 0]
data[data < abs_min] = abs_min

idx = 0
nh = np.log10(data[ : : num_entries])
idx += 1
temp = np.log10(data[idx : : num_entries])
idx += 1
gamma = data[idx : : num_entries]
idx += 1

rates = []
abundances = []
chemical = []
cooling = []

for i in np.arange(num_abundances):
    abundances.append(np.log10(data[idx : : num_entries]))
    idx += 1

for i in np.arange(num_tot_rates):
    rates.append(np.log10(data[idx : : num_entries]))
    idx += 1

for i in np.arange(num_chemical):
    chemical.append(np.log10(data[idx : : num_entries]))
    idx += 1

for i in np.arange(num_cooling):
    cooling.append(np.log10(data[idx : : num_entries]))
    idx += 1

#nh_min = 2
#nh_max = 6

for i in np.arange(nplots):

    ax = fig.add_subplot(2, (nplots + 1) / 2, i + 1)

    if i == 0:
        ax.plot(nh, temp)
        x0, y0 = 20, 3.5
        x = np.array([x0, max(temp)])
        y = x**2./3. + y0 - x0**2./3.
        ax.plot(x, y, 'r--') 
        min_val = np.min(temp)
        max_val = np.max(temp)
        min_val, max_val = 3.5, 4.1
        ylabel = get_label('temp')
#    if i == 1:
#        ax.plot(nh, gamma)
#        min_val = 1
#        max_val = 2
#        ylabel = get_label('gamma')
    if i == 1:
        min_val = max_real_number
        max_val = min_real_number
        linecycler = cycle(lines)
        colorcycler = cycle(colors)
        for j in np.arange(num_abundances):
            line = next(linecycler)
            color = next(colorcycler)
            ax.plot(nh, abundances[j], line, color = color, label = get_str('ab', j))
            min_val = min(np.min(abundances[j]), min_val)
            max_val = max(np.max(abundances[j]), max_val)
        min_val = max(min_val, min_abstol)
        ylabel = r'$y_{\rm X}$'
        ax.legend(labelspacing = 0.2, handlelength = 3)
#    if i > 2 and i < 6:
#        if i == 3: idx = 0
#        min_val = -40
#        max_val = min_real_number
#        linecycler = cycle(lines)
#        colorcycler = cycle(colors)
#        count = 0
#        for j in np.arange(2 * num_rates):
#            rate = rates[idx + j]
#            if np.count_nonzero(rate != np.log10(abs_min)) > 0:
#                line = next(linecycler)
#                color = next(colorcycler)
#                ax.plot(nh, rate, line, color = color, label = str(count))
#                max_val = max(np.max(rate), max_val)
#                count += 1
#        idx += 2 * num_rates
#        ylabel = r'${\rm rate}\,[{\rm s}^{-1}\,{\rm cm}^3]$'
#    if i == 6:
    if i == 2:
        min_val = -10
        max_val = 10
        linecycler = cycle(lines)
        colorcycler = cycle(colors)
#        count = 0
#        for j in np.arange(num_chemical):
#            chem = chemical[j]
#            if np.count_nonzero(chem != np.log10(abs_min)) > 0:
#                line = next(linecycler)
#                color = next(colorcycler)
#                ax.plot(nh, chem, line, color = color, label = str(count))
#                max_val = max(np.max(chem), max_val)
#                count += 1
        count = 0
        for j in np.arange(num_cooling):
            cool = cooling[j]
            print cool
            if np.count_nonzero(cool != np.log10(abs_min)) > 0:
                line = next(linecycler)
#                color = next(colorcycler)
                color = 'k'
                cool = np.log10(10**cool * 10**nh / protonmass)
                thickness = count / nlines + 0.5
                if j in [2, 3, 4]:
                  print 'yes'
                  thickness = 1.5
                  ax.plot(nh, cool, line, color = color, label = get_str('cool', j), linewidth = thickness)
                #max_val = max(np.max(cool), max_val)
                count += 1
        count = 2
        linecycler = cycle(lines)
        line = next(linecycler)
        line = next(linecycler)
        for j in np.arange(num_cooling_omukai):
            if j == 0: 
               data_file = '/n/home02/fbecerra/arepo/sator/data/Omukai01/h2lines.dat'
            elif j == 1:
               data_file = '/n/home02/fbecerra/arepo/sator/data/Omukai01/Hlines.dat'
            elif j == 2:
               data_file = '/n/home02/fbecerra/arepo/sator/data/Omukai01/continuum.dat'
            f=open(data_file, 'r')
            nh = []
            cool = []
            for line in f:
               line=line.split()
               nh.append(float(line[0]))
               cool.append(float(line[1]))
            nh = np.array(nh)
            cool = np.array(cool)
            line = next(linecycler)
            color = 'r'
            thickness = count / nlines + 0.5
            thickness = 1.5
            ax.plot(nh, cool, line, color = color, label = get_str('cool', count), linewidth = thickness)
            count += 1

        ylabel = r'$\Lambda/\rho\,[{\rm erg}\,{\rm s}^{-1}\,{\rm g}^{-1}]$'
        #ax.legend(bbox_to_anchor = [1.3, 1], labelspacing = 0.2, handlelength = 3)
        ax.legend(labelspacing = 0.2, handlelength = 3)

    if max_val == min_real_number:
        max_val = 20.

    ax.set_xlim(nh_min, nh_max)
    ax.set_ylim(min_val, max_val)

    ax.set_xlabel(r'${\rm log}\,n_{\rm H}\,[{\rm cm}^{-3}]$')
    ax.set_ylabel(ylabel)

fig.tight_layout()

plt.savefig('data/out.pdf')

print 'Done!'
