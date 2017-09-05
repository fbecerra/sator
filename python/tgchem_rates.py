
from all import *

figheight = 16.
screen_ratio = 1.
fontsize = 10

nplots = 4

infile = '../src/data/tgchem_rates.dat'

mp.rcParams['font.size'] = fontsize
mp.rc('legend',**{'fontsize': 0.7 * fontsize})

figsize = get_figsize(screen_ratio, figheight)

fig = plt.figure(figsize = figsize)

lines = ["-", "--", "-.", ":"]
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

f = open(infile, mode = 'rb')

h2mode = struct.unpack('i', f.read(4))[0]

temp_min = struct.unpack('d', f.read(8))[0]
temp_max = struct.unpack('d', f.read(8))[0]

num_temp = struct.unpack('i', f.read(4))[0]
num_eq = struct.unpack('i', f.read(4))[0]
num_chem = struct.unpack('i', f.read(4))[0]
num_cool = struct.unpack('i', f.read(4))[0]

temp = np.log10(np.fromfile(f, dtype = 'd', count = num_temp))
eq = np.fromfile(f, dtype = 'd', count = num_temp * num_eq)
chem = np.fromfile(f, dtype = 'd', count = num_temp * num_chem)
cool = np.fromfile(f, dtype = 'd', count = num_temp * num_cool)

if h2mode:
    xval = np.log10(np.fromfile(f, dtype = 'd', count = num_temp))
    fesc_fit = np.log10(np.fromfile(f, dtype = 'd', count = num_temp))
    fesc_sob = np.log10(np.fromfile(f, dtype = 'd', count = num_temp))
else:
    nplots = 3

eq[eq > 0] = np.log10(eq[eq > 0])
chem[chem > 0] = np.log10(chem[chem > 0])
cool[cool > 0] = np.log10(cool[cool > 0])
#fesc[fesc > 0] = np.log10(fesc[fesc > 0])

eq = np.ma.masked_where(eq == 0, eq)
chem = np.ma.masked_where(chem == 0, chem)
cool = np.ma.masked_where(cool == 0, cool)
#fesc = np.ma.masked_where(fesc == 0, fesc)

for i in np.arange(nplots):

    ax = fig.add_subplot(2, (nplots + 1) / 2, i + 1)

    linecycler = cycle(lines)
    colorcycler = cycle(colors)

    count = 0

    if i == 0:
        min_val = -20.
        max_val = min_real_number
        for j in np.arange(num_eq):
            rate = eq[j * num_temp : (j + 1) * num_temp]
            if np.count_nonzero(rate) > 1:
                line = next(linecycler)
                color = next(colorcycler)
                ax.plot(temp, rate, line, color = color, label = str(count))
                max_val = max(np.max(rate), max_val)
                count += 1
        else:
            print 'Eq rate', j, 'has no non-zero entries!'
    if i == 1:
        min_val = -40.
        max_val = 0.
        for j in np.arange(num_chem):
            rate = chem[j * num_temp : (j + 1) * num_temp]
            if np.count_nonzero(rate) > 1:
                line = next(linecycler)
                color = next(colorcycler)
                ax.plot(temp, rate, line, color = color, label = str(count))
                count += 1
        else:
            print 'Chem rate', j, 'has no non-zero entries!'
    if i == 2:
        min_val = -40.
        max_val = min_real_number
        for j in np.arange(num_cool):
            rate = cool[j * num_temp : (j + 1) * num_temp]
            if np.count_nonzero(rate) > 1:
                line = next(linecycler)
                color = next(colorcycler)
                ax.plot(temp, rate, line, color = color, label = str(count))
                max_val = max(np.max(rate), max_val)
                count += 1
        else:
            print 'Cool rate', j, 'has no non-zero entries!'
    if i == 3:
        min_val = -max_real_number
        max_val = min_real_number
        for j in np.arange(num_h2):
            if j == 0:
                rate = fesc_fit
            if j == 1:
                rate = fesc_sob
            if np.count_nonzero(rate) > 1:
                line = next(linecycler)
                color = next(colorcycler)
                ax.plot(xval, rate, line, color = color, label = str(count))
                min_val = max(np.min(rate), min_val)
                max_val = max(np.max(rate), max_val)
                count += 1
            else:
                print 'fesc', j, 'has no non-zero entries!'

    ax.set_ylim(min_val, max_val)

    ax.legend(bbox_to_anchor = [1.13, 1], labelspacing = 0.2, handlelength = 3)

ax.set_xlabel(r'${\rm log}\,T\,[{\rm K}]$')
ax.set_ylabel(r'${\rm log}\,\Lambda\,[{\rm erg}\,{\rm s}^{-1}\,{\rm cm}^{-3}]$')

plt.savefig('data/tgchem_rates.pdf')

print 'Done!'
