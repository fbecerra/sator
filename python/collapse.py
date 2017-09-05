
from all import *

figheight = 6.5
screen_ratio = 8. / 6.5
fontsize = 7.

lbuf = 0.16
rbuf = 0.16
tbuf = 0.02
hspacing = 0.04

vspacing = hspacing * screen_ratio
pwidth = (1. - lbuf - rbuf - hspacing) / 2
pheight = pwidth * screen_ratio

mp.rcParams['font.size'] = fontsize

figsize = get_figsize(screen_ratio, figheight)

fig = plt.figure(figsize = figsize)

lines = ["-", "--", "-.", ":"]
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

infile = 'collapse.dat'

f = open(infile, mode = 'rb')

num_cases = struct.unpack('i', f.read(4))[0]
num_entries = struct.unpack('i', f.read(4))[0]

val = [[] for i in range(num_cases)]

for i in np.arange(num_cases):

    num_iter = struct.unpack('i', f.read(4))[0]

    data = np.fromfile(f, dtype = 'd', count = num_iter * num_entries)

    for j in np.arange(num_entries):
        val[i].append(data[j : : num_entries])

xlabel = get_label('nh')

xval_min = max_real_number
xval_max = min_real_number

for i in np.arange(num_entries - 1):

    pleft = lbuf + (i % 2) * (pwidth + hspacing)
    pbottom = 1 - tbuf - (i / 2 + 1) * pheight - i / 2 * vspacing

    axes = fig.add_axes([pleft, pbottom, pwidth, pheight])

    if i == 0:
        block = 'temp'

    if i == 1:
        block = 'abh2'

    if i == 2:
        block = 'geff'

    if i == 3:
        block = 'cool'

    linecycler = cycle(lines)
    colorcycler = cycle(colors)

    for j in np.arange(num_cases):

        line = next(linecycler)
        color = next(colorcycler)

        axes.plot(val[j][0], val[j][i + 1], line, color = color)

        xval_min = min(min(val[j][0]), xval_min)
        xval_max = max(max(val[j][0]), xval_max)

    if block == 'geff':
        axes.plot([xval_min, xval_max], [1.1, 1.1], 'grey', linestyle = '--')

    if block == 'cool':
        axes.plot([xval_min, xval_max], [0, 0], 'grey', linestyle = '--')

    ylabel = get_label(block)

    if(block == 'temp'):
        ylabel = r'$T\,[{\rm K}]$'

    if i < 2:
        axes.set_xticklabels([])

    if i == 1 or i == 3:
         axes.yaxis.tick_right()
         axes.yaxis.set_label_position('right')

    if i > 1:
        axes.set_xlabel(xlabel, fontsize = fontsize + 1)

    axes.set_ylabel(ylabel, fontsize = fontsize + 1)

    axes.set_xlim(xval_min, xval_max)

    if block == 'temp':
        axes.set_ylim(0, 2.2e3)

    if block == 'abh2':
        axes.set_ylim(-3.5, 0)

    if block == 'geff':
        axes.set_ylim(0.4, 1.4)

    if block == 'cool':
        axes.set_ylim(-0.15, 0.5)


    axes.xaxis.set_major_locator(plt.MaxNLocator(5))
    #axes.locator = plt.MaxNLocator(5 + 1)

    #axes.update_ticks()

#ax.set_xlabel(r'${\rm log}\,n_{\rm H}\,[{\rm cm}^{-3}]$', fontsize = labelsize)
#ax.set_ylabel(r'${\rm log}\,T\,[{\rm K}]$', fontsize = labelsize)

plt.savefig('collapse.pdf')

print 'Done!'
