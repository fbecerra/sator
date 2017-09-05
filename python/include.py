
import imp
#import Image
from itertools import cycle
import numpy as np
import os
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import struct
import sys
import time

from const import *

min_real_number = -sys.float_info.max
max_real_number = sys.float_info.max

#from mpi4py import MPI
#
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
rank = 0
rank_snap = 0
#size = comm.Get_size()
size = 1
size_snap = 1

