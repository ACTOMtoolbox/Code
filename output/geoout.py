import os
import glob
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from textwrap import wrap
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
plt.style.use('seaborn-whitegrid')

files=glob.glob("input/geo/*.nc")

b=".png"

my_dpi=300

if not os.path.isdir('/tmp/geo'):
    os.system('mkdir /tmp/geo')

for a in files:
    c="/tmp/geo/"+a[10:-26]+b
    
    ncfile     = nc.Dataset(a)
    outputList = []

    for v in ncfile.variables:
        outputList += [str(v)]
    x = ncfile.variables[outputList[0]][:]
    y = ncfile.variables[outputList[1]][:]
    mag = ncfile.variables[outputList[2]][:]

    plt.imshow(mag, extent=(x.min(), x.max(), y.max(), y.min()),
        interpolation='nearest')
    title="Synthetic Probability Map from "+a[10:-26]
    titlea='\n'.join(wrap(title,60))
    titlea=titlea+'\n '
    plt.title(titlea)
    plt.xlabel("X")
    plt.ylabel("Y")
    if x.max() > 1000:
      plt.ticklabel_format(axis="X", style="sci", scilimits=(0,0))
    if y.max() > 1000:
      plt.ticklabel_format(axis="Y", style="sci", scilimits=(0,0))
    plt.grid(b=None)
    plt.savefig(c, dpi=my_dpi, facecolor='w', edgecolor='w',
                orientation='portrait', format=None,
                transparent=False, bbox_inches='tight',  pad_inches=0.1,
                metadata=None)
    plt.close()
