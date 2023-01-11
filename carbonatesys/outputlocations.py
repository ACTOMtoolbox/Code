import glob
import netCDF4 as nc
import numpy as np
from textwrap import wrap
from pyproj import Proj

p = Proj(proj='utm',zone=31,ellps='WGS84', preserve_units=False) # use kwargs

def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

print("\n******************************************************")
print("\nCombining Fields Data\n")
print("******************************************************\n")

files=glob.glob("/external/input/*.nc")
files=sorted(files)

count=0
xmin=9e99
ymin=9e99
zmin=9e99
xmax=-9e99
ymax=-9e99
zmax=-9e99
xdiff=9e99
ydiff=9e99

totalplts=sum(a.count("fields") for a in files)
print("------------------------------------------------------\n")

print("Creating new grid from x and y")

for a in files:

    if a[17:23] == "fields":

        ncfile     = nc.Dataset(a)
    
        inx = ncfile.variables['x'][:]
        iny = ncfile.variables['y'][:]


        x = np.mean(inx)
        y = np.mean(iny) 
        x,y = p(x,y,inverse=True)

        print(y,x)

print("\n------------------------------------------------------")
print("\nTada..\n")

