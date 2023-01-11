import glob
import netCDF4 as nc
import numpy as np
import configparser
from scipy import interpolate
from textwrap import wrap

Combined=False

config=configparser.ConfigParser(allow_no_value=True)
config.read('/external/settings/data.ini')

Combined = config['Impacts']['Combined']
            
if not Combined:
    raise SystemExit(0)

print("\n******************************************************")
print("\nCombining Fields Data\n")
print("******************************************************\n")

files=glob.glob("/external/input/*.nc")
files=sorted(files)

c="/external/netcdf/out.jpg"
count=0
xmin=9e99
ymin=9e99
zmin=9e99
xmax=-9e99
ymax=-9e99
zmax=-9e99
xdiff=9e99
ydiff=9e99

inz = []

totalplts=sum(a.count("fields") for a in files)
print("------------------------------------------------------\n")

print("Creating new grid from x and y")

for a in files:

    if "fields" in a and "global" not in a:

        ncfile     = nc.Dataset(a)
    
        inx = ncfile.variables['x'][:]
        iny = ncfile.variables['y'][:]
        intime = ncfile.variables['time'][:]
        inz = ncfile.variables['C'][:]
                
        if inx.max() > xmax:
            xmax = inx.max()
        if iny.max() > ymax:
            ymax = iny.max()
        if inz.max() > zmax:
            zmax = inz.max()
        if inx.min() < xmin:
            xmin = inx.min()
        if iny.min() < ymin:
            ymin = iny.min()
        if inz.min() < zmin:
            zmin = inz.min()

        xsort=np.sort(inx)
        xd=min(n2 - n1 for n1, n2 in zip(xsort, xsort[1:]))
        if xd < xdiff:
            xdiff = xd
        ysort=np.sort(iny)
        yd=min(n2 - n1 for n1, n2 in zip(ysort, ysort[1:]))
        if yd < ydiff:
            ydiff = yd
            
            xmax=xmax+((xmax-xmin)%xdiff)
            ymax=ymax+((ymax-ymin)%ydiff)
            
        tgrid=ncfile.variables['time'][:]
              
xgrid=np.arange(xmin, xmax, xdiff)
ygrid=np.arange(ymin, ymax, ydiff)
zgrid=np.zeros([inz.shape[0],xgrid.shape[0], ygrid.shape[0]])
zgrid1=np.zeros([inz.shape[0],xgrid.shape[0], ygrid.shape[0]])
zgrid2=np.zeros([inz.shape[0],xgrid.shape[0], ygrid.shape[0]])
xgrid2, ygrid2 = np.meshgrid(xgrid,ygrid)
xgrid2=np.swapaxes(xgrid2,0,1)
ygrid2=np.swapaxes(ygrid2,0,1)

print("")

for a in files:
    if "fields" in a and "global" not in a:
    
        ncfile     = nc.Dataset(a)

        print("Interpolating {} to new grid".format(a[16:]))
        inx = ncfile.variables['x'][:]
        iny = ncfile.variables['y'][:]
        inz = ncfile.variables['C'][:]

        inz[:,:,-1]=0.0
        inz[:,:,0]=0.0
        inz[:,-1,:]=0.0
        inz[:,0,:]=0.0

        xin, yin = np.meshgrid(inx,iny)
        xin=np.swapaxes(xin,0,1)
        yin=np.swapaxes(yin,0,1)

#       Nearest node code ------------------------------------------
       
##        for x in inx:
##           xloc=(np.abs(inx - x)).argmin()
##           idx = (np.abs(xgrid - x)).argmin()
##           for y in iny:
##               yloc=(np.abs(iny - y)).argmin()
##               idy = (np.abs(ygrid - y)).argmin()
##               zgrid[:,idx,idy]=zgrid[:,idx,idy]+inz[:,xloc,yloc]

#       ------------------------------------------------------------

        xin=np.reshape(xin,(1,np.product(xin.shape)))
        yin=np.reshape(yin,(1,np.product(yin.shape)))
        xin = np.array([xin for sublist in xin for xin in sublist])
        yin = np.array([yin for sublist in yin for yin in sublist])
#        
        for b in range(inz.shape[0]):

            zin=inz[b,:,:]
            zin=np.reshape(zin,(1,np.product(zin.shape)))
            zin = [zin for sublist in zin for zin in sublist]

#            f0 = interpolate.griddata((xin, yin), zin, (xgrid2, ygrid2), method='nearest')
#            f1 = interpolate.griddata((xin, yin), zin, (xgrid2, ygrid2), method='linear')
            f2 = interpolate.griddata((xin, yin), zin, (xgrid2, ygrid2), method='cubic')

#            f0[np.isnan(f0)] = 0.0
#            f1[np.isnan(f1)] = 0.0
            f2[np.isnan(f2)] = 0.0
#            zgrid[b,:,:]=zgrid[b,:,:]+f0
#            zgrid1[b,:,:]=zgrid1[b,:,:]+f1
            zgrid2[b,:,:]=zgrid2[b,:,:]+f2

#del f0, f1, f2, zin, xin, yin, inx, iny, inz, xgrid2, ygrid2
del f2, zin, xin, yin, inx, iny, inz, xgrid2, ygrid2


print("\nWriting new input file: combined.nc")
             
f = nc.Dataset('/external/input/combined.nc','w', format='NETCDF4')
f.createDimension('x', len(xgrid))
f.createDimension('y', len(ygrid))
f.createDimension('time', None)
x = f.createVariable('x', 'f4', 'x')
y = f.createVariable('y', 'f4', 'y')
time = f.createVariable('time', 'i4', 'time')
x[:] = xgrid
y[:] = ygrid
time[:] = tgrid
f.close()

#f = nc.Dataset('/external/netcdf/combined.nc','a', format='NETCDF4')
#C = f.createVariable('C - nearest', 'f4', ('time', 'x', 'y',))
#C[:,:,:] = zgrid
#f.close()
#
#del C, zgrid

#f = nc.Dataset('/external/netcdf/combined.nc','a', format='NETCDF4')
#C1 = f.createVariable('C - linear', 'f4', ('time', 'x', 'y',))
#C1[:,:,:] = zgrid1
#f.close()
#
#del C1, zgrid1

f = nc.Dataset('/external/input/combined.nc','a', format='NETCDF4')
C2 = f.createVariable('C - cubic', 'f4', ('time', 'x', 'y',))
C2[:,:,:] = zgrid2
f.close()

print("\n------------------------------------------------------")
print("\nTada..\n")

