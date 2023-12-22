import os
import glob
import netCDF4 as nc
import numpy as np
import configparser
import imageio
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
from matplotlib.transforms import Bbox
from textwrap import wrap
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib import ticker, cm
import warnings
warnings.filterwarnings("ignore")
                               
Animation='No'

print("\n******************************************************")
print("\nProcessing Impacts\n")
print("******************************************************\n")

files=glob.glob("/external/input/*.nc")
files=sorted(files)
my_dpi=300
count=0
limit=0.01
fps=24
file = open('temp.txt', 'a')
file.close()
 

config=configparser.ConfigParser(allow_no_value=True)
config.read('/external/settings/data.ini')

mode = config['General']['run'] # brine? 

Output = config['Impacts']['output']
Local = config['Impacts']['local']
if config.has_option('Impacts', 'animation'):
    Animation = config['Impacts']['animation']
Threshold = []

if config.has_option('General', 'threshold-ph'):
    th1=((config['General']['threshold-ph']).split(","))
    th2=[float(x) for x in th1] 
    Threshold.append(th2)
if config.has_option('RateOfChange', 'threshold-ph'):
    th1=((config['RateOfChange']['threshold-ph']).split(","))
    th2=[float(x) for x in th1] 
    Threshold.append(th2)
if config.has_option('CSEEP', 'threshold-ph'):
    th1=((config['CSEEP']['threshold-ph']).split(","))
    th2=[float(x) for x in th1] 
    Threshold.append(th2)

Thresholds = [val for sublist in Threshold for val in sublist]
    
if mode == 'Brine':
   Threshold=[1/i for i in Thresholds]
   Thresholds=Threshold
    
Thresholds.sort()

for a in files:

    if "global" not in a and 'cumulative' not in a:
      if Local=='Global':
        continue

    if "mass" in a or "point_tags" in a:
        continue    
       
#    if "probe" in a or "statistics" in a or "global" not in a:
#        continue         

    print("------------------------------------------------------\n")
    print("Opening File: {}\n".format(a[16:]))

    ncfile     = nc.Dataset(a)
    
    x = ncfile.variables['x'][:]
    y = ncfile.variables['y'][:]
    if 'time' in list(ncfile.variables):
        t = ncfile.variables['time'][:]
    
    dx= np.diff(x)
    np.all(dx[0] == dx)
    dx=np.append(dx,dx[-1])
    
    dy=np.diff(y)
    np.all(dy[0] == dy)
    dy=np.append(dy,dy[-1])
    if "t" in locals():
        dt=np.diff(t)
        np.all(t[0] == t)
        dt=np.append(dt,dt[-1])
    
    # add in height to make volumes when available --------------!
    
    for v in ncfile.variables:

        if v[0:3]=="DIC" and "DIC" in Output or v[0:2]=="pH" and "pH" in Output or mode == 'Brine':
            if mode == 'Brine':
               if v == "u" or v == "v" or v == "x" or v == "y" or \
                  v == "lon" or v == "lat" or v == "time" or v == "indx" or \
                  v == "x_grid" or v == "y_grid" or v == "x_pos" or \
                  v == "location_probability" or v == "x_source" or v == "y_source" or \
                  v == "lon_source" or v == "lat_source" or v == "lat_source" or \
                  v == "y_pos" or v == "dist" or v == "dC" or "delta" in v or  \
                  v == "location probabilities" or v == "source tag" or v == "dt" or \
                  v == "probe" or "var" in v:
                     continue
                     
            if v[0:13]=="pH from delta" or v[0:11]=="pH from var":
                continue
#            if a[16:]=="combined.nc" and v!="pH from C - cubic":
#               continue
            print("Processing impacts for variable: {}".format(v))
            c=a[:-3]+v
            c = c.replace("input", "output")
            maga = ncfile.variables[v][:]
            if "global" in a and mode == 'Brine':
               maga = np.sum(maga, axis=0)
            if v[0:2]=="pH":
               mag = abs(maga-maga.max())
            else:
               mag = abs(maga-maga.min())
            magmax=np.amax(mag, axis=0)
            maxchange = abs(maga.max()-maga.min())

            if "probe" in a:
                time = ncfile.variables['time'][:]
                indx = ncfile.variables['probe'][:]
                fig,ax = plt.subplots()

                for q in range(0,maga.shape[1]):
                    plt.plot(time,maga[:,q],label=indx[q])
                plt.grid(visible=None)
                fig.legend()
                plt.savefig(c, dpi=my_dpi, facecolor='w', edgecolor='w',
                    orientation='portrait', format=None,
                    transparent=False, bbox_inches='tight',  pad_inches=0.1,
                    metadata=None)
                plt.close()

            elif maga.ndim > 1:

                if maga.ndim == 3 and "statitics" not in a: 
                   time = ncfile.variables['time'][:]
                   hours=time/60
                
                if "Global" in a:
                    if v[0:2]=="pH":
                        titlea="\npH change - Full Area\n"
                    elif mode == 'Brine':
                        titlea="\nFull Area\n"
                    else:
                        titlea="\nDIC - Full Area\n"
                else:
                    title="Map from "+a[16:]+": "+v
                    titlea='\n'.join(wrap(title,60))
                    titlea=titlea+'\n'
                
                if maga.ndim > 2 and mode == 'CO2':
                   with open('temp.txt', 'a') as f:
                       outdata =' '.join(['Max Change in',v,a[16:],'-',str(maxchange)])
                       f.write(outdata)
                       f.write('\n')
                
                for val in Thresholds:
                
                    if maga.ndim > 2:
                
                        valarea = 0.0
                        largeareamax = 0.0
                    
                        if maga.max() > abs(val):
                    
                            if maga.ndim ==3:
                            
                                for idxx, dxpoints in enumerate(dx):
                                    for idxy, dypoints in enumerate(dy):
                                         if magmax[idxx,idxy] >= abs(val):
                                             valarea=valarea+(dx[idxx]*dy[idxy])
                            
                            with open('temp.txt', 'a') as f:                                
                                if mode == 'Brine':
                                   outdata =' '.join(['Impact-Area',a[16:],v,str(1/val),str(int(valarea))])
                                else:
                                   outdata =' '.join(['Impact-Area',a[16:],v,str(val),str(valarea)])   
                                f.write(outdata)
                                f.write('\n')
                
#                    if maga.ndim > 1:

                if maga.ndim == 2:
                
                    fig,ax = plt.subplots()
                    plt.set_cmap("jet")
                    cmap = cm.get_cmap('jet')
                    ax.set_facecolor(cmap(0.0))
               
#                        if abs(val) != maxchange:
#                           titleb=titlea+"(pH change of "+str(val)+")"
#                           plt.title(titleb)
#                        else:
                    plt.title(titlea)
                    plt.xlabel("X")
                    plt.ylabel("Y")
                
#                   levels = np.linspace(maga.max()-abs(val), maga.max(), 200)
                    if v[0:2]=="pH":
                       levels = np.linspace(-maxchange, 0, 200)
                    else:
                       levels = np.linspace(0, maxchange, 200)
                    if x.max() > 1000:
                        plt.ticklabel_format(axis="X", style="sci", scilimits=(0,0))
                    if y.max() > 1000:
                        plt.ticklabel_format(axis="Y", style="sci", scilimits=(0,0))
                    plt.grid(False)
                    plt.gca().set_aspect('equal', adjustable='box')
                    axpos=Bbox.from_extents(0.125, 0.1, 0.75, 0.9)
                    ax.set_position(axpos)
                    pos_x = axpos.x0+axpos.width + 0.05
                    pos_y = axpos.y0
                    cax_width = 0.04
                    cax_height = axpos.height
                    #create new axes where the colorbar should go.
                    #it should be next to the original axes and have the same height!
                    pos_cax = fig.add_axes([pos_x,pos_y,cax_width,cax_height])
                               
                    if v[0:2]=="pH":
                        plt.set_cmap("jet_r")
                        cmap = cm.get_cmap('jet_r')
                        ax.set_facecolor(cmap(1.0))
                    fmt = FormatStrFormatter("%.2f")
                    magaout1=np.swapaxes(mag,0,1)
                    if v[0:2]=="pH":
                        contourf=ax.contourf(x,y,-magaout1,levels=levels)
                    else: 
                        contourf=ax.contourf(x,y,magaout1,levels=levels)
                    plt.colorbar(contourf, cax=pos_cax, format=fmt)
                    plt.savefig(c+".png", dpi=my_dpi,
                        facecolor='w', edgecolor='w',
                        orientation='portrait', format=None,
                        transparent=False, bbox_inches='tight',
                        pad_inches=0.1, metadata=None)                                
                    plt.close()
                
                    fig,ax = plt.subplots()
                    plt.set_cmap("Greys")
                    #cmap = cm.get_cmap('Greys')
                        
                    plt.title(titlea)
                    plt.xlabel("X")
                    plt.ylabel("Y")
                    levels2 = np.array(Thresholds)
                    np.insert(levels2,0,maxchange)
#                    levels2=np.insert(levels2,0,maxchange)
                    if v[0:2]=="pH":
                       levels2=np.append(levels2,0.001)
                       ax.set_facecolor(cmap(0.001))
                    else:
                       levels2=np.append(levels2,0)
                       ax.set_facecolor(cmap(0))
                    levels2 = np.sort(levels2)
                      
                    leveltitle=1/levels2       
                    leveltitle[0]=0
#                    leveltitle[-1]=np.inf        
                                        
                    if x.max() > 1000:
                        plt.ticklabel_format(axis="X", style="sci", scilimits=(0,0))
                    if y.max() > 1000:
                        plt.ticklabel_format(axis="Y", style="sci", scilimits=(0,0))
                    plt.grid(False)
                    plt.gca().set_aspect('equal', adjustable='box')
                    axpos=Bbox.from_extents(0.125, 0.1, 0.75, 0.9)
                    ax.set_position(axpos)
                    pos_x = axpos.x0+axpos.width + 0.05
                    pos_y = axpos.y0
                    cax_width = 0.04
                    cax_height = axpos.height
                    #create new axes where the colorbar should go.
                    #it should be next to the original axes and have the same height!
                    pos_cax = fig.add_axes([pos_x,pos_y,cax_width,cax_height])
                    if v[0:2]=="pH":
                        plt.set_cmap("Greys")
                        cmap = cm.get_cmap('Greys')
                        ax.set_facecolor(cmap(0.001))
                    fmt = FormatStrFormatter("%.2f") 
                    magaout1=np.swapaxes(mag,0,1)
                        
                    if v[0:2]=="pH":
                        contourf=ax.contourf(x,y,abs(magaout1),locator=ticker.LogLocator(),\
                    levels=levels2)
                        plt.colorbar(contourf, cax=pos_cax, format=fmt)
                    else:
                        contourf=ax.contourf(x,y,abs(magaout1),levels=levels2,extend="max")
                        cbar=plt.colorbar(contourf, cax=pos_cax, format=fmt) 
                        cbar.ax.set_yticklabels(leveltitle) 
                        cbar.cmap.set_over('k')
                        
                    plt.savefig(c+"Thresholds"+".png", dpi=my_dpi,
                        facecolor='w', edgecolor='w',
                        orientation='portrait', format=None,
                        transparent=False, bbox_inches='tight',
                        pad_inches=0.1, metadata=None)
                    plt.close()
                
                if maga.ndim == 3:
                    fig,ax = plt.subplots()
                    plt.set_cmap("jet")
                    cmap = cm.get_cmap('jet')
                    ax.set_facecolor(cmap(0.0))
                    plt.title(titlea)
                    plt.xlabel("X")
                    plt.ylabel("Y")
                    if v[0:2]=="pH":
                        levels1 = np.linspace(-maxchange, 0, 200)
                    else:
                        levels1 = np.linspace(0, maxchange, 200)
                    if x.max() > 1000:
                        plt.ticklabel_format(axis="X", style="sci", scilimits=(0,0))
                    if y.max() > 1000:
                        plt.ticklabel_format(axis="Y", style="sci", scilimits=(0,0))
                    plt.grid(False)
                    plt.gca().set_aspect('equal', adjustable='box')
                    axpos=Bbox.from_extents(0.125, 0.1, 0.75, 0.9)
                    ax.set_position(axpos)
                    pos_x = axpos.x0+axpos.width + 0.05
                    pos_y = axpos.y0
                    cax_width = 0.04
                    cax_height = axpos.height
                    #create new axes where the colorbar should go.
                    #it should be next to the original axes and have the same height!
                    pos_cax = fig.add_axes([pos_x,pos_y,cax_width,cax_height])
                    if v[0:2]=="pH":
                        plt.set_cmap("jet_r")
                        cmap = cm.get_cmap('jet_r')
                        ax.set_facecolor(cmap(1.0))
                    fmt = FormatStrFormatter("%.3f") 
                    magaout1=np.swapaxes(magmax,0,1)
                    if v[0:2]=="pH":
                       magaout1=-magaout1
                    contourf=ax.contourf(x,y,magaout1,levels=levels1)
                    plt.colorbar(contourf, cax=pos_cax, format=fmt)
                    plt.savefig(c+".png", dpi=my_dpi,
                        facecolor='w', edgecolor='w',
                        orientation='portrait', format=None,
                        transparent=False, bbox_inches='tight',
                        pad_inches=0.1, metadata=None)
                    plt.close()
                        
                    fig,ax = plt.subplots()
                    plt.set_cmap("Greys")
                    #cmap = cm.get_cmap('Greys')
                        
                    plt.title(titlea)
                    plt.xlabel("X")
                    plt.ylabel("Y")
                    levels2 = np.array(Thresholds)
                    np.insert(levels2,0,maxchange)
#                    levels2=np.insert(levels2,0,maxchange)
                    if v[0:2]=="pH":
                       levels2=np.append(levels2,0.001)
                       ax.set_facecolor(cmap(0.001))
                    else:
                       levels2=np.append(levels2,0)
                       ax.set_facecolor(cmap(0))
                    levels2 = np.sort(levels2)
                      
                    leveltitle=1/levels2       
                    leveltitle[0]=0
#                    leveltitle[-1]=np.inf        
                                        
                    if x.max() > 1000:
                        plt.ticklabel_format(axis="X", style="sci", scilimits=(0,0))
                    if y.max() > 1000:
                        plt.ticklabel_format(axis="Y", style="sci", scilimits=(0,0))
                    plt.grid(False)
                    plt.gca().set_aspect('equal', adjustable='box')
                    axpos=Bbox.from_extents(0.125, 0.1, 0.75, 0.9)
                    ax.set_position(axpos)
                    pos_x = axpos.x0+axpos.width + 0.05
                    pos_y = axpos.y0
                    cax_width = 0.04
                    cax_height = axpos.height
                    #create new axes where the colorbar should go.
                    #it should be next to the original axes and have the same height!
                    pos_cax = fig.add_axes([pos_x,pos_y,cax_width,cax_height])
                    if v[0:2]=="pH":
                        plt.set_cmap("Greys")
                        cmap = cm.get_cmap('Greys')
                        ax.set_facecolor(cmap(0.001))
                    fmt = FormatStrFormatter("%.2f") 
                    magaout1=np.swapaxes(magmax,0,1)
                        
                    if v[0:2]=="pH":
                        contourf=ax.contourf(x,y,abs(magaout1),locator=ticker.LogLocator(),\
                    levels=levels2)
                        plt.colorbar(contourf, cax=pos_cax, format=fmt)
                    else:
                        contourf=ax.contourf(x,y,abs(magaout1),levels=levels2,extend="max")
                        cbar=plt.colorbar(contourf, cax=pos_cax, format=fmt) 
                        cbar.ax.set_yticklabels(leveltitle) 
                        cbar.cmap.set_over('k')
                        
                    plt.savefig(c+"Thresholds"+".png", dpi=my_dpi,
                        facecolor='w', edgecolor='w',
                        orientation='portrait', format=None,
                        transparent=False, bbox_inches='tight',
                        pad_inches=0.1, metadata=None)
                    plt.close()
                    
                if maga.ndim == 3 and 'No' not in Animation and 'fields' in a:
                    if "Global" not in Animation or 'global' in a:
                        print("Producing animation for variable: {}".format(v))
                        folder = c#.replace("/external/output/", "/")
                        folder = folder.replace(" ", "")
                        os.mkdir(folder)
                       
                        for i, hour in enumerate(hours):
                            fig,ax = plt.subplots()
                            plt.set_cmap("jet")
                            cmap = cm.get_cmap('jet')
                            ax.set_facecolor(cmap(0.0))
                            plt.title(titlea)
                            plt.xlabel("X")
                            plt.ylabel("Y")
                            if v[0:2]=="pH":
                               levels = np.linspace((-maxchange)*0.7, (-maxchange)*0.0001, 200)
                            else:
                               levels = np.linspace(0, maxchange, 200)
                            if x.max() > 1000:
                                plt.ticklabel_format(axis="X", style="sci", scilimits=(0,0))
                            if y.max() > 1000:
                                plt.ticklabel_format(axis="Y", style="sci", scilimits=(0,0))
                            plt.grid(False)
                            plt.gca().set_aspect('equal', adjustable='box') 
                            axpos=Bbox.from_extents(0.125, 0.1, 0.75, 0.9)
                            ax.set_position(axpos)
                            pos_x = axpos.x0+axpos.width + 0.05
                            pos_y = axpos.y0
                            cax_width = 0.04
                            cax_height = axpos.height
                            #create new axes where the colorbar should go.
                            #it should be next to the original axes and have the same height!
                            pos_cax = fig.add_axes([pos_x,pos_y,cax_width,cax_height])
                            if v[0:2]=="pH":
                                plt.set_cmap("jet_r")
                                cmap = cm.get_cmap('jet_r')
                                ax.set_facecolor(cmap(1.0))
                            fmt = FormatStrFormatter("%.3f")
                            magaout=np.swapaxes(mag,1,2)
                            if v[0:2]=="pH":
                               magaout=-magaout
                               
                            if "Global" in a:
                                fig.set_size_inches(10.24, 7.68, True)
                            filenm = str(i)+'.png'
                            if v[0:2]=="pH":
                                contourf=ax.contourf(x,y,magaout[i,:,:],levels=levels,extend="both",vmin=(-maxchange)*0.7, vmax=(-maxchange)*0.0001)
                            else:
                                contourf=ax.contourf(x,y,magaout[i,:,:],levels=levels)
                            timeframe="Time: {:.2f} hours".format(hours[i])
                            plt.title(timeframe)
                            cmap.set_over('white')
                            plt.colorbar(contourf, cax=pos_cax, format=fmt)
                            plt.savefig(folder+'/'+filenm, dpi=my_dpi)
                            plt.close("all")
                            
                        images = []
                            
                        files2=os.listdir(folder)
                        files2 = [files3[:-4] for files3 in files2]
                        files3= [int(i) for i in files2]
                        files3=sorted(files3)
                            
                        for filename in files3:
                            filename2=folder+'/'+str(filename)+".png"
                            images.append(imageio.v2.imread(filename2))
                        imageio.mimsave(c+'.gif',images,duration=0.1) # duration sets the speed of the gif
    
    print("\n")
print("------------------------------------------------------")
lines_seen = set() # holds lines already seen
outfile = open('/external/output/Impact.txt', "w")
for line in open('temp.txt', "r"):
    if line not in lines_seen: # not a duplicate
        outfile.write(line)
        lines_seen.add(line)
outfile.close()
print("\nTada..\n")
