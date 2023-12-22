#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Optimal cover

Each set is translated to another center of the claster.

@author: aol053
"""

import numpy as np
import matplotlib.pyplot as plt

import PyCO2SYS as pyco2
from mip import Model, xsum, minimize, BINARY, OptimizationStatus

from netCDF4 import Dataset
import xarray as xr

import configparser
import os

from distutils.util import strtobool

#%% Reading ini file

parent_dir = os.getcwd()
parent_dir






th=[]
thtit=[]

config=configparser.ConfigParser(allow_no_value=True)
config.read('/external/settings/data.ini')

mode = config['General']['run'] # brine? 

q = float(config['General']['Rate']) # flux 
qunits = config['General']['Rate-Units'] # flux units

if mode == 'CO2':
    if config.has_option('General', 'ta'):
        TAmean = float(config['General']['ta'])
    else:
        TAmean = float(config['CSEEP']['ta-mean'])
    
    if config.has_option('General', 'dic'):
        DICmean = float(config['General']['dic'])
    else:
        DICmean = float(config['CSEEP']['dic-mean'])
    
    if config.has_option('General', 'dens'):
        dens = float(config['General']['dens'])
    else:
        dens = float(config['CSEEP']['dens-mean'])
    
    if config.has_option('General', 'threshold-pH'):
        th1=((config['General']['threshold-ph']).split(","))
        for i in range(len(th1)):
          kwargs = dict(
            par1 = DICmean,
            par2 = TAmean,
            par1_type = 2,
            par2_type = 1,
            )
          outpH=pyco2.sys(**kwargs)
          phthres=outpH["pH_total"]
          phthres=phthres+float(th1[i])
          kwargs = dict(
            par1 = phthres,
            par2 = TAmean,
            par1_type = 3,
            par2_type = 1,
            )
          outDIC=pyco2.sys(**kwargs)
          DICthres=outDIC["dic"] 
          DICthres=abs(DICthres-DICmean)
          DICthres=DICthres*dens*0.0440095/1000000 
          thtit=thtit+['USER']
          th.append(DICthres)
    
    if config.has_option('RateOfChange', 'threshold'):
        th1=((config['RateOfChange']['threshold']).split(","))
        for i in range(len(th1)):
            thtit=thtit+['ROC']
            th.append(float(th1[i]))

    if config.has_option('CSEEP', 'threshold'):
        th1=((config['CSEEP']['threshold']).split(","))
        for i in range(len(th1)):
            thtit=thtit+['CSEEP']
            th.append(float(th1[i]))
            
#Rate Conversion

    RateUnitsTimeIn=qunits.partition('/')
    TimeIn=RateUnitsTimeIn[2]
    MassIn=RateUnitsTimeIn[0]
       
    if MassIn[0] == "k":
       q = q
    if MassIn[0] == "g":
       q = q/1000
    if MassIn[0] == "t":
       q = q*1000
    if MassIn[0] == "m":
       q = q*0.04401
    
    if TimeIn[0] == "m":
        if TimeIn[1] == "i":
            q = q/60
        if TimeIn[1] == "o":
            q = q/2629800
    if TimeIn[0] == "h":
        q = q/3600
    if TimeIn[0] == "d":
        q = q/86400
    if TimeIn[0] == "w":
        q = q/604800
    if TimeIn[0] == "y":
        q = q/31557600            
            
elif mode == 'Brine':
    th1=((config['General']['threshold-pH']).split(","))
    for i in range(len(th1)):
        thtit=thtit+['Brine']
        th.append(1.0/float(th1[i]))
        
        

include_cluster_points = bool(strtobool(config['OptCover']['include_cluster_points']))
Timelim = int(config['OptCover']['Time_limit']) # optimization Time in seconds
cost_function = int(config['OptCover']['cost_function'])

#%%

def cost_prob(A,p):
    p=p/sum(p); # normized to be sure
    P=np.diag(p);
    z=1-sum(np.dot(P,A));
    return z    

#%%

def find_ind(x0,x):
    k0=np.argmin(np.abs(x-x0));
    return k0

def formingS(S0, x0, x1, x, y0, y1, y):
    #S0 is an (Nx x Ny)  of 0s and 1s
    Nx, Ny = S0.shape;
    S1 = np.zeros((Nx,Ny));
    i0 = find_ind(x1,x)-find_ind(x0,x); 
    j0 = find_ind(y1,y)-find_ind(y0,y);
    # an ugly way, that I do not like!
    for i in range(Nx):
        for j in range(Ny):
            if ((i-i0)>-1)& ((i-i0)<Nx)&((j-j0)>-1)& ((j-j0)<Ny):
                S1[i,j]=S0[i-i0,j-j0]
    
    return S1      
     
#test (works well!)

# x=np.array([1,2,3,4]);
# y=np.array([1,2,3,4]);
# x0=3; y0=2;
# S0=np.zeros((4,4)); S0[2,1]=1; S0[2,2]=1; S0[1,1]=1;


# x1=1; y1=4;
# S1=formingS(S0,x0,x1,x,y0,y1,y)

# fig, ax = plt.subplots(nrows=1,ncols=2)
# ax[0].pcolor(np.transpose(S0))
# ax[1].pcolor(np.transpose(S1))


#%% Input Nl x Nx x Ny nonegative values


#d = Dataset('Input/from_A_to_B.nc');
d = Dataset('Input/statistics_global.nc');

S=np.array(d.variables['max'][:]) # in (max,y,x)

# in km
xx=np.array(d.variables['x'][:])
yy=np.array(d.variables['y'][:])


x_source=np.array(d.variables['x_source'][:]);
y_source=np.array(d.variables['y_source'][:]);

#
pvec=d.variables['location_probability'][:];

## 
X,Y=np.meshgrid(xx,yy);
Xflat=X.flatten();
Yflat=Y.flatten();

#S=np.swapaxes(Sin,1,2)

#print(S.shape)
#print(xx.shape)
#print(yy.shape)
#print(x_source.shape)
#print(y_source.shape)
#print(pvec.shape)
#print(X.shape)

for thres in range(len(th)):

    #%%
    Sets=1.*(S>float(th[thres])/q);
    
    (Nl,Nx,Ny)=Sets.shape;
    
#    th2 = [elem[:6] for elem in th]

    if mode == 'Brine':
       titl=' '.join([thtit[thres],'- dilution factor of',str(1/th[thres])])
    if mode == 'CO2':    
       titl=' '.join([thtit[thres],'-',str(round(th[thres], 6)),'kg/m^2'])
    
    print(titl)
    
    #%% Translate sets for each point in the cluster
    
    if include_cluster_points : # True = include
        cluxr = xr.open_dataset('Input/point_tags_reduced.nc');
        gr = cluxr.groupby(cluxr.source);
        BigSets=np.empty((0,Ny,Nx),float);
        prob=np.empty((0,),float);
        for i in gr.groups.keys():
            print(i)
            S0=Sets[i,:,:];
            BigSets=np.append(BigSets,[S0],axis=0);
            prob = np.append(prob, pvec[i]);
            x_in_ith=np.array(cluxr.x[gr.groups[i]])
            y_in_ith=np.array(cluxr.y[gr.groups[i]])
            
            n_in_ith=len(x_in_ith);
            for k in range(len(x_in_ith)):
                Sk=formingS(S0,x_source[i],x_in_ith[k],xx,y_source[i],y_in_ith[k],yy)
                BigSets=np.append(BigSets,[Sk],axis=0);
                prob = np.append(prob, pvec[i]);
        Nl=BigSets.shape[0];
        Sets=BigSets;
        pvec=prob;
    
    #%%
    
    
    A = np.reshape(Sets,[Nl,Nx*Ny])
    
    C, iA, iC= np.unique(A,return_index=True,return_inverse=True, axis=1); #remove repeated columns
    
    m=C.shape[1];
    
    
    f_1=np.ones((m,1)); #minimal number
    
    f_over=(len(C)-np.sum(C,axis=0)).reshape(m,1); #overdetection
    
    pvec=pvec/sum(pvec);
    f_p=cost_prob(C,pvec).reshape(m,1); #max probability
    
    
    #%% Choose the cost function: f1, fp, f_over
    
    f=f_over
    #%% Optimization model
    I = range(m);
    
    mo = Model();
    
    x = [mo.add_var(var_type=BINARY) for i in I]
    
    Cl=C.tolist();
    fl=f.tolist();
    
    mo.objective = minimize(xsum(fl[i][0]*x[i] for i in I))
    
    for k in range(Nl):
        mo += xsum(Cl[k][i] * x[i] for i in I) >= 1;
    
    
    #mo.emphasis=1
    mo.verbose = 0
    status = mo.optimize(max_seconds=Timelim);
    if status == OptimizationStatus.OPTIMAL:
        print('optimal solution cost {} found'.format(mo.objective_value))
    elif status == OptimizationStatus.FEASIBLE:
        print('sol.cost {} found, best possible: {}'.format(mo.objective_value, m.objective_bound))
    elif status == OptimizationStatus.NO_SOLUTION_FOUND:
        print('no feasible solution found, lower bound is: {}'.format(mo.objective_bound))
    elif status == OptimizationStatus.INFEASIBLE:
        print('no feasible solution found, lower bound is: {}'.format(mo.objective_bound))
    if status == OptimizationStatus.OPTIMAL or status == OptimizationStatus.FEASIBLE:
        print('solution:')
        for v in mo.vars:
           if abs(v.x) > 1e-6: # only printing non-zeros
              print('{} : {}'.format(v.name, v.x))
    
        selected = [i for i in I if x[i].x >= 0.99]
        print("selected items: {}".format(selected))
    
        values = [x[i].x for i in I]
    
    #%%
    
        xsol=np.array(values).reshape(m,1); # solution vector
        Ns=int(sum(xsol[:])[0]); # number of sensors
        C2=C[:,selected];
        leaks_per_sensor = np.sum(C2,axis=0).astype(int)
    
        iS=iA[selected]; # in the original matrix A
    
        #p=np.matmul(C,xsol)
    
        #Nover=Ns*Nl-sum(f_over*xsol)-Nl #total overdetection
        #print('Total over-detection: ', Nover[0])
    
    
        #%% visualization
    
        list_of_sources=[];
        list_of_pins=[];
        D=[]; Sums=np.zeros([Ny,Nx]);
        
        for i in range(len(leaks_per_sensor)):
            indS=np.nonzero(C2[:,i]==1)[0];
            list_of_sources.append(indS)
            aux=np.prod(A[indS,:],axis=0).reshape(Nx,Ny);
            aux1=np.swapaxes(aux,0,1)
            
            D.append(aux1);
            Sums=Sums+aux1;
        
            xx=Xflat[aux1.flatten()==1];
            yy=Yflat[aux1.flatten()==1];
            xpin=xx.mean();
            ypin=yy[np.argmin(abs(xx-xpin))];
            list_of_pins.append([xpin,ypin]);
        
        
        given=dict(facecolor='white', boxstyle='circle')
        
        
        
        #%% FIGURE 1
        eps=0/10;
        sc=1000; # in km
        
        fig,ax=plt.subplots()
        aux=np.reshape(np.sum(A,axis=0),(Nx,Ny));
        
        #img=ax.imshow(aux,origin='lower')
        aux1=np.swapaxes(aux,0,1)
        img=ax.pcolor(X/sc,Y/sc,aux1,shading='auto')
        fig.colorbar(img)
        iS=iA[selected]; # indices of x ans y, in the original matrix A
        
        if include_cluster_points:
            ax.plot(cluxr.x/sc,cluxr.y/sc,'b.')
        ax.plot(x_source/sc,y_source/sc,'r.')
        
        #ax.plot(Xflat[iS]/sc,Yflat[iS]/sc,'r*')
        
        titl=' '.join(['Stationary sensor placement, N=',str(float(sum(xsol)))])
        ax.set_title(titl)
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        for i in range(len(leaks_per_sensor)):
            #ax.text(Xflat[iS[i]]/sc+eps,Yflat[iS[i]]/sc+eps, int(leaks_per_sensor[i]),backgroundcolor='white', size=6,
            #        bbox=given)
            ax.text(list_of_pins[i][0]/sc+eps,list_of_pins[i][1]/sc+eps, int(leaks_per_sensor[i]),backgroundcolor='white', size=6,
                    bbox=given)


        
        
        if mode == 'Brine':
           titl=' '.join([thtit[thres],'-',str(1/th[thres])])
        if mode == 'CO2':    
           titl=''.join([thtit[thres],'-',str(round(th[thres], 6))])
        output=''.join(['Output/',titl,'-map.png'])
        output.replace(" ", "")     
        plt.savefig(output)
        
        #%% FIGURE 2
        eps=2/10;
        sc=1000; # in km
        
        fig,ax=plt.subplots()
        
        ax.pcolor(X/sc,Y/sc,Sums,shading='auto')
        for i in range(Ns):
            indx=D[i].flatten()==1;
            ax.plot(Xflat[indx]/sc,Yflat[indx]/sc,'.')
            # ax.text(Xflat[iS[i]]/sc+eps,Yflat[iS[i]]/sc+eps, int(leaks_per_sensor[i]),backgroundcolor='white', size=6,
            #         bbox=given)
            ax.text(list_of_pins[i][0]/sc+eps,list_of_pins[i][1]/sc+eps, int(leaks_per_sensor[i]),backgroundcolor='white', size=6,
                    bbox=given)
            
        
        # if include_cluster_points=='yes':
        #     ax.plot(cluxr.x/sc,cluxr.y/sc,'b.')
        
        # ax.plot(x_source/sc,y_source/sc,'r.')
        
        titl=' '.join(['Placement areas for each of the',str(Ns),'sensors'])
        ax.set_title(titl)
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        
        if mode == 'Brine':
           titl=' '.join([thtit[thres],'-',str(1/th[thres])])
        if mode == 'CO2':    
           titl=''.join([thtit[thres],'-',str(round(th[thres], 6))])
        output=''.join(['Output/',titl,'-map2.png'])
        output.replace(" ", "")     
        plt.savefig(output)
        
        #%% Creating NC file, what to include??
        #list_of_sources=np.array(list_of_sources,dtype='object');
        #list_of_sources.astype('object')
        
        x_all=x_source;
        y_all=y_source;
        
        if include_cluster_points:
            x_all=np.append(x_all,cluxr.x);
            y_all=np.append(y_all,cluxr.y);
                
        data_vars = {'x_sensor' : (['sensor'], Xflat[iS]),
                     'y_sensor' : (['sensor'], Yflat[iS]),
                      #'detected_per_sensor': (['sensor'], list_of_sources),
                      'ns_per_sensor': (['sensor'], leaks_per_sensor),
                      'x_source':(['source'], x_all),
                      'y_source':(['source'], y_all),
                      'N_over':([], int(sum(leaks_per_sensor)-Nl)),}
        
        coords = {'sensor'    : (['sensor'], list(range(1,Ns+1)) ),
                  'source': (['source'],list(range(1,Nl+1)))}
        
        dataset = xr.Dataset(data_vars=data_vars,coords=coords)
        
        output=''.join(['Output/',titl,'-sensor_locations.nc'])
        output.replace(" ", "")
        
        dataset.to_netcdf(output)# allow_pickle=True)
