import netCDF4 as nc
import sys
import csv
import os
import configparser
import PyCO2SYS as pyco2
import numpy as np
import seawater as sw

config=configparser.ConfigParser(allow_no_value=True)
config.read('/external/settings/data.ini')
TA = float(config['General']['TA']) # Sensor measurement frequency
J=1

fn = 'input/test.nc'
ds = nc.Dataset(fn, 'w', format='NETCDF4')
time = ds.createDimension('time', None)

if os.path.isdir('input/'):

    filename = open('input/STEMM-CCS.csv','r')
    Pout=[]
    Tout=[]
    Sout=[]
    pHout=[]

        
    file = csv.DictReader(filename)
    for col in file:
        Pout.append(col['pres'])
        Tout.append(col['temp'])
        Sout.append(col['sal'])
#            pH1in = list(zip(*reader))[10]
#            pH2in = list(zip(*reader))[11]
        pHout.append(col['pH'])
    
    timein=list(range(1,len(pHout)))
    timeout=l = [x * 20 for x in timein]
    
    
    Pout = [float(x) for x in Pout]
    Tout = [float(x) for x in Tout]
    Sout = [float(x) for x in Sout]
    pHout = [float(x) for x in pHout]
    
#    for I in range(len(Pout)):
#                if pH1in[I] not = 0 and pH2in[I] not = 0:
#                if pHin[I] is not 0:
#                    Tout[J]=Tin[I]
#                    Sout[J]=Sin[I]
#                    Pout[J]=Pin[I]
#                    pH1out[J]=pH1in[I]
#                    pH2out[J]=pH2in[I]
#                    pHout[J]=pHin[I]
#        timeout[I]=(I+1)/86400
#                    J=J+1
#
#            pHout = (pH1out+pH2out)/2
        
    kwargs = dict(
        par1 = pHout,
        par2 = TA,
        par1_type = 3,
        par2_type = 1,
        temperature = Tout,
        salinity = Sout,
        pressure = Pout,
        )
        
    output=pyco2.sys(**kwargs)
    DICout=output["dic"]
    
    Dens=sw.dens(Sout,Tout,Pout)
    
    DICout2=DICout*Dens
    DICout2=[x * 0.0440095/1000000 for x in DICout2]
        
    time = ds.createVariable('time', 'f4', ('time',))
    pH = ds.createVariable('pH', 'f4', ('time',))
    DIC = ds.createVariable('DIC','f4', ('time',))
    DIC.units = 'micromol/kg'
    DIC.standard_name = 'Dissolved Inorganic Carbon'

    DIC2 = ds.createVariable('DIC-kgm3','f4', ('time',))
    DIC2.units = 'kg/m3'
    DIC2.standard_name = 'Dissolved Inorganic Carbon'

    DIC2[:] = DICout2      
    DIC[:] = DICout
    pH[:] = pHout
    time[:] = timeout

                
ds.close()


                          
                  


            
            
            
          
