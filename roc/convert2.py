import netCDF4 as nc
import sys
import os
import configparser
import PyCO2SYS as pyco2
import numpy as np

old=1
TA=0

config=configparser.ConfigParser(allow_no_value=True)
config.read('/external/settings/data.ini')
TA = float(config['General']['TA']) # Sensor measurement frequency
J=1


if os.path.isdir('input/'):
    filelist = os.listdir('input/')
    filelist = sorted(filelist)

    for x in filelist:
        filetype = x[-3:]
        if filetype == ".nc":
       
            fn = 'input/' + x
            ds = nc.Dataset(fn,'a')
            variable=ds.variables.keys()
            variables=list(variable)
           
            if "pH" in variables and "DIC" not in variables:

                if old==1:
                    print("\n******************************************************")
                    print("\nProcessing Carbonate System\n")
                    print("******************************************************\n")
                    old=old+1
        
                print( "Opening File: {}".format(x))
                
                Carba=ds.variables["pH"]
                Carb=ds["pH"][:] 

                print("Solving the carbonate system for pH")
                print("  This may take multiple minutes depending on the file size")
                print(Carb)
                print(TA)
                kwargs = dict(
                    par1 = Carb,
                    par2 = TA,
                    par1_type = 3,
                    par2_type = 1,
                    )

                output=pyco2.sys(**kwargs)
                DICout=output["dic"]
                
                DIC = ds.createVariable('DIC',np.float64,Carba.dimensions)
                DIC.units = 'micromol/kg'
                DIC.standard_name = 'Dissolved Inorganic Carbon'
                
                if DICout.ndim == 1:
                    DIC[:] = DICout
                if DICout.ndim == 2:
                    DIC[:,:] = DICout
                if DICout.ndim == 3:
                    DIC[:,:,:] = DICout
                if DICout.ndim == 4:
                    DIC[:,:,:,:] = DICout
                print(" ")
                print("Tada..")


                          
                  


            
            
            
          
