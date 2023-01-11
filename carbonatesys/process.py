import netCDF4 as nc
import sys
import os
import PyCO2SYS as pyco2
import numpy as np
import configparser

print("\n******************************************************")
print("\nProcessing Carbonate System\n")
print("******************************************************\n")

config=configparser.ConfigParser(allow_no_value=True)
config.read('/external/settings/data.ini')

if config.has_option('General', 'ta'):
    TA = float(config['General']['ta'])
else:
    TA = float(config['CSEEP']['ta-mean'])
    
if config.has_option('General', 'dic'):
    DICback = float(config['General']['dic'])
else:
    DICback = float(config['CSEEP']['dic-mean'])

Rate = float(config['General']['rate'])
RateUnits = config['General']['rate-units']
Dens = float(config['CSEEP']['dens-mean'])
CseepThreshold = float(config['CSEEP']['threshold_dic_micromol-kg'])
localglobal = config['Impacts']['local']
            
#            print(" ")
#            print("------------------------------------------------------\n")
#            print("Solving the carbonate system for Cseep threshold")
#            CseepThreshold=CseepThreshold+DICback
#            output=pyco2.sys(par1=CseepThreshold, par2=TA, par1_type=2,par2_type=1)
#            output2=pyco2.sys(par1=DICback, par2=TA, par1_type=2,par2_type=1)
#            CseepThresholdpH=(output2["pH_total"]-output["pH_total"])
#            print("giving {} units ".format(CseepThresholdpH))
#            print( " ")
#            del output
#            del output2
# fix later !!!!!!!!!!!!!!!!!           


print( "Background TA = {}".format(TA))
print( "Background DIC = {}".format(DICback))
print( "Seawater Density = {}".format(Dens))
print( "Leakage Rate = {} {}".format(Rate, RateUnits))

# Multiply output by leakage rate ! if not using varied rate
    
# convert units to mol/Kg first

Factor1 = 0
Factor = 0

RateUnitsTime=RateUnits.partition('/')
Time=RateUnitsTime[2]
Mass=RateUnitsTime[0]

if Mass[0] == "k":
   Factor1 = Rate/(0.04401*Dens)
elif Mass[0] == "g":
   Factor1 = Rate/(1000*0.04401*Dens)
elif Mass[0] == "t":
   Factor1 = Rate*1000/(0.04401*Dens)
elif Mass[0] == "m":
   Factor1 = Rate/Dens

# convert units to micromol/Kg

Factor1=Factor1*1000000

if Factor1 == 0:
    print("Invalid leakage rate \"{}\" or mass units \"{}\"".format\
    (Rate, Mass))
    exit() 

if Time[0] == "s":
   Factor = Factor1
elif Time[0] == "m":
    if Time[1] == "i":
        Factor = Factor1/60
    if Time[1] == "o":
        Factor = Factor1/2628
elif Time[0] == "h":
    Factor = Factor1/3600
elif Time[0] == "d":
    Factor = Factor1/86400
elif Time[0] == "w":
    Factor = Factor1/604800
elif Time[0] == "y":
    Factor = Factor1/31557600

if Factor == 0:
    print("Invalid leakage rate \"{}\" or time units \"{}\"".format(Rate, Time))
    exit() 

# Open simulation output files in a loop, input = 1.0 (units/s)/kg water
if os.path.isdir('/external/input'):
    filelist = os.listdir('/external/input')
    filelist = sorted(filelist)
    print( " ")
    print("------------------------------------------------------\n")
    for x in filelist:

        if "mass" in x or "point_tags" in x:
#        if "mass" in x or "point_tags" in x or "fields_global" in x:
            continue
            
        if "Global" in localglobal and "global" not in x:
            if "cumulative" not in x:
               continue
        
        os.system("cp /external/input/\'{}\' /external/output/\'{}\'".format(x, x))  
        filetype = x[-3:]
        if filetype == ".nc":
           print( "Opening File: {}".format(x))
# Open file
           Carb=0
           fn = '/external/output/' + x
           ds = nc.Dataset(fn,'a')
           variable=ds.variables.keys()
           variables=list(variable) 
           print("Variables in file: {}".format(variables))
           
           for y in variables:
           
               if y == "u" or y == "v" or y == "x" or y == "y" or \
               y == "lon" or y == "lat" or y == "time" or y == "indx" or \
               y == "x_grid" or y == "y_grid" or y == "x_pos" or \
               y == "location_probability" or y == "x_source" or y == "y_source" or \
               y == "lon_source" or y == "lat_source" or y == "lat_source" or \
               y == "y_pos" or y == "dist" or y == "dC" or "delta" in y or  \
               y == "location probabilities" or y == "source tag" or y == "dt" or \
               y == "probe":
                   continue
               else:
                   Carba=ds.variables[y]
                   Carb=ds[y][:]    
                   Carb=Carb.clip(min=0)
                   
# caculate pH (total scale)
                   print("Solving the carbonate system for \'{}\'".format(y))
                   print("  This may take multiple minutes depending on the file size")

                   if "global" in x:
                       Carb2 = np.sum(Carb, axis=0)
                       Carba1 = Carba.dimensions[1:]
# Add concentration to DIC
                       DICout=(Carb2*Factor)+DICback
                   else:
# Add concentration to DIC
                       DICout=(Carb*Factor)+DICback
                       Carba1 = Carba.dimensions
                       
                   pHout=np.zeros(DICout.shape)
                   

                   if DICout.ndim == 4:
                       for countdown in range(DICout.shape[0]):                
                           output=pyco2.sys(par1=DICout[countdown,:,:,:],\
                           par2=TA, par1_type=2, par2_type=1)
                           pHout[countdown,:,:,:]=output["pH_total"]
                           del output
                   
                   elif DICout.ndim == 3:
                       for countdown in range(DICout.shape[0]):
                           output=pyco2.sys(par1=DICout[countdown,:,:],\
                           par2=TA, par1_type=2, par2_type=1)
                           pHout[countdown,:,:]=output["pH_total"]
                           del output
                           
                   else:
                           output=pyco2.sys(par1=DICout, par2=TA, par1_type=2,par2_type=1)
                           pHout=output["pH_total"]
                           del output
# print output
                   DICName='DIC from '+y
                   pHName='pH from '+y
                   
                   print( "Saving pH and DIC data to {}".format(x))
                   
                   DIC = ds.createVariable(DICName,np.float64,Carba1)
                   DIC.units = 'micromol/kg'
                   DIC.standard_name = 'Dissolved Inorganic Carbon'
                   pH = ds.createVariable(pHName,np.float64,Carba1)
                   pH.standard_name = 'Total Scale pH'
                   if DICout.ndim == 1:
                       DIC[:] = DICout
                       pH[:] = pHout
                   if DICout.ndim == 2:
                       DIC[:,:] = DICout
                       pH[:,:] = pHout
                   if DICout.ndim == 3:
                       DIC[:,:,:] = DICout
                       pH[:,:,:] = pHout
                   if DICout.ndim == 4:
                       DIC[:,:,:,:] = DICout
                       pH[:,:,:,:] = pHout
                   print(" ")
        else:
           print( "Ignoring non-NetCDF File/Dir: {}".format(x))
           print( " ")
           
        print("------------------------------------------------------\n")
else:
# use test value
    print( "using test values")
    Carb=10
# Add concentration to DIC
    DICout=(Carb*Factor)+DIC
# caculate pH (total scale)
    output=pyco2.sys(par1=DICout, par2=TA, par1_type=2, par2_type=1)
    pHout=output["pH_total"]
    del output
# print output
    print( "Saving ph data to {}".format(pHout))
    
print("Tada..")

