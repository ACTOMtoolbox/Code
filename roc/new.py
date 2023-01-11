import netCDF4 as nc
import os
     
diff=10**20      
err=0    
                
print("\n******************************************************")
print("\nProcessing Rate of Change Anomally\n")
print("******************************************************")

# Open simulation output files in a loop, input = 1.0 (units/s)/kg water

dataPath1='input/external/'
if len(os.listdir(dataPath1)) == 0:
   dataPath1='input/netcdf/'

if os.path.isdir(dataPath1):
    filelist = os.listdir(dataPath1)
    filelist = sorted(filelist)
    print(" ")
    print("------------------------------------------------------\n")
    for x in filelist:
        filetype = x[-3:]
        if filetype == ".nc":
        
           print( "Opening File: {}".format(x))
           fn = dataPath1 + x
           ds = nc.Dataset(fn,'a')
           TP=ds["time"][:]    
           TP=TP.clip(min=0)
           PV=ds["DIC"][:]    
           PV=PV.clip(min=0)
           PV1=ds["DIC-kgm3"][:]    
           PV1=PV1.clip(min=0)
           PV2=ds["pH"][:]    
           PV2=PV2.clip(min=0)
           
           TP=(TP-TP[0])
           BinCount=[0]
           MaxInBin=[0]
           MaxInBinD=[0]
           MaxInBinpH=[0]
           Total=[0]
           
           n = len(TP)
           
           diff = abs(TP[1] - TP[0])

           for i in range(n-1):
               for j in range(i+1,n):
                   if abs(TP[i]-TP[j]) < diff:
                       diff = abs(TP[i] - TP[j])
           
           Binwidth = int(diff+1)*2
           
           for I in range(n-1):
               for J in range(I+1,n):
                   DeltaTP = TP[J] - TP[I]
                   if DeltaTP<(86400/2):     # 86400 = number of seconds in a day, prevents comparison of values more than a day apart
                       DeltaPV = max(abs((PV[J] - PV[I])),0)
                       DeltaPV1 = max(abs((PV1[J] - PV1[I])),0)
                       DeltaPV2 = max(abs((PV2[J] - PV2[I])),0)
# Chuck into bins
                       Bin = int(DeltaTP/Binwidth)	
                       if len(BinCount)-1 < Bin:
                           while len(BinCount)-1 < Bin:
                               BinCount.append(0)
                               Total.append(0)
                               MaxInBin.append(0)
                               MaxInBinD.append(0)
                               MaxInBinpH.append(0)
                              
                       BinCount[Bin] = BinCount[Bin] + 1
                       Total[Bin] = Total[Bin] + DeltaPV
                       MaxInBin[Bin] = max(MaxInBin[Bin], DeltaPV)
                       MaxInBinD[Bin] = max(MaxInBinD[Bin], DeltaPV1)
                       MaxInBinpH[Bin] = max(MaxInBinpH[Bin], DeltaPV2)

           fileout="/external/output/ROCAC-"+x[:-3]

           with open(fileout, 'w') as f:
               for I in range(len(BinCount)):
                   f.write('{} {} {} {} \n'.format(I+1,(I+1)*Binwidth,BinCount[I],MaxInBin[I],))
                   
           with open(fileout+'-kgm3', 'w') as f:
               for I in range(len(BinCount)):
                   f.write('{} {} {} {} \n'.format(I+1,(I+1)*Binwidth,BinCount[I],MaxInBinD[I],))
           
           with open(fileout+"-pH", 'w') as f:
               for I in range(len(BinCount)):
                   f.write('{} {} {} {} \n'.format(I+1,(I+1)*Binwidth,BinCount[I],MaxInBinpH[I],))


