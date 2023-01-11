#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 21:20:03 2022

@author: jean
"""

import glob
import os
import numpy as np
import pandas as pd
import netCDF4 as nc
import pickle
import gsw

def readFile(dataFile):
    import time
    fid = nc.Dataset(dataFile)
    lat = fid.variables['lat']
    lon = fid.variables['lon']
    
    [tmpyy,tmpxx] = np.meshgrid(lon,lat)
    
    tmptime = fid.variables['time']
    tmpdyp = fid.variables['depth']
    
    tmppH = fid.variables['ph']
    tmpphosphate = fid.variables['po4']
    tmpsalinity = fid.variables['so']
    tmpoxygen = fid.variables['o2']
    tmptemperature = fid.variables['thetao']
    tmpnitrate = fid.variables['no3']
    tmpchla = fid.variables['chl']
    tmpta = fid.variables['O3_TA']
    tmpdic = fid.variables['O3_c']
    tmpsilicate = fid.variables['si']
    tmphor_spd = fid.variables['uo']
    tmpvrt_spd = fid.variables['vo']
    tmpnpp = fid.variables['nppv']
    tmpphyc = fid.variables['phyc']
    tmpppco2 = fid.variables['pco2']
    
    pH = []
    phosphate = []
    salinity = []
    oxygen = []
    temperature = []
    nitrate = []
    chla = []
    ta = []
    dic = []
    silicate = []
    hor_spd = []
    vrt_spd = []
    npp = []
    phyc = []
    ppco2 = []
    
    xx = []
    yy = []
    tid = []
    year = []
    month = []
    day = []
    dyp = []
    
    for d in range(len(tmpdyp)):
        for t in range(len(tmptime)):
            pH.append(tmppH[t][d].reshape(-1))
            phosphate.append(tmpphosphate[t][d].reshape(-1))
            salinity.append(tmpsalinity[t][d].reshape(-1))
            oxygen.append(tmpoxygen[t][d].reshape(-1))
            temperature.append(tmptemperature[t][d].reshape(-1))
            nitrate.append(tmpnitrate[t][d].reshape(-1))
            chla.append(tmpchla[t][d].reshape(-1))
            ta.append(tmpta[t][d].reshape(-1))
            dic.append(tmpdic[t][d].reshape(-1))
            silicate.append(tmpsilicate[t][d].reshape(-1))
            hor_spd.append(tmphor_spd[t][d].reshape(-1))
            vrt_spd.append(tmpvrt_spd[t][d].reshape(-1))
            npp.append(tmpnpp[t][d].reshape(-1))
            phyc.append(tmpphyc[t][d].reshape(-1))
            ppco2.append(tmpppco2[t][d].reshape(-1))
            
            xx.append(tmpxx.reshape(-1))
            yy.append(tmpyy.reshape(-1))
            
            tmptid = np.ones(len(tmpppco2[t][d].reshape(-1)), dtype=int)*int(tmptime[t])
            # tmptid = [time.gmtime(t) for t in tmptid]
            tid.append(tmptid)
            year.append([time.gmtime(t).tm_year for t in tmptid])
            month.append([time.gmtime(t).tm_mon for t in tmptid])
            day.append([time.gmtime(t).tm_mday for t in tmptid])
            # year.append(np.ones(len(tmpppco2[t][d].reshape(-1)))*int(tmptime[t]))
            # month.append(np.ones(len(tmpppco2[t][d].reshape(-1)))*int(tmptime[t]))
            # day.append(np.ones(len(tmpppco2[t][d].reshape(-1)))*int(tmptime[t]))
            dyp.append(np.ones(len(tmpppco2[t][d].reshape(-1)))*tmpdyp[d])
    
    rho = gsw.rho_t_exact(salinity, temperature, dyp)/1000
    if (fid.variables['O3_c'].getncattr('units') == 'mmol C/m^3' or
            fid.variables['O3_c'].getncattr('units') == 'umol/kg'):
        dic = dic / rho
    rho=rho*1000
    
    df = pd.DataFrame(np.array([np.array(tid).reshape(-1),
                                np.array(year).reshape(-1),
                                np.array(month).reshape(-1),
                                np.array(day).reshape(-1),
                                np.array(xx).reshape(-1),
                                np.array(yy).reshape(-1),
                                np.array(dyp).reshape(-1),
                                np.array(pH).reshape(-1),
                                np.array(phosphate).reshape(-1),
                                np.array(salinity).reshape(-1),
                                np.array(oxygen).reshape(-1),
                                np.array(temperature).reshape(-1),
                                np.array(nitrate).reshape(-1),
                                np.array(chla).reshape(-1),
                                np.array(ta).reshape(-1),
                                np.array(dic).reshape(-1),
                                np.array(silicate).reshape(-1),
                                np.array(hor_spd).reshape(-1),
                                np.array(vrt_spd).reshape(-1),
                                np.array(npp).reshape(-1),
                                np.array(phyc).reshape(-1),
                                np.array(ppco2).reshape(-1),
                                np.array(rho).reshape(-1),]).T,
                      columns=['tid','year','month','day','xx','yy','dyp','pH',
                               'phosphate','salinity','oxygen','temperature',
                               'nitrate','chla','ta','dic','silicate','hor_spd',
                               'vrt_spd','npp','phyc','ppco2','Dens'])
    fid.close()
    return df

def readData(dataPath):
    fileList = glob.glob(dataPath + '*.nc')
    fileList.sort()
    dataOut = None
    for file in fileList:
        print('Now processing ' + file.split('/')[-1])
        if dataOut is None:
            dataOut = readFile(file)
        else:
            dataOut = dataOut.merge(readFile(file), how='outer')
    return dataOut
    
    # data = None
    # fid = nc.Dataset(dataFile + '*.nc')
    # if not fid is None:
    #     lat = fid.variables['lat'][:]
    #     lon = fid.variables['lon'][:]
    #     if ('CanESM5' in dataFile.split('_')[2]) or ('CMCC-ESM2' in dataFile.split('_')[2]) or ('CESM2' in dataFile.split('_')[2]) or ('GFDL' in dataFile.split('_')[2]) or ('NorESM' in dataFile.split('_')[2]) :
    #         date = timeConverter(fid.variables['time'][:], fid.variables['time'].units, 365)
    #     elif ('UKESM1' in dataFile.split('_')[2]) :
    #         date = timeConverter(fid.variables['time'][:], fid.variables['time'].units, 360)
    #     else:
    #         date = timeConverter(fid.variables['time'][:], fid.variables['time'].units)
    #     if 'IPSL' in dataFile.split('_')[2]:
    #         depth = fid.variables['olevel'][:]
    #     else:
    #         depth = fid.variables['lev'][:]
    #     exp = dataFile.split('_')[3]
    #     if exp=='piControl':
    #         data = fid.variables[variable][-600:][:][:][:]
    #         date = date[-600:]
    #     else:
    #         data = fid.variables[variable][(np.array(date > np.array(T[0])) & np.array(date < np.array(T[-1])))][:][:][:]
    #         date = date[(np.array(date > np.array(T[0])) & np.array(date < np.array(T[-1])))]
    #     if len(date)!=600 and dataFile.split('_')[3] == 'piControl':
    #         print('Time vector error!!')
    # #    time = [dt.datetime.strptime(t, '%Y-%m-%d').date() for t in time]
    #     time = []
    #     for t in date:
    #         time.append(dt.datetime.strptime(t, '%Y-%m-%d').year)# + (dt.datetime.strptime(t, '%Y-%m-%d').month-1)/12)
    #     m = min(time)
    #     time = [x - m for x in time]
    #     FV = fid.variables[variable].getncattr('_FillValue')
    #     MV = fid.variables[variable].getncattr('missing_value')
    #     data = xr.DataArray(
    #         data = data.squeeze(),
    #         coords=[time, lat, lon],
    #         dims=['time', 'lat', 'lon'],
    #         attrs=dict(
    #             var=variable,
    #             exp = exp,
    #             name=dataFile.split('_')[2],
    #             units=fid.variables[variable].units,
    #             realDepth = depth[0],
    #             depth=D
    #             )
    #         )
    #     fid.close()
    #     if len(data['time'].values) > len(np.unique(data['time'].values)):
    #         attrs = data.attrs
    #         data = data.groupby('time').mean()
    #         data.attrs = attrs
    #     # data = data.sel(depth=D , method='nearest')
    #     data = data.sortby(data['time'])
    #     data.values[data.values==FV] = None
    #     data.values[data.values==MV] = None
    #     data.values[data.values==1e99] = None
    # return data

#dataPath1='/Users/abom/Documents/data/Goldeneye model data/'
#dataPath2='/Users/abom/Documents/Arbeid/DemoCode_Cseep_python/Input_data_Cseep/'
#dataPath = 'Data/'
#dataPath = 'Data/'
dataPath1='input/external/'
if len(os.listdir(dataPath1)) == 0:
   dataPath1='input/netcdf/'
modData = readData(dataPath1)
pickle.dump(modData, open('modData.dat', 'wb'))
