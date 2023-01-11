#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 10:44:27 2022

@author: jean
"""

## This is for the demonstarion of anomaly criteria of the ACTOM toolbox. 
# 
# Input: data (model or obsevations) from the Goldeneye site including the parameters:
# dissolved inorganic carbon (DIC), total alkanlinity (AT), nutrinets (PO4
# and NO3), temperature, and salinity
#
# Dependency: runs the function 'minimize_DIC_var_Stemm.m'
#
#Ourput: 
# 1) natral C-variations predicted from non-co2 parameters as well as total (natural and seepage) C-variations 
# 2) plots total variations vs predicts natral C-variations
## load data 

from scipy.io import loadmat # to load matlab files
import numpy as np
import matplotlib.pyplot as plt
import pickle
import configparser
import PyCO2SYS as pyco2
from datetime import date

def readData(datafile):
    tmp = loadmat(datapathin + 'obs_data.mat')['data_o'] #this is the stemm 2018 cruise
    obs_data = {
        'Cruise': tmp['Cruise'][0][0][:].squeeze(),
        'DICmolkg': tmp['DICmolkg'][0][0][:].squeeze(),
        'Day': tmp['Day'][0][0][:].squeeze(),
        'Depth': tmp['Depth'][0][0][:].squeeze(),
        'Hour': tmp['Hour'][0][0][:].squeeze(),
        'Latitude': tmp['Latitude'][0][0][:].squeeze(),
        'Longitude': tmp['Longitude'][0][0][:].squeeze(),
        'Minute': tmp['Minute'][0][0][:].squeeze(),
        'Month': tmp['Month'][0][0][:].squeeze(),
        'Niskin': tmp['Niskin'][0][0][:].squeeze(),
        'NitrateM': tmp['NitrateM'][0][0][:].squeeze(),
        'NitriteM': tmp['NitriteM'][0][0][:].squeeze(),
        'PhosphateM': tmp['PhosphateM'][0][0][:].squeeze(),
        'Salinity': tmp['Salinity'][0][0][:].squeeze(),
        'SilicateM': tmp['SilicateM'][0][0][:].squeeze(),
        'Station': tmp['Station'][0][0][:].squeeze(),
        'TAmolkg': tmp['TAmolkg'][0][0][:].squeeze(),
        'TDNM': tmp['TDNM'][0][0][:].squeeze(),
        'Temperature': tmp['Temperature'][0][0][:].squeeze(),
        'Year': tmp['Year'][0][0][:].squeeze(),
        'pCO2atminsitu': tmp['pCO2atminsitu'][0][0][:].squeeze(),
        # 'pCO2_Merbach': tmp['pCO2_Merbach'],
        'pHat25': tmp['pHat25'][0][0][:].squeeze(),
        'pHinsitu': tmp['pHinsitu'][0][0][:].squeeze(),
        'Press_dbar': tmp['Press_dbar'][0][0][:].squeeze(),
        'Tpot': tmp['Tpot'][0][0][:].squeeze(),
        'Dens': tmp['Dens'][0][0][:].squeeze(),
        'Nitrateumolkg': tmp['Nitrateumolkg'][0][0][:].squeeze(),
        'Phosphateumolkg': tmp['Phosphateumolkg'][0][0][:].squeeze(),
        'Silicateumolkg': tmp['Silicateumolkg'][0][0][:].squeeze(),
        'nDIC': tmp['nDIC'][0][0][:].squeeze(),
        'MonthYear': tmp['MonthYear'][0][0][:].squeeze(),
        # 'Time': tmp['Time'][0][0][:].squeeze()
        }
    return obs_data

# Function only used to define plots sizes
def cm_to_inch(value):
    return value/2.54

def computeModCb(data, depth):
## Adaptation of minimize_DIC_var_Stemm to modData
    
        # data = data[data['dyp']>90]
        # data.loc[data.dyp>90]
    
    ## Constatnts and refeference values
    # Redfield ratios
    rCP=117 # ?13.7 #
    # annual increase of anthropogenic carbon in the North Sea (Omar et al, 2019 JGR)
    Cant=1.2 #? 0.3  
    # compute potential allkalinity (PTA) from TA and nitrate (if not available use phosphate) 
    data['PTA'] = data['ta'] + 1.36 * data['nitrate'] #Potential TA followindata = g Carter et al. (2014) BG    
    # slope and intercept of the PTA vs S relatioship based on all Baseline data. 
    PTA_S0 = -46    # ? 261 
    PTA_a = 67.5 # ? 7.4
    # Arbitrary chosen reference values. Here we use the mean values of the
    # historic baseline data
    S_refGoldeneye = 35.0994 #nanmean(data.Salinity_psu)
    PO4_refGoldeneye = 0.6241 #nanmean(data.PO4_molkg)
    NO3_refGoldeneye = 7.4414 #nanmean(data.NO3_molkg). Note: nanmean(data.NO3_NO2_molkg) = 7.8110! 
    AT_refGoldeneye = 2.3232e+03 # determined from PTA_a*S_refGoldeneye+PTA_S0
    
    # Use the same reference values for STEMM-CCS and historic baseline data
    S_refBaseline = S_refGoldeneye
    PO4_refBaseline = PO4_refGoldeneye
    NO3_refBaseline = NO3_refGoldeneye
    AT_refBaseline = AT_refGoldeneye
    AT_refBaselineGoldeneye = PTA_a * S_refBaseline+PTA_S0
    #Calculating S contribution in TA varaibility and eliminate it (two ways):
    # 1. by salinity normalization
    data['nPTA_nor'] = (( data['PTA'] - PTA_S0 ) * ( S_refBaseline / data['salinity'] )) + PTA_S0
    # 2. by difference. This is used in this implementation
    data['DS'] = ( data['salinity'] - S_refBaseline ) * PTA_a
    data['nPTA'] = data['PTA'] - data['DS']
    ## Calculate contributions of the natural DIC variability (DCnat) and eliminate each one of them from Cm 
    # DCbio_org contribution: formation/decay of organic matter
    data['Corg_P'] = rCP * ( data['phosphate'] - PO4_refBaseline )
    # DCbio_calc contribution: formation/dissolution of CaCO3
    data['Ccalc'] = 0.5 * ( data['nPTA'] - AT_refBaseline )
    #Eliminating the influence of ogrganic matter
    data['C_abio_P'] = data['dic'] - data['Corg_P']
    #Eliminating the influence of CaCO3
    data['C_acalc'] = data['dic'] - data['Ccalc']
    
    #Eliminating both organic matter and CaCO3
    data['Cseep_P'] = data['dic'] - data['Corg_P'] - data['Ccalc']
    # Compute and eliminate influence of salinity changes on DIC (two alternative ways):
    # 1. by salinity normalization
    data['nCseep_P_nor'] = (( data['Cseep_P']-PTA_S0 ) * S_refBaseline / data['salinity'] ) + PTA_S0
    #2. by using the slope PTA_a and salnity difference (see DS above). 
    data['nCseep_P'] = data['Cseep_P'] - data['DS']
    # Eliminating influence of anthropogenic carbon (Cant) differences
    # correct all data to 2010 assuminng 1.2 umolkg-1/yr increase (see Omar et al, 2019 JGR)
    data['Canth'] = ( data['year'] - 2010 ) * Cant
    
    data['Cb'] = data['nCseep_P'] - data['Canth']
    return data

def minimize_DIC_var_Stemm(data):
## Computes Baseline DIC (Cb) containing minimal spatiotemportal variability in the Goldeneye area
# Input: STEMM-CCS baseline OR monitoring data including the parameters:
# dissolved inorganic carbon (DIC), total alkanlinity (TA), nutrinets (PO4
# and NO3), temperature, and salinity
# Output: Out=[data.Cb data.nCseep_P_nor data.Corg_P data.Ccalc data.DS data.Canth data.Cseep_P]. NB: Cb= Cb+Cseep for monitoring data

    ## Constatnts and refeference values
    # Redfield ratios
    rCP=117 # ?13.7 #
    # annual increase of anthropogenic carbon inimport PyCO2SYS as pyco2 the North Sea (Omar et al, 2019 JGR)
    Cant=1.2 #? 0.3  
    # compute potential allkalinity (PTA) from TA and nitrate (if not available use phosphate) 
    data['PTA'] = data['TAmolkg'] + 1.36 * data['Nitrateumolkg'] #Potential TA followindata = g Carter et al. (2014) BG    
    # slope and intercept of the PTA vs S relatioship based on all Baseline data. 
    PTA_S0 = -46    # ? 261 
    PTA_a = 67.5 # ? 7.4
    # Arbitrary chosen reference values. Here we use the mean values of the
    # historic baseline data
    S_refGoldeneye = 35.0994 #nanmean(data.Salinity_psu)
    PO4_refGoldeneye = 0.6241 #nanmean(data.PO4_molkg)
    NO3_refGoldeneye = 7.4414 #nanmean(data.NO3_molkg). Note: nanmean(data.NO3_NO2_molkg) = 7.8110! 
    AT_refGoldeneye = 2.3232e+03 # determined from PTA_a*S_refGoldeneye+PTA_S0
    
    # Use the same reference values for STEMM-CCS and historic baseline data
    S_refBaseline = S_refGoldeneye
    PO4_refBaseline = PO4_refGoldeneye
    NO3_refBaseline = NO3_refGoldeneye
    AT_refBaseline = AT_refGoldeneye
    AT_refBaselineGoldeneye = PTA_a * S_refBaseline+PTA_S0
    #Calculating S contribution in TA varaibility and eliminate it (two ways):
    # 1. by salinity normalization
    data['nPTA_nor'] = (( data['PTA'] - PTA_S0 ) * ( S_refBaseline / data['Salinity'] )) + PTA_S0
    # 2. by difference. This is used in this implementation
    data['DS'] = ( data['Salinity'] - S_refBaseline ) * PTA_a
    data['nPTA'] = data['PTA'] - data['DS']
    ## Calculate contributions of the natural DIC variability (DCnat) and eliminate each one of them from Cm 
    # DCbio_org contribution: formation/decay of organic matter
    data['Corg_P'] = rCP * ( data['Phosphateumolkg'] - PO4_refBaseline )
    # DCbio_calc contribution: formation/dissolution of CaCO3
    data['Ccalc'] = 0.5 * ( data['nPTA'] - AT_refBaseline )
    #Eliminating the influence of ogrganic matter
    data['C_abio_P'] = data['DICmolkg'] - data['Corg_P']
    #Eliminating the influence of CaCO3
    data['C_acalc'] = data['DICmolkg'] - data['Ccalc']
    
    #Eliminating both organic matter and CaCO3
    data['Cseep_P'] = data['DICmolkg'] - data['Corg_P'] - data['Ccalc']
    # Compute and eliminate influence of salinity changes on DIC (two alternative ways):
    # 1. by salinity normalization
    data['nCseep_P_nor'] = (( data['Cseep_P']-PTA_S0 ) * S_refBaseline / data['Salinity'] ) + PTA_S0
    #2. by using the slope PTA_a and salnity difference (see DS above). 
    data['nCseep_P'] = data['Cseep_P'] - data['DS']
    # Eliminating influence of anthropogenic carbon (Cant) differences
    # correct all data to 2010 assuminng 1.2 umolkg-1/yr increase (see Omar et al, 2019 JGR)
    data['Canth'] = ( data['Year'] - 2010 ) * Cant
    
    data['Cb'] = data['nCseep_P'] - data['Canth']
    return data

def makePlot(xObs, yObs, xMod, yMod, p, err, month):
    
    x1 = min([min(xObs), min(xMod)])
    x2 = max([max(xObs), max(xMod)])
    X = [x1,x2]
    Y = np.polyval([0.76, -3.1],[x1,x2])
    # res = np.sqrt(res)*2
    
    ax = plt.figure(figsize=(cm_to_inch(15),cm_to_inch(10)), dpi=200)
    plt.grid(linestyle=':')
    plt.fill_between([x1,x2], Y-err, Y+err, color='y', alpha=.1) #std curves.
    # plt.errorbar(xMod, yMod, yerr=np.abs(yMod-xMod), xerr=np.abs(yMod-xMod),
    #              fmt = '.', color='black', alpha=.3)
    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    colors = plt.cm.jet(np.linspace(0, 1, 12))
    for n in range(12):
        plt.plot(xMod[month==n], yMod[month==n], '.', label=labels[n], color=colors[n], alpha=.3)
    ax.legend(loc='center right')
#    plt.plot(xObs, yObs, 'o', color='blue')
#    plt.plot(X,Y)
    plt.xlabel('Predicted natural C-variability (umolkg)')
    plt.ylabel('Observed C-variability (umolkg)')
    # plt.plot([min(X),max(X)], yfit, color='g')
    plt.savefig(datapathout+'Fig.png')
    plt.show()
    plt.close()

def detectLeak(modData, obsData, plot = True):
    modData.dropna(inplace=True)
    GE_CB = 2148 #from Omar et al 2021.
    xMod = modData['Corg_P'] + modData['Ccalc'] + modData['DS'] + modData['Canth']
    yMod = modData['dic'] - 1.2 * (modData['year'] - 2010) - GE_CB
    # yMod = modData['dic'] - 1.2 * (modData['year'] - 2010) - np.mean(modData['dic'] - 1.2 * (modData['year'] - 2010))
    p, res, _, _, _ = np.polyfit(xMod, yMod, 1, full=True)
    print('Correlation coeffiscient: ' + '{:.2f}'.format(np.corrcoef(xMod, yMod)[0,1]) + 
          ', slope: ' + '{:.2f}'.format(p[0]) + ', y-offset: ' + '{:.2f}'.format(p[1]))
          
    err = 24 #max(np.abs(yMod-xMod))
    DICcb=1.2*(float(date.today().year)-2010)+GE_CB
    DICmean=np.mean(modData['dic'])
    TAmean=np.mean(obsData['TAmolkg'])
    Densmean=np.mean(obsData['Dens'])

    with open(datapathout+'output.txt', 'w') as f:
        f.write('Threshold:' + str(err) +', Correlation coefficient: ' + '{:.2f}'.format(np.corrcoef(xMod, yMod)[0,1]) + 
          ', slope: ' + '{:.2f}'.format(p[0]) + ', y-offset: ' + '{:.2f}'.format(p[1]))


    config=configparser.ConfigParser(allow_no_value=True)
    config.read('/external/settings/data.ini')
    config.add_section('CSEEP')
    config.set('CSEEP', 'Threshold_DIC_micromol-kg', str(err))
    DICout2=err*1024*0.0440095/1000000 
    config.set('CSEEP', 'Threshold', str(DICout2))
    kwargs = dict(
        par1 = DICmean-(0.5*err),
        par2 = TAmean,
        par1_type = 2,
        par2_type = 1,
        )
    outDIC=pyco2.sys(**kwargs)
    pHout1=outDIC["pH_total"]
    kwargs = dict(
        par1 = DICmean+(0.5*err),
        par2 = TAmean,
        par1_type = 2,
        par2_type = 1,
        )
    outDIC=pyco2.sys(**kwargs)    
    pHout2=outDIC["pH_total"]    
    pHout=abs(pHout2-pHout1)
    config.set('CSEEP', 'threshold-ph', str(pHout))
#    config.set('CSEEP', 'DICcb', str(DICcb))
    config.set('CSEEP', 'DIC-mean', str(DICmean))
    config.set('CSEEP', 'TA-mean', str(TAmean))
    config.set('CSEEP', 'dens-mean', str(Densmean))               
    with open('/external/settings/data.ini', 'w') as configfile:
        config.write(configfile)
    
    xObs = obsData['Corg_P'] + obsData['Ccalc'] + obsData['DS'] + obsData['Canth']
    yObs = obsData['DICmolkg'] - 1.2 * (obsData['Year'] - 2010) - GE_CB
    yFit = np.polyval([0.8,-3.1],xObs)
    if sum(yObs>yFit+err)>0:
        print('Possible anomalies found!')
    else:
        print('No leak found')
    if plot:
        makePlot(xObs, yObs, xMod, yMod, p, err, modData['month'])


# datapath = '/media/Titan/Boulot/2021_NORCE/Code/Python/Abdir/Docker/Data/'
#datapath = 'Data/'
#datapath = '/Users/abom/Documents/Arbeid/DemoCode_Cseep_python/Input_data_Cseep/'
datapathin='input/'
datapathout='output/'
depth = 'bottom'
obsData = readData([datapathin +'obs_data.mat']) #this is the stemm 2018 cruise
modData = pickle.load(open('modData.dat', 'rb'))
if depth == 'bottom':
    modData = modData.loc[modData['dyp']>90]
modData = computeModCb(modData, depth)
data = minimize_DIC_var_Stemm(obsData)
# [data2.Cb data2.nCseep_P_nor data2.Corg_P data2.Ccalc data2.DS data2.Canth data2.Cseep_P]
GE_CB = 2148 #from Omar et al 2021.
# Plotting the first figure
detectLeak(modData, obsData)
