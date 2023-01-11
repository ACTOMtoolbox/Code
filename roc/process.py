import glob
import csv
import os
import numpy as np
import configparser
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

# define the true objective function
def objective(x, a, b, c):
#  return a * np.exp(-b * x) + c
#  return  a*(x**b)+c
  return a * np.log(b*x) + c

plt.style.use('seaborn-whitegrid')

# print(glob.glob("/tmp/*.dat"))
files=glob.glob("/external/output/*")

# --- we want an input of intime & units

config=configparser.ConfigParser(allow_no_value=True)
config.read('/external/settings/data.ini')

intime = float(config['RateOfChange']['Sampling_Frequency']) # Sensor measurement frequency
intimeunit = config['RateOfChange']['Sampling_Frequency_Units'] # Sensor measurement frequency units

e="/external/output/"
h = "output.dat"
f = e+h

g = open(f, "w")
writer = csv.writer(g)
i = ['name', 'sample time', 'units', 'anomaly criteria']
writer.writerow(i)

outROC=1

for a in files:
    if "png" in a or "dat" in a:
        continue
    with open(a, "r") as csv_file:
        reader = csv.reader(csv_file, delimiter=' ', skipinitialspace=True)
        Yin = list(zip(*reader))[3]
        Yl = [float(i) for i in Yin]
        csv_file.close()
    with open(a, "r") as csv_file:
        reader = csv.reader(csv_file, delimiter=' ', skipinitialspace=True)
        Xin = list(zip(*reader))[1]
        Xin_list = [s.replace(".00", "") for s in Xin]
        integer_map = map(int, Xin_list)
        Xl = list(integer_map)
        XL  = np.array(Xl)
        XL = np.delete(XL, np.where(XL >= 18000))
        Xlength=len(XL)
        Yl = Yl[:Xlength]
        YL  = np.array(Yl)
        limit = 18000

        if intimeunit[0] == "m":
            if intimeunit[1] == "i":
                XL = XL/60
                limit=limit/60
            if intimeunit[1] == "o":
                XL = XL/2629800
                limit=limit/2629800
        if intimeunit[0] == "h":
            XL = XL/3600
            limit=limit/3600
        if intimeunit[0] == "d":
            XL = XL/86400
            limit=limit/86400
        if intimeunit[0] == "w":
            XL = XL/604800
            limit=limit/604800
        if intimeunit[0] == "y":
            XL = XL/31557600
            limit=limit/31557600

        if intime > max(XL):
          print("Error: the baseline data is too small for the chosen sampling rate")
          exit()

        csv_file.close()
         
        plt.plot(XL,YL,'go',markersize=1)
        
        popt, _ = curve_fit(objective, XL, YL)
        x_line = np.arange(min(XL), max(XL), 1)
        y_line = objective(x_line, *popt)
        outROC=objective(intime, *popt)
        
#        d=1
#        for x in XL:
#          if intime > x:
#            d=d+1
#        dx = XL[d]-XL[d-1]
#        dy = YL[d]-YL[d-1]
#        outROC=YL[d-1]+((dy/dx)*(intime-XL[d-1]))
#        outROC=trendlineeq(intime)

        plt.plot(x_line, y_line,'b-')
        plt.plot(intime,outROC,'ro')
        output=str(round(outROC,6))+" at "+str(intime)+" "+intimeunit+" sampling rate"
        plt.annotate(output, (intime, outROC))

        title="Anomaly Criteria from "+a[17:]
        plt.title(title)
        k="Time ("
        l=k+intimeunit+')'
        plt.xlabel(l)
        plt.xlim((0, limit))
        if "pH" in a:
            plt.ylabel("pH change")
        else:
            plt.ylabel("DIC change")
        plt.grid(visible=None)
        b=".png"
        c=a+b
        print("-- PRINTING OUTPUT FROM %s --" % (a[17:]))
        plt.savefig(c, dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            metadata=None)
        plt.close()
        if "pH" in a:
            config.set('RateOfChange', 'Threshold-pH', str(outROC))
        elif "kgm3" in a:
            config.set('RateOfChange', 'Threshold', str(outROC))
        else:
            config.set('RateOfChange', 'Threshold-DIC', str(outROC))              
        with open('/external/settings/data.ini', 'w') as configfile:
            config.write(configfile)
        j=[a[17:], intime, intimeunit, outROC]
        writer.writerow(j)


g.close()
