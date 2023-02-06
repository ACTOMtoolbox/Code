import pandas as pd
import numpy as np
import netCDF4 as nc
import jinja2
import os
import glob
import csv
import sys
from termcolor import colored, cprint
import configparser

env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=''))

filelocations="input/"

print("\n******************************************************")
print("\nProcessing Outputs\n")
print("******************************************************\n")

os.system('cp -r storage/files/ input/')

# -- Checking for tools and printing appropriate documentation-- #

tool_list=''
tool_list2=''
OUTPUTS=''
REF=[]

#print('checking for the Sub-Surface CO₂ Seep and Rate Mapping Tool:')
#if os.path.isdir('/tmp/geo'):
#  print('found')
#  tool_list=tool_list+'<b><li>The Sub-Surface CO<sub>2</sub> Seep and Rate Mapping Tool</li></b>'
#  tool_list2=tool_list2+'<b>The Sub-Surface CO<sub>2</sub> Seep and Rate Mapping Tool</b>'
#else:
#  cprint("Missing!", 'red', attrs=['bold'])
#  tool_list=tool_list+'<font color="red"><li>The Sub-Surface CO<sub>2</sub> Seep and Rate Mapping Tool</li></font>'
#  tool_list2=tool_list2+'<font color="red">The Sub-Surface CO<sub>2</sub> Seep and Rate Mapping Tool</font>'

print('checking for the Tracer Transport Model:')
if os.path.isdir(filelocations+'Advdiff'):
  print('found')
  tool_list=tool_list+'<b><li>The Tracer Transport Model</li></b>'
  tool_list2=tool_list2+'<b><br><br>The Tracer Transport Model</b>'
else:
  cprint("Missing!", 'red', attrs=['bold'])
  tool_list=tool_list+'<font color="red"><li>The Tracer Transport Model</li></font>'
  tool_list2=tool_list2+'<font color="red"><br><br>The Tracer Transport Model</font>'

print('checking for Cseep:')
if os.path.isdir(filelocations+'Cseep'):
  print('found')
  tool_list=tool_list+'<b><li>C<sub>seep</sub></li></b>'
  tool_list2=tool_list2+'<b><br><br>C<sub>seep</sub></b>'
else:
  cprint("Missing!", 'red', attrs=['bold'])
  tool_list=tool_list+'<font color="red"><li>C<sub>seep</sub></li></font>'
  tool_list2=tool_list2+'<font color="red"><br><br>C<sub>seep</sub></font>'
  
print('checking for the Rate of Change Anomaly Criteria:')
if os.path.isdir(filelocations+'ROC'):
  print('found')
  tool_list=tool_list+'<b><li>The Rate of Change Anomaly Criteria</li></b>'
  tool_list2=tool_list2+'<b><br><br>The Rate of Change Anomaly Criteria</b>'
else:
  cprint("Missing!", 'red', attrs=['bold'])
  tool_list=tool_list+'<font color="red"><li>The Rate of Change Anomaly Criteria</li></font>'
  tool_list2=tool_list2+'<font color="red"><br><br>The Rate of Change Anomaly CriteriaTool</font>'
  
print('checking for the Carbonate System:')
if os.path.isdir(filelocations+'carbon'):
  print('found')
  tool_list=tool_list+'<b><li>The Carbonate System</li></b>'
  tool_list2=tool_list2+'<b><br><br>The Carbonate System</b>'
else:
  cprint("Missing!", 'red', attrs=['bold'])
  tool_list=tool_list+'<font color="red"><li>The Carbonate System</li></font>'
  tool_list2=tool_list2+'<font color="red"><br><br>The Carbonate System</font>'

print('checking for Impact Analysis:')
if os.path.isdir(filelocations+'impacts'):
  print('found')
  tool_list=tool_list+'<b><li>Impact Analysis</li></b>'
  tool_list2=tool_list2+'<b><br><br>Impact Analysis</b>'
else:
  cprint("Missing!", 'red', attrs=['bold'])
  tool_list=tool_list+'<font color="red"><li>Impact Analysis</li></font>'
  tool_list2=tool_list2+'<font color="red"><br><br>Impact Analysis</font>'
  
print('checking for the Optimal Cover Tool:')
if os.path.isdir(filelocations+'OptCover'):
  print('found')
  tool_list=tool_list+'<b><li>The Optimal Cover Tool</li></b>'
  tool_list2=tool_list2+'<b><br><br>The Optimal Cover Tool</b>'
else:
  cprint("Missing!", 'red', attrs=['bold'])
  tool_list=tool_list+'<font color="red"><li>The Optimal Cover Tool</li></font>'
  tool_list2=tool_list2+'<font color="red"><br><br>The Optimal Cover Tool</font>'

template = env.get_template('storage/template.html')
a = template.render(Tool_List=tool_list,Documentation=': Monitoring Plan',OUTPUTS='{{ OUTPUTS }}',RISKMAP='{{ RISKMAP }}',HYDRODYNAMICS='{{ HYDRODYNAMICS }}',REF='{{ REF }}')
with open('storage/template-wt.html', 'w') as f:
  f.write(a)

figure_count = 1

# -- fill template with text if geo data from AdvDiff exisit --#
if os.path.isdir(filelocations+'Advdiff'):
  geo_loc = os.listdir(filelocations+'Advdiff/output')
  geolen=len(geo_loc)
  
  sourcelen = 0
  for x in geo_loc:
    if x[0:6]=='fields' and x[7] !='g':
      sourcelen=sourcelen+1
      
  RISKMAP2="<h2>Site Specific Input</h2>"
  
  RISKMAP="The riskmap for the study area is built upon "
  
  config=configparser.ConfigParser(allow_no_value=True)
  #config.read('/external/settings/Advdiff/AdvDiff.ini')
  config.read(filelocations+'Advdiff/AdvDiff.ini')
  source = config['setup']['get_source_from'].lower()

  RISKMAP=RISKMAP+source+" from the file"
  
  sourcesplit=source.split(" ")
  
  addsource=[]
  
  if len(sourcesplit) > 1:
    sourcesplit.remove('and')
  sourcesplit=list(map(lambda x: x.replace('riskmap', 'riskmaps'),sourcesplit))

  for x in sourcesplit:
  
    if x == "sources":
       y = "source_file"
    if x == "riskmaps":
       y = "risk_file"
    if x == "wells":
       y = "well_file"
     
    addsource.append(config[x][y])
    
  if len(addsource)==1:
     addsourceout=': '+addsource[0]
  else:
     addsourceout='s: '+' and '.join(addsource)
    
  RISKMAP=RISKMAP+addsourceout
  figure_count=figure_count+1
  
  RISKMAP=RISKMAP+", with the total of "+str(sourcelen)+" sources"
  
  RISKMAP2=RISKMAP2+RISKMAP

  RISKMAP=RISKMAP+", as shown on Figure "+str(figure_count)+"."
    
  imageRISKMAP="<p style=\"text-align:center;\"><img src=\'Advdiff/Figures/sources_local_grids_velocity.png\' height=\"400\"></p>"
  
  imageRISKMAP=imageRISKMAP+"<p style=\"text-align:center;\">Figure "+str(figure_count)+". Source locations (black circles) and local grids (dashed line boxes) superimposed on the u-velocity field (input data enclosed in red boundary, and extrapolated beyond).</p>"
  
  RISKMAP=RISKMAP+imageRISKMAP

  HYDRODYNAMICS="The local velocities for the region are collected from the file "
  source = config['velocity']['velocity_file']
  
  OUTPUT2=RISKMAP2+". The local velocities for the region are also collected from the file "+source+", with the example of the 'u' velocity shown in Figure "+str(figure_count)+"."
  
  HYDRODYNAMICS=HYDRODYNAMICS+source+", with the example of the 'u' velocity shown in Figure "+str(figure_count)+" above."
  
  HYDRODYNAMICS2="<br><br>....."
  
  OUTPUT2=OUTPUT2+imageRISKMAP+"<br><h2>Advection-Diffusion Module (Tracer Transport Simulator)</h2>"

  figstats=''
  if config['visualizer']['plot_stats']:
    figure_count=figure_count+1
    AD1_figs=figure_count
    stats="In Figure "+str(figure_count)+". we give an example of the worst case impact region from the "+str(sourcelen)+ " sources with equal weight of 1.0 units/m<sup>2</sup>s over the full simulation covering various tidal cycles. "
    
    figstats=figstats+'<p style=\"text-align:center\"><img src=\'Advdiff/Figures/statistics_max.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. Map of the largest impacted regions based on the sources provided and a release of 1.0 units/m<sup>2</sup>s, with the release sources as cricles and position of probes as stars.</p>'
  
  if config['visualizer']['plot_global_fields']:
    figure_count=figure_count+1
    if 'stats' not in locals():
      stats="With equal weight of 1.0 units/m<sup>2</sup>s time"
    else:
      stats=stats+"Time"
    stats=stats+" series for the fluctuating concentrations are provided (when viewed in html format) in Figure "+str(figure_count)+". "
    figstats=figstats+'<p style=\"text-align:center\"><img src=\'Advdiff/Figures/Global_field/Cmovie.gif\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. Time series for the fluctuating concentrations (when viewed in html format) based on the sources provided and a release of 1.0 units/m<sup>2</sup>s, with the release sources as cricles and position of probes as stars.</p>'

  if config['visualizer']['plot_probes']:
    if 'stats' not in locals():
      stats="With equal weight of 1.0 units/m<sup>2</sup>s probes are"
    else:
      stats=stats+"Probes are also"
    figure_count=figure_count+1
    stats=stats+" provided at various locations, such as potential monitoring stations can detect varying levels of fluctuations over time as shown in Figure "+str(figure_count)+". "
    
    figstats=figstats+'<p style=\"text-align:center\"><img src=\'Advdiff/Figures/probe_total_C.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. Time series for detecting fluctuation carbon levels at each probe based on the sources provided and a release of 1.0 units/m<sup>2</sup>s.</p>'
  
  with open('storage/AdvDiff.html') as f:
    add = f.read()
    
  OUTPUT2=OUTPUT2+stats+figstats
    
  OUTPUTS=OUTPUTS+add+stats+figstats

with open('storage/AC.html', 'r') as f:
  AC = f.read()    
OUTPUTS=OUTPUTS+AC
  
if os.path.isdir(filelocations+'Cseep'):    
  with open('storage/CSEEP.html', 'r') as f:
    CSEEP = f.read()  
  OUTPUTS=OUTPUTS+CSEEP
  
  OUTPUT2=OUTPUT2+"<h2>Anomaly Criteria Identification</h2><h3><i>C<sub>seep</sub></i></h3>"

  CSEEPOUT='With the input data provided to <i>C<sub>seep</sub></i>, a detection threshold of +/- '
  with open(filelocations+'Cseep/output/output.txt', 'r') as f:
    CSEEP = f.read()
  CSEEP=CSEEP.split(',', -1)
  figure_count=figure_count+1
  CSEEPOUT=CSEEPOUT+CSEEP[0].split(":",1)[1]+' (μmol/kg) is found, as shown in Figure '+str(figure_count)+'.'
  
  CSEEPOUT=CSEEPOUT+'<p style=\"text-align:center\"><img src=\'Cseep/output/Fig.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. <i>C<sub>seep</sub></i> output showing the observed variability in <i>C<sub>b</sub></i> against the predicted natural variability in <i>C<sub>b</sub></i> for the given data, The range is highlighted, where anything outside of this area of +/- '+CSEEP[0].split(":",1)[1]+' (μmol/kg) could potentially be an anomaly.</p>'

  OUTPUTS=OUTPUTS+CSEEPOUT
  OUTPUT2=OUTPUT2+CSEEPOUT+"<h3>Rate of Change Anomaly Criteria</h3>"
  
# -- Rate of Change Anomaly start -- #

if os.path.isdir(filelocations+'ROC'):

# -- fill template with text if Rate of Change Anomaly results exisit --#

  with open('storage/ROC.html', 'r') as f:
    ROC = f.read()

  col_list = ["name", "sample time", "units", "anomaly criteria"]
  roc = pd.read_csv("input/ROC/output/output.dat", usecols=col_list)
  roc_loc=list(roc["name"])
  roc_sample_time=list(roc["sample time"])
  roc_units=list(roc["units"])
  roc_anomaly_criteria=list(roc["anomaly criteria"])

  index_min = min(range(len(roc_loc)), key=roc_loc.__getitem__)

  figure_count=figure_count+1
  ROC2='From the input baseline data and set monitoring frequency of measurements of '+str(roc_sample_time[index_min])+' '+roc_units[index_min]+' DIC changes larger than '+str(round(roc_anomaly_criteria[index_min],2))+' (μmol/kg) are predicted to potentially be an anomaly that would require further analysis to determine if a potential release has occurred as shown in Figure '+str(figure_count)+'.'
  
  ROC2=ROC2+'<p style=\"text-align:center\"><img src=\'ROC/output/'+roc_loc[index_min]+'.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. The rate of change anomaly criteria tool determining that is a monitoring sample is taken every '+str(roc_sample_time[index_min])+' '+roc_units[index_min]+' a change in DIC of +/- '+str(round(roc_anomaly_criteria[index_min],2))+' (μmol/kg) is predicted to potentially be an anomaly (red dot).'

  OUTPUTS=OUTPUTS+ROC+ROC2
  OUTPUT2=OUTPUT2+ROC2

# -- The Carbonate System -- #


if os.path.isdir(filelocations+'carbon'):
  config=configparser.ConfigParser(allow_no_value=True)
  #config.read('/external/settings/data.ini')
  config.read(filelocations+'data.ini')
  rate = config['General']['rate']
  rate_units = config['General']['rate-units']
  rate = rate+' '+rate_units
  dic_mean = config['CSEEP']['dic-mean']
  ta_mean = config['CSEEP']['ta-mean']
  carbon='<h3>The Carbonate System</h3>The fields of tracer concentration are scaled according to the desired range of rates as set by the user. This means different release rates may be predicted without having to re-run the tracer transport simulator. In this case the leakage rate chosen was '+rate+'. Taking Total Alkalinity (<i>TA</i>): '+ta_mean[0:6]+' μmol/kg, and theoretical baseline DIC concentration (<i>C<sub>b</sub></i>) '+dic_mean[0:6]+' μmol/kg from <i>C<sub>seep</sub></i>, along with the scaled concentration from the tracer transport simulator, we can predict the DIC from the release. A carbonate system tool then converts these DIC and TA values into pH or pCO<sub>2</sub>, along with changes from background values based on CO2SYS (<a href="#References">Humphreys et al. 2022</a>).   The carbonate system also converts the thresholds provided by the anomaly criteria tools.'

  OUTPUTS=OUTPUTS+carbon

# -- Impacts -- #

if os.path.isdir(filelocations+'impacts'):
  with open('storage/Impacts.html', 'r') as f:
    impacts = f.read()

  OUTPUTS=OUTPUTS+impacts+'<br><br>'
  
  config=configparser.ConfigParser(allow_no_value=True)
  #config.read('/external/settings/Advdiff/AdvDiff.ini')
  config.read(filelocations+'data.ini')
  
  with open(filelocations+'impacts/output/Impact.txt') as f:
    lines = f.readlines()
  
  sub1 = [sub1.replace('Max Change in pH from C fields_global.nc - ', '') for sub1 in lines]
  sub2 = [sub2.replace('Impact-Area fields_global.nc pH from C ', '') for sub2 in sub1]
  sub3 = [sub3.replace('\n', '') for sub3 in sub2]

  impacts15='Setting each release to a rate of '+rate+', we find that the maximum pH change predicted by the combined leakages is '+sub3[0][0:5]+'. '
  
  sub4 = sub3[1:]
  
  
  if config.has_option('General', 'threshold-ph'):
    Thres=((config['General']['threshold-ph']).split(","))
    Thres=sorted(Thres, key=float, reverse=True)
       
  
  
  Impactarea=''
  impacts2=''
  if config.has_option('General', 'threshold-ph'):
    impacts2='With pH change thresholds based on user set values'
    init=0
    Thres1=((config['General']['threshold-ph']).split(","))
    Thres=[Thres.replace(' ', '') for Thres in Thres1]
    Thres=sorted(Thres, key=float, reverse=True)
    impacts2=impacts2+': '
    if len(Thres)>1:
      for x in range(len(Thres)):
            
        for y in range(len(sub4)):
          text=sub4[y]
          text=text.split(' ')
          if text[0]==Thres[x]:
            outtext=text[1]
            outtext=outtext.split('.')
            if len(sub4)>1 and len(Thres)>1:
              if x+1<len(Thres):
                Impactarea=Impactarea+outtext[0]+', '
              else:
                Impactarea=Impactarea[:-2]+' and '+outtext[0]+' m<sup>2</sup>'
            else:
              Impactarea=Impactarea+outtext[0]+' m<sup>2</sup>'
      
        if x+1<len(Thres):
          impacts2=impacts2+Thres[x][0:5]+', '
        else:
          impacts2=impacts2[:-2]+' and '+Thres[x][0:5]
    else:
      impacts2=impacts2+Thres[0][0:5]
      
      for y in range(len(sub4)):
        text=sub4[y]
        text=text.split(' ')
        if text[0]==Thres[0]:
          outtext=text[1]
          outtext=outtext.split('.')
          Impactarea=Impactarea+outtext[0]+' m<sup>2</sup>'
          
    if init==0:
       impacts2=impacts2+', the total area impacted is '+Impactarea
    else:
       impacts2=impacts2+', where the total area impacted is '+Impactarea
    if len(Thres)>1:
      impacts2=impacts2+', respectively. '
    else:
      impacts2=impacts2+'. '

  Impactarea='' 
  if config.has_option('CSEEP', 'threshold-ph'):
    Thres=((config['CSEEP']['threshold-ph']).split(","))
    Thres=sorted(Thres, key=float, reverse=True)
    if impacts2=='':
      impacts2=impacts2+'With pH change thresholds from <i>C<sub>seep</sub></i>'
      init=0
    else:
      impacts2=impacts2+'Along with pH change thresholds from <i>C<sub>seep</sub></i>'
      init=1
    impacts2=impacts2+': '
    if len(Thres)>1:
      for x in range(len(Thres)):
      
        for y in range(len(sub4)):
          text=sub4[y]
          text=text.split(' ')
          
          print(text[0])
          print(Thres[x])
          if text[0]==Thres[x]:
            outtext=text[1]
            outtext=outtext.split('.')
            if len(sub4)>1 and len(Thres)>1:
              if x+1<len(Thres):
                Impactarea=Impactarea+outtext[0]+', '
              else:
                Impactarea=Impactarea[:-2]+' and '+outtext[0]+' m<sup>2</sup>'
            else:
              Impactarea=Impactarea+outtext[0]+' m<sup>2</sup>'
      
        if x+1<len(Thres):
          impacts2=impacts2+Thres[x][0:5]+', '
        else:
          impacts2=impacts2[:-2]+' and '+Thres[x][0:5]
    else:
      impacts2=impacts2+Thres[0][0:5]
      
      for y in range(len(sub4)):
        text=sub4[y]
        text=text.split(' ')
        if text[0]==Thres[0]:
          outtext=text[1]
          outtext=outtext.split('.')
          Impactarea=Impactarea+outtext[0]+' m<sup>2</sup>'
    
    if init==0:
       impacts2=impacts2+', the total area impacted is '+Impactarea
    else:
       impacts2=impacts2+', where the total area impacted is '+Impactarea
    if len(Thres)>1:
      impacts2=impacts2+' respectively. '
    else:
      impacts2=impacts2+'. '
      
  Impactarea=''
  if config.has_option('RateOfChange', 'threshold-ph'):
    Thres=((config['RateOfChange']['threshold-ph']).split(","))
    Thres=sorted(Thres, key=float, reverse=True)
    if impacts2=='':
      impacts2=impacts2+'With pH change thresholds the Rate of Change Anomaly Criteria'
      init=0
    else:  
      impacts2=impacts2+'And pH change thresholds from the Rate of Change Anomaly Criteria'
      init=1
    impacts2=impacts2+': '
    if len(Thres)>1:
      for x in range(len(Thres)):
      
        for y in range(len(sub4)):
          text=sub4[y]
          text=text.split(' ')
          if text[0]==Thres[x]:
            outtext=text[1]
            outtext=outtext.split('.')
            if len(sub4)>1 and len(Thres)>1:
              if x+1<len(Thres):
                Impactarea=Impactarea+outtext[0]+', '
              else:
                Impactarea=Impactarea[:-2]+' and '+outtext[0]+' m<sup>2</sup>'
            else:
              Impactarea=Impactarea+outtext[0]+' m<sup>2</sup>'

        if x+1<len(Thres):
          impacts2=impacts2+Thres[x][0:5]+', '
        else:
          impacts2=impacts2[:-2]+' and '+Thres[x][0:5]
    else:
      impacts2=impacts2+Thres[0][0:5]
      
      for y in range(len(sub4)):
        text=sub4[y]
        text=text.split(' ')
        if text[0]==Thres[0]:
          outtext=text[1]
          outtext=outtext.split('.')
          Impactarea=Impactarea+outtext[0]+' m<sup>2</sup>'
      
    if init==0:
       impacts2=impacts2+', the total area impacted is '+Impactarea
    else:
       impacts2=impacts2+', where the total area impacted is '+Impactarea
    if len(Thres)>1:
      impacts2=impacts2+' respectively. '
    else:
      impacts2=impacts2+'. '
  
  rate = config['General']['rate']
  rateno=rate
  rate_units = config['General']['rate-units']
  rate = rate+' '+rate_units

  impactstats=''

  figure_count=figure_count+1
  impacts2=impacts2+'<br><br>The map of the largest impacted regions from the '+str(sourcelen)+ " sources in terms of pH changes over the full simulation are shown in Figure "+str(figure_count)+", "
  impactstats=impactstats+'<p style=\"text-align:center\"><img src=\'impacts/output/statistics_globalpH from max.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. Map of the largest impacted regions as pH changes at a release of '+rate+'.</p>'
  
  figure_count=figure_count+1
  impacts2=impacts2+'with the impact exceeding each threshold shown in Figure '+str(figure_count)+'. '
  impactstats=impactstats+'<p style=\"text-align:center\"><img src=\'impacts/output/fields_globalpH from CThresholds.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. Map of the largest impacted regions exceeding each pH change threshold at a release rate of '+rate+'.</p>'

  figure_count=figure_count+1
  impacts2=impacts2+' The time series for the fluctuating pH changes are also provided (when viewed in html format) in Figure '+str(figure_count)
  impactstats=impactstats+'<p style=\"text-align:center\"><img src=\'impacts/output/fields_globalpH from C.gif\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. Time series for the fluctuating pH changes (when viewed in html format) based on the sources provided and a release rate of '+rate+'.</p>'
  
  figure_count=figure_count+1
  impacts2=impacts2+' and probes provided at various locations, such as potential monitoring stations can detect varying levels of fluctuations from the releases at '+rate+' over time as shown in Figure '+str(figure_count)+'. '
  
  impactstats=impactstats+'<p style=\"text-align:center\"><img src=\'impacts/output/probe_cumulativepH from C.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. Time series for detecting fluctuation carbon levels at each probe based on the sources provided and a release of '+rate+'.</p>'

  OUTPUTS=OUTPUTS+impacts15+impacts2+impactstats
  
  OUTPUT2=OUTPUT2+"<h2>Impact Potential</h2>"+impacts2+impactstats

# -- The Optimal Cover Tool -- #


if os.path.isdir(filelocations+'OptCover'):

# -- fill template with text if Optimal Cover Tool exisits --#

  with open('storage/OptCover.html', 'r') as f:
    OptCover = f.read()  
    
  config=configparser.ConfigParser(allow_no_value=True)
  config.read(filelocations+'data.ini')
  
  numberofsensors=0
  OptCover2=''
  OptCover3=''


  
  if config.has_option('General', 'threshold-ph'):
    file1=glob.glob(filelocations+'OptCover/output/USER*.nc')
    file1.sort(reverse=True)
    Thres=((config['General']['threshold-ph']).split(","))
    Thres=sorted(Thres, key=float, reverse=True)
    for x in range(len(file1)):
      ncfile = nc.Dataset(file1[x])
      numberofsensors=ncfile.dimensions['sensor'].size
      if x == 0:
        OptCover2=OptCover2+'With a user set pH change threshold of '+ Thres[x][:5]
      else:
        OptCover2=OptCover2+'With another user set pH change threshold of '+ Thres[x][:5]
      figure_count=figure_count+1
      OptCover2=OptCover2+', Figure '+str(figure_count)+'. shows the best locations to place sensors (and the areas that they cover) to maximise the detection, using a minimum number of sensors: '+str(numberofsensors)+'. '
      
      OptCover3=OptCover3+'<p style=\"text-align:center\"><img src=\''+file1[x][6:-19]+'map2.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. The best locations to place sensors (and the areas that they cover) to maximise detection from a user set pH change threshold of '+Thres[x][:5]+', using a minimum number of sensors: '+str(numberofsensors)+'.</p>'
      
  if config.has_option('RateOfChange', 'threshold-ph'):
  
    file1=glob.glob(filelocations+'OptCover/output/ROC*.nc')
    file1.sort(reverse=True)
    Thres=((config['RateOfChange']['threshold-ph']).split(","))
    Thres=sorted(Thres, key=float, reverse=True)
    for x in range(len(file1)):
      ncfile = nc.Dataset(file1[x])
      numberofsensors=ncfile.dimensions['sensor'].size
      if x == 0:
        if OptCover2=='':       
          OptCover2=OptCover2+'With a pH change threshold from the Rate of Change Anomaly Criteria of '+ Thres[x][:5]
        else:
          OptCover2=OptCover2+'Finally looking at a pH change threshold from the Rate of Change Anomaly Criteria of '+ Thres[x][:5]
      else:
        OptCover2=OptCover2+'With another pH change threshold from the Rate of Change Anomaly Criteria of '+ Thres[x][:5]
      figure_count=figure_count+1
      OptCover2=OptCover2+', Figure '+str(figure_count)+'. shows the best locations to place sensors (and the areas that they cover) to maximise the detection, using a minimum number of sensors: '+str(numberofsensors)+'. '
      OptCover3=OptCover3+'<p style=\"text-align:center\"><img src=\''+file1[x][6:-19]+'map2.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. The best locations to place sensors (and the areas that they cover) to maximise detection from a Rate of Change Anomaly Criteria pH change threshold of '+Thres[x][:5]+', using a minimum number of sensors: '+str(numberofsensors)+'.</p>' 
     
  if config.has_option('CSEEP', 'threshold-ph'):
    
    file1=glob.glob(filelocations+'OptCover/output/CSEEP*.nc')
    file1.sort(reverse=True)
    Thres=((config['CSEEP']['threshold-ph']).split(","))
    Thres=sorted(Thres, key=float, reverse=True)
    for x in range(len(file1)):
      ncfile = nc.Dataset(file1[x])
      numberofsensors=ncfile.dimensions['sensor'].size
      if x == 0:
        if OptCover2=='':       
          OptCover2=OptCover2+'With a pH change threshold from <i>C<sub>seep</sub></i> of '+ Thres[x][:5]
        else:
          OptCover2=OptCover2+'Then looking at a pH change threshold from <i>C<sub>seep</sub></i> of '+ Thres[x][:5]
      else:
        OptCover2=OptCover2+'With another pH change threshold from <i>C<sub>seep</sub></i> of '+ Thres[x][:5]
      figure_count=figure_count+1
      OptCover2=OptCover2+', Figure '+str(figure_count)+'. shows the best locations to place sensors (and the areas that they cover) to maximise the detection, using a minimum number of sensors: '+str(numberofsensors)+'. '
      OptCover3=OptCover3+'<p style=\"text-align:center\"><img src=\''+file1[x][6:-19]+'map2.png\' height=\"400\"></p><p style=\"text-align:center\">Figure '+str(figure_count)+'. The best locations to place sensors (and the areas that they cover) to maximise detection from a <i>C<sub>seep</sub></i> pH change threshold of '+Thres[x][:5]+', using a minimum number of sensors: '+str(numberofsensors)+'.</p>'
    

  OUTPUTS=OUTPUTS+OptCover+OptCover2+OptCover3

  OUTPUT2=OUTPUT2+"<h2>Deployment strategies - Optimal Cover</h2>"+OptCover2+OptCover3

with open('storage/summary.html', 'r') as f:
  summary = f.read()  
    
OUTPUTS = OUTPUTS+summary

# Print output

template = env.get_template('storage/template-wt.html')
a = template.render(OUTPUTS=OUTPUTS, RISKMAP=RISKMAP, HYDRODYNAMICS=HYDRODYNAMICS, REF='{{ REF }}')

with open('storage/template-all.html', 'w') as f:
  f.write(a)

# -- References -- #

with open('storage/Scenarios_REF.html', 'r') as f:
  ref=f.readlines()
ref = [x.strip() for x in ref] 
REF=REF+ref

if os.path.isdir(filelocations+'geo'):
  with open('storage/GEO_REF.html', 'r') as f:
    ref=f.readlines()
  ref = [x.strip() for x in ref] 
  REF=REF+ref

if os.path.isdir(filelocations+'Advdiff'):
  with open('storage/ADVDIFF_REF.html', 'r') as f:
    ref=f.readlines()
  ref = [x.strip() for x in ref] 
  REF=REF+ref

if os.path.isdir(filelocations+'ROC'):
  with open('storage/ROC_REF.html', 'r') as f:
    ref=f.readlines()
  ref = [x.strip() for x in ref] 
  REF=REF+ref
  
if os.path.isdir(filelocations+'carbon'):
  with open('storage/carbon_REF.html', 'r') as f:
    ref=f.readlines()
  ref = [x.strip() for x in ref] 
  REF=REF+ref
  
if os.path.isdir(filelocations+'impacts'):
  with open('storage/impacts_REF.html', 'r') as f:
    ref=f.readlines()
  ref = [x.strip() for x in ref] 
  REF=REF+ref

if os.path.isdir(filelocations+'OptCover'):
  with open('storage/OptCover_REF.html', 'r') as f:
    ref=f.readlines()
  ref = [x.strip() for x in ref] 
  REF=REF+ref
  
REF = sorted(REF)
REF = list(dict.fromkeys(REF))
REF[0]=' '.join(REF)
template = env.get_template('storage/template-all.html')
a = template.render(REF=REF[0])
with open(filelocations+'Report.html', 'w') as f:
  f.write(a)
  
template = env.get_template('storage/template-ts.html')
a = template.render(Documentation=': Technical Summary',OUTPUTS=OUTPUT2)
with open(filelocations+'Technical-Summary.html', 'w') as f:
  f.write(a)

if sys.argv[1] == "lin":
  os.system('cp -r storage/Re-Run-lin.sh '+filelocations+'Re-Run.sh')
elif sys.argv[1] == "mac":
  os.system('cp -r storage/Re-Run-mac.sh '+filelocations+'Re-Run.sh')
elif sys.argv[1] == "win":
  os.system('cp -r storage/Re-Run.bat '+filelocations+'Re-Run.bat')
