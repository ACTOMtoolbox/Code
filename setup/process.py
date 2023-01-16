import configparser
import os
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from PIL import ImageTk,Image

#def select_dir1():
#    global Cseepdir
#    Cseepdir=fd.askdirectory(initialdir='/',
#    title="Select Directory")
    
#def select_dir2():
#    global ROCdir
#    ROCdir=fd.askdirectory(initialdir='/',
#    title="Select Directory")

class ScrollableFrame(tk.Frame):
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)
        
        canvas = tk.Canvas(self, width = 1008, height = 800, background="#ffffff")
        self.viewPort = tk.Frame(canvas, background="#ffffff")   
        scrollbar = tk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.scrollable_frame = tk.Frame(canvas,background="#ffffff")

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")         

print("\n******************************************************")
print("\nRequesting inputs for ACTOM Toolbox Setup\n")
print("******************************************************\n")

Cseepdir=''
ROCdir=''

window=tk.Tk()
window.title('ACTOM Toolbox')
window.geometry("1024x800")
window.configure(bg='white')

frame = ScrollableFrame(window)

canvas = tk.Canvas(frame.scrollable_frame, width = 984, height = 170, bg='white', bd=0, highlightthickness=0)  
canvas.pack()  
img = ImageTk.PhotoImage(Image.open("storage/ACTOM.png"))  
canvas.create_image(0,0, anchor=tk.NW, image=img)

la = tk.Label(frame.scrollable_frame, text = "\nAct on Offshore Monitoring Toolbox Setup", background="white")
la.pack(side=tk.TOP, anchor="w", padx=12)
la.config(font=("Times New Roman (serif)", 24, 'bold'))

OUTPUTS='\nThe research in ACTOM works for the advancement of offshore monitoring to ensure alignment of CO\u2082 storage projects with national and international regulations and societal concerns. The ACTOM toolbox (shown in Figure 1.) equips operators with the ability to plan strategies under site-specific conditions and provides regulators a reliable and independent assessment of proposed monitoring strategies from license applicants. In addition to the tool assisting in the technical design of monitoring programs, it can be used to enhance communication with governments and the public in view of Marine Spatial Planning and Responsible Research and Innovation.\n\nThe toolbox is designed to provide value over a range of field cases with diverse subsurface geology and environmental marine characteristics. The end-product of the toolbox is to aid users in defining a monitoring plan that will satisfy local stakeholders.\n\nThe ACTOM toolbox contains algorithms for determining how the main monitoring aims required by guidance and regulations can be achieved at each site. Namely the algorithm optimizes approaches for detecting anomalies that could signal leakage in the marine environment, and then determining if these anomalies represent leakage. This is not an easy task because CO\u2082 is already an integral and dynamic part of the marine environment. Which is already undergoing shifts due to climate change. Through using geophysical data on shallow features such as chimneys or fractures, the likely areas where leakage would emanate can be predicted. Knowing the hydrodynamics and bio-chemical processes in marine sediments and seawater allow strategies for location and attribution, setting trigger points for when more action is needed and what instruments can be used for measurements, can be optimized for finding and attributing leakage emissions accounting and environmental protection.\n\n'

lb = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lb.config(font=("Times New Roman (serif)", 12,))
lb.pack(anchor="w", padx=12)

canvas2 = tk.Canvas(frame.scrollable_frame, bg='white', width = 800, height = 435, bd=0, highlightthickness=0)  
img2 = ImageTk.PhotoImage(Image.open("storage/toolbox.png"))  
canvas2.create_image(400,217,image=img2,anchor=tk.CENTER)
canvas2.pack(anchor='center')

lb2 = tk.Label(frame.scrollable_frame, text = 'Figure 1. Schematic of the toolbox, with tools interacting along with the site specific and user programmable information.', background="white", wraplength=984, justify="center")
lb2.config(font=("Times New Roman (serif)", 12,))
lb2.pack(anchor='center', padx=12)

lc = tk.Label(frame.scrollable_frame, text = "\nThe ACTOM toolbox:", background="white")
lc.pack(side=tk.TOP, anchor="w", padx=12)
lc.config(font=("Times New Roman (serif)", 18, 'bold'))

OUTPUTS='\n • Enables regulators to quantifiably assess that a proposed monitoring strategy delivers an acceptable standard of \n    assurance.\n • Enables operators to properly plan, cost and adapt monitoring strategies to site specific circumstances.\n • Enables regulators and operators to communicate to the effectiveness of proposed monitoring strategies to enable\n    informed societal consensus in view of marine spatial planning.'

ld = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
ld.config(font=("Times New Roman (serif)", 12))
ld.pack(anchor="w", padx=12)

le = tk.Label(frame.scrollable_frame, text = "\nSetup:", background="white")
le.pack(side=tk.TOP, anchor="w", padx=12)
le.config(font=("Times New Roman (serif)", 18, 'bold'))

OUTPUTS='\nThe ACTOM Toolbox has been designed so that with a few settings and input files it can run for a specifc region.'

lf = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lf.config(font=("Times New Roman (serif)", 12))
lf.pack(anchor="w", padx=12)

le = tk.Label(frame.scrollable_frame, text = "\nGeneral:", background="white")
le.pack(side=tk.TOP, anchor="w", padx=12)
le.config(font=("Times New Roman (serif)", 14, 'bold'))

OUTPUTS='\nFor general settings used by many of the tools, we first choose if we want to use debugging mode. This provides on screen information and warnings as the toolbox is running.\n\n The leakage rate from each source and any detection thresholds are also required to be set (beyond those predicted within the toolbox itself).\n'

lf = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lf.config(font=("Times New Roman (serif)", 12))
lf.pack(anchor="w", padx=12)

c4='Yes'
framefour3 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
c4_v1=tk.StringVar()
c4 = tk.Checkbutton(framefour3, text='Do we want to use debugging mode?', variable=c4_v1,
	onvalue='Yes',offvalue='No',bg="white",highlightthickness=0,bd=0)
c4.select()
c4.config(font=("Times New Roman (serif)", 12))
c4.pack(side=tk.LEFT)
framefour3.pack(anchor="w", padx=12)

frameone = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
OUTPUTS='\nEnter a chosen leakage rate:\n'
lg = tk.Label(frameone, text = OUTPUTS, background="white", wraplength=984, justify="left")
lg.config(font=("Times New Roman (serif)", 12))
lg.pack(side=tk.LEFT, anchor="w", padx=12)
T = tk.StringVar()
T.set('')
rateBox = tk.Entry(frameone,textvariable = T, width = 8)
rateBox.config(font=("Times New Roman (serif)", 12))
rateBox.pack(side=tk.LEFT, anchor="w", padx=12)
masschoices = ['tonnes', 'kilograms', 'grams', 'moles']
massoptions = tk.StringVar(frameone)
massoptions.set('kilograms')
mass = tk.OptionMenu(frameone, massoptions, *masschoices)
mass.config(font=("Times New Roman (serif)", 12))
mass.pack(side=tk.LEFT, anchor="w", padx=2)
lg = tk.Label(frameone, text = '/', background="white", wraplength=984, justify="left")
lg.config(font=("Times New Roman (serif)", 12))
lg.pack(side=tk.LEFT, anchor="w", padx=0)
timechoices = ['year', 'month', 'week', 'day', 'hour', 'minute', 'second']
timeoptions = tk.StringVar(frameone)
timeoptions.set('second')
time = tk.OptionMenu(frameone, timeoptions, *timechoices)
time.config(font=("Times New Roman (serif)", 12))
time.pack(side=tk.LEFT, anchor="w", padx=2);
frameone.pack(anchor="w", padx=12)

frametwo= tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
OUTPUTS='\nEnter any set pH detection thresholds, seperated by a comma:\n'
lh = tk.Label(frametwo, text = OUTPUTS, background="white", wraplength=984, justify="left")
lh.config(font=("Times New Roman (serif)", 12))
lh.pack(side=tk.LEFT, anchor="w", padx=12)
Thr = tk.StringVar()
Thr.set('')
ThrBox = tk.Entry(frametwo,textvariable = Thr, width = 20)
ThrBox.config(font=("Times New Roman (serif)", 12))
ThrBox.pack(side=tk.LEFT, anchor="w", padx=2)
frametwo.pack(anchor="w", padx=12)

lea = tk.Label(frame.scrollable_frame, text = "C\u209B\u2091\u2091\u209A:", background="white")
lea.pack(side=tk.TOP, anchor="w", padx=12)
lea.config(font=("Times New Roman (serif)", 14, 'bold'))

frametwofive = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
OUTPUTS='\nEnter the directory of input files (if not using built in data sets):\n'
lia = tk.Label(frametwofive, text = OUTPUTS, background="white", wraplength=984, justify="left")
lia.config(font=("Times New Roman (serif)", 12))
lia.pack(side=tk.LEFT, anchor="w", padx=12)
Thrdir = tk.StringVar()
Thrdir.set('')
ThrdirBox = tk.Entry(frametwofive,textvariable = Thrdir, width = 20)
ThrdirBox.config(font=("Times New Roman (serif)", 12))
ThrdirBox.pack(side=tk.LEFT, anchor="w", padx=2)
frametwofive.pack(anchor="w", padx=12)

le = tk.Label(frame.scrollable_frame, text = "Rate of Change Anomaly Criteria:", background="white")
le.pack(side=tk.TOP, anchor="w", padx=12)
le.config(font=("Times New Roman (serif)", 14, 'bold'))

OUTPUTS='\nFor the Rate of Change Anomaly Criteria, the time between pH measurements in a single location is required.'
lf = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lf.config(font=("Times New Roman (serif)", 12))
lf.pack(anchor="w", padx=12)

framethree = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
OUTPUTS='\nEnter a sensor sampling frequency:\n'
li = tk.Label(framethree, text = OUTPUTS, background="white", wraplength=984, justify="left")
li.config(font=("Times New Roman (serif)", 12))
li.pack(side=tk.LEFT, anchor="w", padx=12)
T = tk.StringVar()
T.set('')
sampfreq = tk.Entry(framethree,textvariable = T, width = 8)
sampfreq.config(font=("Times New Roman (serif)", 12))
sampfreq.pack(side=tk.LEFT, anchor="w", padx=2)
sampfreqtimechoices = ['hours', 'minutes', 'seconds']
sampfreqtimeoptions = tk.StringVar(framethree)
sampfreqtimeoptions.set('seconds')
sampfreqtime = tk.OptionMenu(framethree, sampfreqtimeoptions, *sampfreqtimechoices)
sampfreqtime.config(font=("Times New Roman (serif)", 12))
sampfreqtime.pack(side=tk.LEFT, anchor="w", padx=2);
framethree.pack(anchor="w", padx=12)

framethreepfive = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
OUTPUTS='\nSelect the directory of input files (if not using built in data sets):\n'
lia = tk.Label(framethreepfive, text = OUTPUTS, background="white", wraplength=984, justify="left")
lia.config(font=("Times New Roman (serif)", 12))
lia.pack(side=tk.LEFT, anchor="w", padx=12)
Thrdir2 = tk.StringVar()
Thrdir2.set('')
ThrdirBox2 = tk.Entry(framethreepfive,textvariable = Thrdir2, width = 20)
ThrdirBox2.config(font=("Times New Roman (serif)", 12))
ThrdirBox2.pack(side=tk.LEFT, anchor="w", padx=2)
framethreepfive.pack(anchor="w", padx=12)


le = tk.Label(frame.scrollable_frame, text = "Optimal Cover Tool:", background="white")
le.pack(side=tk.TOP, anchor="w", padx=12)
le.config(font=("Times New Roman (serif)", 14, 'bold'))

OUTPUTS='\nFor the Optimal Cover tool, we need to know if we are looking at all points in a leakage cluster, or just a single point at the leakage source. Are we looking at the maximum total number of leaks detected, or detect those the leaks with higher probability proportionally many times, and the maximum computational time we allow to compute the data.\n'
lj = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lj.config(font=("Times New Roman (serif)", 12))
lj.pack(anchor="w", padx=12)

c1='False'
framefour = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
c1_v1=tk.StringVar()
c1 = tk.Checkbutton(framefour, text='Include leakage clusters?', variable=c1_v1,
	onvalue='True',offvalue='False',bg="white",highlightthickness=0,bd=0)
c1.config(font=("Times New Roman (serif)", 12))
c1.deselect()
c1.pack(side=tk.LEFT)
framefour.pack(anchor="w", padx=12)

li1 = tk.Label(frame.scrollable_frame, text = '', background="white", wraplength=984, justify="left")
li1.config(font=("Times New Roman (serif)", 12))
li1.pack()

c2='0'
framefour1 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
c2_v1=tk.StringVar()
c2 = tk.Checkbutton(framefour1, text='Detect leaks with higher probability proportionally many more times?', variable=c2_v1,
	onvalue='1',offvalue='0',bg="white",highlightthickness=0,bd=0)
c2.config(font=("Times New Roman (serif)", 12))
c2.deselect()
c2.pack(side=tk.LEFT)
framefour1.pack(anchor="w", padx=12)

#framefive = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
#OUTPUTS='\nMaximum computational time we allow to compute the data (seconds):\n'
#lo = tk.Label(framefive, text = OUTPUTS, background="white", wraplength=984, justify="left")
#lo.config(font=("Times New Roman (serif)", 12))
#lo.pack(side=tk.LEFT, anchor="w", padx=12)
#T = tk.StringVar()
#T.set('')
#comptime = tk.Entry(framefive,textvariable = T, width = 8)
#comptime.config(font=("Times New Roman (serif)", 12))
#comptime.pack(side=tk.LEFT, anchor="w", padx=12)
#framefive.pack(anchor="w", padx=12)
comptime=120

le = tk.Label(frame.scrollable_frame, text = "Impact Analysis:", background="white")
le.pack(side=tk.TOP, anchor="w", padx=12)
le.config(font=("Times New Roman (serif)", 14, 'bold'))

OUTPUTS='\nFor the impact analysis, we need to know if we want to analyse local or global data, for example local data is usually higher resolution data surrounding each leakage source, but this slows down the toolbox, so is only for detailed analysis and testing. In the same aspect, are animations required for the global findings, all findings or none? Only the global results are analysed in the report.\n'
lja = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lja.config(font=("Times New Roman (serif)", 12))
lja.pack(anchor="w", padx=12)

c3='Global'
framefour2 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
c3_v1=tk.StringVar()
c3 = tk.Checkbutton(framefour2, text='Do we want to compute localised data not analysed in the report?', variable=c3_v1,
	onvalue='All',offvalue='Global',bg="white",highlightthickness=0,bd=0)
c3.deselect()
c3.config(font=("Times New Roman (serif)", 12))
c3.pack(side=tk.LEFT)
framefour2.pack(anchor="w", padx=12)

frameseven = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
OUTPUTS='\nDo we also want animations?\n'
lp = tk.Label(frameseven, text = OUTPUTS, background="white", wraplength=984, justify="left")
lp.config(font=("Times New Roman (serif)", 12))
lp.pack(side=tk.LEFT, anchor="w", padx=12)
anichoices = ['Yes all', 'Global only', 'No']
anioptions = tk.StringVar(frameseven)
anioptions.set('Global only')
ani = tk.OptionMenu(frameseven, anioptions, *anichoices)
ani.config(font=("Times New Roman (serif)", 12))
ani.pack(side=tk.LEFT, anchor="w", padx=12);
frameseven.pack(anchor="w", padx=12)


OUTPUTS='\nA new directory will be created for the toolbox run data based on the start date and time.\n'
lm = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lm.config(font=("Times New Roman (serif)", 12))
lm.pack(anchor="w", padx=12)
MyButton1 = tk.Button(frame.scrollable_frame, text="Run", width=10, command=lambda: window.quit())
MyButton1.pack(anchor="c", padx=12)

lk = tk.Label(frame.scrollable_frame, text = "\nDisclaimer", background="white")
lk.pack(anchor="w", padx=12)
lk.config(font=("Times New Roman (serif)", 18))

OUTPUTS = 'All figures and data are indication based on simulations. The precision and quality will vary based on the quality, precision and range of data inputted to each tool.'

ll = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
ll.pack(anchor="w", padx=12)
ll.config(font=("Times New Roman (serif)", 12))

lm = tk.Label(frame.scrollable_frame, text = "\nAknowledgements", background="white")
lm.pack(anchor="w", padx=12)
lm.config(font=("Times New Roman (serif)", 18))

OUTPUTS='\nThis project, ACTOM, is funded through the ACT programme (Accelerating CCS Technologies, Horizon2020 Project No 294766). Financial contributions made from; The Research Council of Norway, (RCN), Norway, Netherlands Enterprise Agency (RVO), Netherlands, Department for Business, Energy & Industrial Strategy (BEIS) together with extra funding from NERC and EPSRC research councils, United Kingdom, US-Department of Energy (US-DOE), USA. In-kind contributions from the University of Bergen are gratefully acknowledged.\n'

ln = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
ln.pack(anchor="w", padx=12)
ln.config(font=("Times New Roman (serif)", 12))

canvas3 = tk.Canvas(frame.scrollable_frame, bg='white', width = 960, height = 51, bd=0, highlightthickness=0)  
img3 = ImageTk.PhotoImage(Image.open("storage/partners.png"))  
canvas3.create_image(480,26,image=img3,anchor=tk.CENTER)
canvas3.pack(anchor='center')

lm = tk.Label(frame.scrollable_frame, text = '\n', background="white", wraplength=984, justify="left")
lm.pack(anchor="w", padx=12)
lm.config(font=("Times New Roman (serif)", 12))

window.resizable(False, False)
frame.pack()
window.mainloop()

print("Writing settings for toolbox run...\n")

rateun=''

if massoptions.get() == "kilograms":
   rateun='kg'
if massoptions.get() == "grams":
   rateun='g'
if massoptions.get() == "tonnes":
   rateun='t'
if massoptions.get() == "moles":
   rateun='m'
   
rateun=rateun+'/'

if timeoptions.get() == "second":
   rateun=rateun+'s'
if timeoptions.get() == "minute":
   rateun=rateun+'min'
if timeoptions.get() == "hour":
   rateun=rateun+'hr'
if timeoptions.get() == "day":
   rateun=rateun+'day'
if timeoptions.get() == "week":
   rateun=rateun+'week'
if timeoptions.get() == "month":
   rateun=rateun+'mon'
if timeoptions.get() == "year":
   rateun=rateun+'year'

config = configparser.ConfigParser()
config.add_section('General')
config.set('General', 'rate', rateBox.get())
config.set('General', 'rate-units', rateun)
config.set('General', 'threshold-ph', ThrBox.get())
config.add_section('RateOfChange')
config.set('RateOfChange', 'Sampling_Frequency', sampfreq.get())
config.set('RateOfChange', 'Sampling_Frequency_Units', sampfreqtimeoptions.get())
config.add_section('OptCover')
config.set('OptCover', 'include_cluster_points', c1_v1.get())
config.set('OptCover', 'cost_function', c2_v1.get())
config.set('OptCover', 'time_limit', str(comptime))
config.add_section('Impacts')
config.set('Impacts', 'animation', anioptions.get())
config.set('Impacts', 'local', c3_v1.get())
config.set('Impacts', 'output', 'pH')

with open('input/data.ini', 'w') as configfile:
    config.write(configfile)

if c4_v1.get()=='Yes':
   OUTPUTS='Debugging mode is on.\n\n'
else:
   OUTPUTS='Debugging mode is off.\n\n'

OUTPUTS='The leakage rate of has been set as: '+rateBox.get()+' '+rateun+', and the user set pH detection thresholds have been set as: '+ThrBox.get()+'.'

if ThrdirBox.get() !='':
   OUTPUTS=OUTPUTS+'\n\nThe directory for C\u209B\u2091\u2091\u209A input data has been set as '+ThrdirBox.get()+'.'

OUTPUTS=OUTPUTS+'\n\nThe Rate of Change Anomaly Criteria sampling frequency has been set as: '+sampfreq.get()+' '+sampfreqtimeoptions.get()

if ThrdirBox2.get() !='':
   OUTPUTS=OUTPUTS+', with the directory for input data set as '+ThrdirBox2.get()

OUTPUTS=OUTPUTS+'.\n\nThe Optimal Cover tool has been set to '

if c1_v1.get()=='True':
   OUTPUTS=OUTPUTS+'accept '
else:
   OUTPUTS=OUTPUTS+'ignore '
   
OUTPUTS=OUTPUTS+'cluster points, and '

if c2_v1.get()=='1':
   OUTPUTS=OUTPUTS+'Detect leaks with higher probability proportionally many more times.'
else:
   OUTPUTS=OUTPUTS+'analyse the maximum total number of leaks detected.'

OUTPUTS=OUTPUTS+'\n\nThe impacts will analyse'

if c3_v1.get()=='All':
   OUTPUTS=OUTPUTS+' all the outputs,'
elif c3_v1.get()=='Global':
   OUTPUTS=OUTPUTS+' only the global outputs,'

OUTPUTS=OUTPUTS+' and '
if anioptions.get()=='All':
   OUTPUTS=OUTPUTS+'show all'
elif anioptions.get()=='Global only':
   OUTPUTS=OUTPUTS+'show only global'
elif anioptions.get()=='No':
   OUTPUTS=OUTPUTS+'show no'
   
OUTPUTS=OUTPUTS+' animations.'

print(OUTPUTS)

with open('storage/Run-All.sh', 'r') as original: data = original.read()


if ThrdirBox.get() !='':
    data = data.replace("#insertcseepinputhere", "--mount type=bind,source="+ThrdirBox.get()+",target=/srv/actom-app/input/external \\")
else:
    data = data.replace("\n          #insertcseepinputhere", "")
if ThrdirBox2.get() !='':
    data = data.replace("#insertrocinputhere", "--mount type=bind,source="+ThrdirBox2.get()+",target=/srv/actom-app/input/external \\")
else:
    data = data.replace("\n          #insertrocinputhere", "")
    
#print(data)

if c4_v1.get()=='Yes':
   first_line="#!/bin/bash\noptions=\'\'\n"
else:
   first_line="#!/bin/bash\noptions=\'-a STDERR\'\n"
with open('input/Run-All.sh', 'w') as modified: modified.write(first_line+data)

print('\nThe run files have been generated...\n\n')

