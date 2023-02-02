import jinja2
import configparser
import os
import tkinter as tk
from tkinter import ttk
from PIL import ImageTk,Image

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
print("\nRequesting inputs for re-run\n")
print("******************************************************\n")

window=tk.Tk()
window.title('ACTOM Toolbox')
window.geometry("1024x800")
window.configure(bg='white')

frame = ScrollableFrame(window)

canvas = tk.Canvas(frame.scrollable_frame, width = 984, height = 170, bg='white', bd=0, highlightthickness=0)  
canvas.pack()  
img = ImageTk.PhotoImage(Image.open("storage/ACTOM.png"))  
canvas.create_image(0,0, anchor=tk.NW, image=img)

la = tk.Label(frame.scrollable_frame, text = "\nAct on Offshore Monitoring Toolbox Re-Run", background="white")
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

le = tk.Label(frame.scrollable_frame, text = "\nRe-Run:", background="white")
le.pack(side=tk.TOP, anchor="w", padx=12)
le.config(font=("Times New Roman (serif)", 18, 'bold'))

OUTPUTS='\nThe ACTOM Toolbox has been designed so that the Tracer Transport Model, C\u209B\u2091\u2091\u209A and the Rate of Change Anomaly Criteria tools all only require to be run once for a specific region which speeds up the process dramatically when anlysing different releases. Therefore to re-run the analysis with different leakage rates or detection thresholds, set these below. Submitting these details will re-run the required tools automatically.'

lf = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lf.config(font=("Times New Roman (serif)", 12))
lf.pack(anchor="w", padx=12)

config=configparser.ConfigParser(allow_no_value=True)
config.read('input/data.ini')
Thresgen=config['General']['threshold-ph']
rate = config['General']['rate']
rateno=rate
rate_units = config['General']['rate-units']
RateUnitsTime=rate_units.partition('/')
rateTime=RateUnitsTime[2]
rate = rate+' '+rate_units

OUTPUTS='\nThe leakage rate set for this report was given as: '+rate+'.\n'

lg = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lg.config(font=("Times New Roman (serif)", 12))
lg.pack(anchor="w", padx=12)

current_value = tk.DoubleVar()
current_value1 = tk.StringVar()
def get_current_value():
    return '{: .2f}'.format(current_value.get())
def slider_changed(event):
    value_label.configure(text=get_current_value())
    current_value1.set(str(get_current_value()))
    rateBox.configure(textvariable=current_value1)
frameone = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
slider_label = tk.Label(frameone, text = '• Enter a new rate (from each source) using the slider, or override the value below:\n', background="white", wraplength=984, justify="left")
slider_label.config(font=("Times New Roman (serif)", 12))
slider_label.pack(side=tk.LEFT, anchor="w", padx=12)
frameone.pack(anchor="w", padx=12)
frameonep5 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
slider = ttk.Scale(frameonep5,from_=0,to=1000,orient='horizontal',command=slider_changed, variable=current_value, length=960)
slider.pack(side=tk.LEFT, anchor="w", padx=12)
value_label = ttk.Label(frameonep5,text=get_current_value())
frameonep5.pack(anchor="w", padx=12)
frameonep6 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
rateBox = ttk.Entry(frameonep6,textvariable=current_value1, width = 8)
rateBox.config(font=("Times New Roman (serif)", 12))
rateBox.pack(side=tk.LEFT, anchor="w", padx=2)
masschoices = ['tonnes', 'kilograms', 'grams', 'moles']
massoptions = tk.StringVar(frameonep6)
massoptions.set('tonnes')
mass = tk.OptionMenu(frameonep6, massoptions, *masschoices)
mass.config(font=("Times New Roman (serif)", 12))
mass.pack(side=tk.LEFT, anchor="w", padx=2)
lg = tk.Label(frameonep6, text = '/', background="white", wraplength=984, justify="left")
lg.config(font=("Times New Roman (serif)", 12))
lg.pack(side=tk.LEFT, anchor="w", padx=0)
timechoices = ['year', 'month', 'week', 'day', 'hour', 'minute', 'second']
timeoptions = tk.StringVar(frameonep6)
timeoptions.set('year')
time = tk.OptionMenu(frameonep6, timeoptions, *timechoices)
time.config(font=("Times New Roman (serif)", 12))
time.pack(side=tk.LEFT, anchor="w", padx=2);
frameonep6.pack(expand=True)

OUTPUTS='\nThe user set pH detection threholds for this report were given as: '+Thresgen+'.'
li = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
li.config(font=("Times New Roman (serif)", 12))
li.pack(anchor="w", padx=12)

frametwo= tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
OUTPUTS='\n• Enter new set pH detection thresholds, seperated by a comma:\n'
lh = tk.Label(frametwo, text = OUTPUTS, background="white", wraplength=984, justify="left")
lh.config(font=("Times New Roman (serif)", 12))
lh.pack(side=tk.LEFT, anchor="w", padx=12)
Thr = tk.StringVar()
Thr.set(Thresgen)
ThrBox = tk.Entry(frametwo,textvariable = Thr, width = 20)
ThrBox.config(font=("Times New Roman (serif)", 12))
ThrBox.pack(side=tk.LEFT, anchor="w", padx=2)
frametwo.pack(anchor="w", padx=12)

OUTPUTS='A new directory will be created for the updated run based on the start date and time.\n'
lj = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lj.config(font=("Times New Roman (serif)", 12))
lj.pack(anchor="w", padx=12)

MyButton1 = tk.Button(frame.scrollable_frame, text="Re-Run", width=10, command=lambda: window.quit())
MyButton1.pack(anchor="c", padx=12)

lk = tk.Label(frame.scrollable_frame, text = "\nDisclaimer", background="white")
lk.pack(anchor="w", padx=12)
lk.config(font=("Times New Roman (serif)", 18, 'bold'))

OUTPUTS = 'All figures and data are indication based on simulations. The precision and quality will vary based on the quality, precision and range of data inputted to each tool.'

ll = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
ll.pack(anchor="w", padx=12)
ll.config(font=("Times New Roman (serif)", 12))

lm = tk.Label(frame.scrollable_frame, text = "\nAknowledgements", background="white")
lm.pack(anchor="w", padx=12)
lm.config(font=("Times New Roman (serif)", 18, 'bold'))

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

print("Writing settings for re-run...")

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

config.set('General', 'rate', rateBox.get())
config.set('General', 'rate-units', rateun)
config.set('General', 'threshold-ph', ThrBox.get())

with open('input/data.ini', 'w') as configfile:
    config.write(configfile)

OUTPUTS='\nThe leakage rate of '+rate+' has been replaced by: '+rateBox.get()+' '+rateun+'.\nThe user set pH detection threholds of '+Thresgen+' have been replaced by: '+ThrBox.get()+'.'

print(OUTPUTS)

config.read('input/data.ini')
genopt=config['General']['options']

if os.path.isfile('input/logs/Run-All.sh'):
  with open('input/logs/Run-All.sh', 'r') as runallfile: first_line = runallfile.readlines()[0:1]
  first_line=''.join(first_line)
  with open('storage/Re-Run-Toolbox.sh', 'r') as original1: data = original1.read()
  if genopt=='':
    data=data.replace("$options","")
  else:
    data=data.replace("$options","-a STDERR")
  with open('input/Re-Run-Toolbox.sh', 'w') as modified1: modified1.write(first_line+data)
if os.path.isfile('input/logs/Run-All.bat'):
  with open('storage/Re-Run-Toolbox.bat', 'r') as original2: data2 = original2.read()
  if genopt=='':
    data2=data2.replace("$options","")
  else:
    data2=data2.replace("$options","-a STDERR")
  with open('input/Re-Run-Toolbox.bat', 'w') as modified2: modified2.write(data2)

print('The re-run files have been generated')

