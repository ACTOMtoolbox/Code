import configparser
import os
import re
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from PIL import ImageTk,Image

os.system('cp storage/Setup.ini input/Setup.ini')

off_color = "red"
on_color = "black"

def TTMframe():
    if ttmbut_v1.get() == "True":ttmbut["fg"] = on_color     # if (get current checkbutton state) is "1" then....
    else:ttmbut["fg"] = off_color
    if ttmbut_v1.get() == "True":
       OUTPUTS='\nDue to the large number of settings within the Tracer Transport Model, a manual to adjust the settings is available from the following address:\n\nhttps://htmlpreview.github.io/?https://github.com/ACTOMtoolbox/Code/blob/main/advdiff/Doc/docs/_build/html/about.html\n\nPlease adjust all settings as required BEFORE clicking on \'Run\' below.\n\nThe Setup.ini file has been placed in your working directory that has been generated in the folder this program has run from, titled: \'RUN-\' and then date and time.\n'
       lia21 = tk.Label(framefour51, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lia21.config(font=("Times New Roman (serif)", 12))
       lia21.pack(anchor="w", padx=12)
       ca1='False'
       
       frametwofive22 = tk.Frame(framefour51, bg='white',width = 1008)
       OUTPUTS='\n• Enter the directory of input files (if not using built in data sets):\n'
       lia22 = tk.Label(frametwofive22, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lia22.config(font=("Times New Roman (serif)", 12))
       lia22.pack(side=tk.LEFT, anchor="w", padx=12)
       Thrdir22 = tk.StringVar()
       Thrdir22.set('')
       Thrdir22Box = tk.Entry(frametwofive22,textvariable = Thrdir22, width = 20)
       Thrdir22Box.config(font=("Times New Roman (serif)", 12))
       Thrdir22Box.pack(side=tk.LEFT, anchor="w", padx=2)
       frametwofive22.pack(anchor="w", padx=12)
    if ttmbut_v1.get()=='False':
       for widgets in framefour51.winfo_children():
           widgets.pack_forget()
       framefour51.config(height=0)
       
              
def cseepframe():
    if cseepbut_v1.get() == "True":cseepbut["fg"] = on_color     # if (get current checkbutton state) is "1" then....
    else:cseepbut["fg"] = off_color
    if cseepbut_v1.get()=='True':
       frametwofive = tk.Frame(frametwofivea, bg='white',width = 1008)
       OUTPUTS='\n• Enter the directory of input files (if not using built in data sets):\n'
       lia = tk.Label(frametwofive, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lia.config(font=("Times New Roman (serif)", 12))
       lia.pack(side=tk.LEFT, anchor="w", padx=12)
       Thrdir = tk.StringVar()
       Thrdir.set('')
       ThrdirBox = tk.Entry(frametwofive,textvariable = Thrdir, width = 20)
       ThrdirBox.config(font=("Times New Roman (serif)", 12))
       ThrdirBox.pack(side=tk.LEFT, anchor="w", padx=2)
       frametwofive.pack(anchor="w", padx=12)
       globals().update(locals())
    if cseepbut_v1.get()=='False':
       for widgets in frametwofivea.winfo_children():
           widgets.pack_forget()
       frametwofivea.config(height=0)
       
def rocframe():
    if rocbut_v1.get() == "True":rocbut["fg"] = on_color     # if (get current checkbutton state) is "1" then....
    else:rocbut["fg"] = off_color
    if rocbut_v1.get()=='True':       
       OUTPUTS='\nFor the Rate of Change Anomaly Criteria, the time between pH measurements in a single location is required.'
       lf = tk.Label(frametwofiveb, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lf.config(font=("Times New Roman (serif)", 12))
       lf.pack(anchor="w", padx=12)

       framethree = tk.Frame(frametwofiveb, bg='white',width = 1008)
       OUTPUTS='\n• Enter a sensor sampling frequency:\n'
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
       framethreepfive = tk.Frame(frametwofiveb, bg='white',width = 1008)
       OUTPUTS='\n• Enter the directory of input files (if not using built in data sets):\n'
       lia = tk.Label(framethreepfive, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lia.config(font=("Times New Roman (serif)", 12))
       lia.pack(side=tk.LEFT, anchor="w", padx=12)
       Thrdir2 = tk.StringVar()
       Thrdir2.set('')
       ThrdirBox2 = tk.Entry(framethreepfive,textvariable = Thrdir2, width = 20)
       ThrdirBox2.config(font=("Times New Roman (serif)", 12))
       ThrdirBox2.pack(side=tk.LEFT, anchor="w", padx=2)
       framethreepfive.pack(anchor="w", padx=12)
       globals().update(locals())
    if rocbut_v1.get()=='False':
       for widgets in frametwofiveb.winfo_children():
           widgets.pack_forget()
       frametwofiveb.config(height=0)

def optframe():
    if optbut_v1.get() == "True":optbut["fg"] = on_color     # if (get current checkbutton state) is "1" then....
    else:optbut["fg"] = off_color
    if optbut_v1.get()=='True':
       OUTPUTS='\nFor the Optimal Cover tool, we need to know if we are looking at all points in a leakage cluster, or just a single point at the leakage source. Are we looking at the maximum total number of leaks detected, or detect those the leaks with higher probability proportionally many times, and the maximum computational time we allow to compute the data.\n'
       lj = tk.Label(frametwofivec, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lj.config(font=("Times New Roman (serif)", 12))
       lj.pack(anchor="w", padx=12)
       c1='False'
       framefour = tk.Frame(frametwofivec, bg='white',width = 1008)
       c1_v1=tk.StringVar()
       c1 = tk.Checkbutton(framefour, text='Include leakage clusters?', variable=c1_v1,
	onvalue='True',offvalue='False',bg="white",highlightthickness=0,bd=0)
       c1.config(font=("Times New Roman (serif)", 12))
       c1.deselect()
       c1.pack(side=tk.LEFT)
       framefour.pack(anchor="w", padx=12)
	
       li1 = tk.Label(frametwofivec, text = '', background="white", wraplength=984, justify="left")
       li1.config(font=("Times New Roman (serif)", 12))
       li1.pack()
       c2='0'
       framefour1 = tk.Frame(frametwofivec, bg='white',width = 1008)
       c2_v1=tk.StringVar()
       c2 = tk.Checkbutton(framefour1, text='Detect leaks with higher probability proportionally many more times?', variable=c2_v1,
	onvalue='1',offvalue='0',bg="white",highlightthickness=0,bd=0)
       c2.config(font=("Times New Roman (serif)", 12))
       c2.deselect()
       c2.pack(side=tk.LEFT)
       framefour1.pack(anchor="w", padx=12)

       OUTPUTS='\n• Enter the directory of input files (if not using the Tracer Transport Model):\n'
       lia2 = tk.Label(frametwofivec, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lia2.config(font=("Times New Roman (serif)", 12))
       lia2.pack(side=tk.LEFT, anchor="w", padx=12)
       Thrdir3 = tk.StringVar()
       Thrdir3.set('')
       ThrdirBox3 = tk.Entry(frametwofivec,textvariable = Thrdir3, width = 20)
       ThrdirBox3.config(font=("Times New Roman (serif)", 12))
       ThrdirBox3.pack(side=tk.LEFT, anchor="w", padx=2)
       
       li2 = tk.Label(frametwofivec, text = '', background="white", wraplength=984, justify="left")
       li2.config(font=("Times New Roman (serif)", 12))
       li2.pack()
       comptime=120
       globals().update(locals())
    if optbut_v1.get()=='False':
       for widgets in frametwofivec.winfo_children():
           widgets.pack_forget()
       frametwofivec.config(height=0)

def carbframe():
    if carbbut_v1.get() == "True":carbbut["fg"] = on_color     # if (get current checkbutton state) is "1" then....
    else:carbbut["fg"] = off_color
    if carbbut_v1.get()=='True':
       frametwofivenew = tk.Frame(frametwofived, bg='white',width = 1008)
       OUTPUTS='\n• Enter the directory of input files (if not using the Tracer Transport Model):\n'
       lia4 = tk.Label(frametwofivenew, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lia4.config(font=("Times New Roman (serif)", 12))
       lia4.pack(side=tk.LEFT, anchor="w", padx=12)
       Thrdir4 = tk.StringVar()
       Thrdir4.set('')
       ThrdirBox4 = tk.Entry(frametwofivenew,textvariable = Thrdir4, width = 20)
       ThrdirBox4.config(font=("Times New Roman (serif)", 12))
       ThrdirBox4.pack(side=tk.LEFT, anchor="w", padx=2)
       frametwofivenew.pack(anchor="w", padx=12)
       globals().update(locals())
    if carbbut_v1.get()=='False':
       for widgets in frametwofived.winfo_children():
           widgets.pack_forget()
       frametwofived.config(height=0)

def impframe():
    if impbut_v1.get() == "True":impbut["fg"] = on_color     # if (get current checkbutton state) is "1" then....
    else:impbut["fg"] = off_color
    if impbut_v1.get()=='True':
       OUTPUTS='\nFor the impact analysis, we need to know if we want to analyse local or global data, for example local data is usually higher resolution data surrounding each leakage source, but this slows down the toolbox, so is only for detailed analysis and testing. In the same aspect, are animations required for the global findings, all findings or none? Only the global results are analysed in the report.\n'
       lja = tk.Label(frametwofivee, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lja.config(font=("Times New Roman (serif)", 12))
       lja.pack(anchor="w", padx=12)
       c3='Global'
       framefour2 = tk.Frame(frametwofivee, bg='white',width = 1008)
       c3_v1=tk.StringVar()
       c3 = tk.Checkbutton(framefour2, text='Do we want to compute localised data not analysed in the report?', variable=c3_v1,
	onvalue='All',offvalue='Global',bg="white",highlightthickness=0,bd=0)
       c3.deselect()
       c3.config(font=("Times New Roman (serif)", 12))
       c3.pack(side=tk.LEFT)
       framefour2.pack(anchor="w", padx=12)

       frameseven = tk.Frame(frametwofivee, bg='white',width = 1008)
       OUTPUTS='\n• Do we also want animations?\n'
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
       
       frameseven1 = tk.Frame(frametwofivee, bg='white',width = 1008)
       OUTPUTS='\n• Enter the directory of input files (if not using the Carbonate System):\n'
       lia2 = tk.Label(frameseven1, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lia2.config(font=("Times New Roman (serif)", 12))
       lia2.pack(side=tk.LEFT, anchor="w", padx=12)
       Thrdir5 = tk.StringVar()
       Thrdir5.set('')
       ThrdirBox5 = tk.Entry(frameseven1,textvariable = Thrdir5, width = 20)
       ThrdirBox5.config(font=("Times New Roman (serif)", 12))
       ThrdirBox5.pack(side=tk.LEFT, anchor="w", padx=2)
       frameseven1.pack(anchor="w", padx=12)
       
       globals().update(locals())
    if impbut_v1.get()=='False':
       for widgets in frametwofivee.winfo_children():
           widgets.pack_forget()
       frametwofivee.config(height=0)

def repframe():
    if repbut_v1.get() == "True":repbut["fg"] = on_color     # if (get current checkbutton state) is "1" then....
    else:repbut["fg"] = off_color
    if repbut_v1.get()=='True':
       frametwofivenew1 = tk.Frame(frametwofivef, bg='white',width = 1008)
       OUTPUTS='\n• Enter the directory of input files (if only reporting old datasets):\n'
       lia6 = tk.Label(frametwofivenew1, text = OUTPUTS, background="white", wraplength=984, justify="left")
       lia6.config(font=("Times New Roman (serif)", 12))
       lia6.pack(side=tk.LEFT, anchor="w", padx=12)
       Thrdir6 = tk.StringVar()
       Thrdir6.set('')
       ThrdirBox6 = tk.Entry(frametwofivenew1,textvariable = Thrdir6, width = 20)
       ThrdirBox6.config(font=("Times New Roman (serif)", 12))
       ThrdirBox6.pack(side=tk.LEFT, anchor="w", padx=2)
       frametwofivenew1.pack(anchor="w", padx=12)
       globals().update(locals())
    if repbut_v1.get()=='False':
       for widgets in frametwofivef.winfo_children():
           widgets.pack_forget()
       frametwofivef.config(height=0)


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

OUTPUTS='\nThe ACTOM Toolbox has been designed so that with a few settings and input files it can run for a specifc region. Any settings surrounded by a box are optional in certain cases. However, please read the instructions carefully.'

lf = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lf.config(font=("Times New Roman (serif)", 12))
lf.pack(anchor="w", padx=12)

le = tk.Label(frame.scrollable_frame, text = "\nGeneral:", background="white")
le.pack(side=tk.TOP, anchor="w", padx=12)
le.config(font=("Times New Roman (serif)", 18, 'bold'))

OUTPUTS='\nFor general settings used by many of the tools.\n\nThe leakage rate from each source and any detection thresholds are required to be set (beyond those predicted within the toolbox itself). Optional background values may be added, especially if not using C\u209B\u2091\u2091\u209A to process biogeochemical data.\n'

lf = tk.Label(frame.scrollable_frame, text = OUTPUTS, background="white", wraplength=984, justify="left")
lf.config(font=("Times New Roman (serif)", 12))
lf.pack(anchor="w", padx=12)

current_value = tk.DoubleVar()
current_value1 = tk.StringVar()
def get_current_value():
    return '{: .2f}'.format(current_value.get())
def slider_changed(event):
    value_label.configure(text=get_current_value())
    current_value1.set(str(get_current_value()))
    rateBox.configure(textvariable=current_value1)

#need to add link to discliamer on values#---------------------------------

frameone = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
slider_label = tk.Label(frameone, text = '• Select a chosen leakage rate (from each source) using the slider, or override the value below:\n', background="white", wraplength=984, justify="left")
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

frametwo= tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
OUTPUTS='\n• Enter any set pH detection thresholds, seperated by a comma:\n'
lh = tk.Label(frametwo, text = OUTPUTS, background="white", wraplength=984, justify="left")
lh.config(font=("Times New Roman (serif)", 12))
lh.pack(side=tk.LEFT, anchor="w", padx=12)
Thr = tk.StringVar()
Thr.set('')
ThrBox = tk.Entry(frametwo,textvariable = Thr, width = 20)
ThrBox.config(font=("Times New Roman (serif)", 12))
ThrBox.pack(side=tk.LEFT, anchor="w", padx=2)
frametwo.pack(anchor="w", padx=12)

c4='Yes'

frametwoaa= tk.Frame(frame.scrollable_frame, bg='white',width = 1008, highlightbackground="grey", highlightthickness=1)
frametwoa= tk.Frame(frametwoaa, bg='white',width = 1008)
OUTPUTS='\n• Enter background values for TA, DIC and seawater density*:'
lha = tk.Label(frametwoa, text = OUTPUTS, background="white", wraplength=984, justify="left")
lha.config(font=("Times New Roman (serif)", 12))
lha.pack(side=tk.LEFT, anchor="w", padx=12)
frametwoa.pack(anchor="w")
frametwoc= tk.Frame(frametwoaa, bg='white',width = 1008)
OUTPUTS='     *optional and overrides background data values, but is required if C\u209B\u2091\u2091\u209A is not in use.\n'
lhb = tk.Label(frametwoc, text = OUTPUTS, background="white", wraplength=984, justify="left")
lhb.config(font=("Times New Roman (serif)", 8))
lhb.pack(side=tk.LEFT, anchor="w", padx=12)
frametwoc.pack(anchor="w")
frametwob= tk.Frame(frametwoaa, bg='white',width = 984)
ThrTA = tk.StringVar()
ThrTA.set('TA')
ThrTABox = tk.Entry(frametwob,textvariable = ThrTA, width = 8)
ThrTABox.config(font=("Times New Roman (serif)", 12))
ThrTABox.pack(side=tk.LEFT, anchor="w", padx=2)
lhb1 = tk.Label(frametwob, text = 'µmol/kg   ', background="white", wraplength=984, justify="left")
lhb1.config(font=("Times New Roman (serif)", 12))
lhb1.pack(side=tk.LEFT, anchor="w", padx=2)
ThrDIC = tk.StringVar()
ThrDIC.set('DIC')
ThrDICBox = tk.Entry(frametwob,textvariable = ThrDIC, width = 8)
ThrDICBox.config(font=("Times New Roman (serif)", 12))
ThrDICBox.pack(side=tk.LEFT, anchor="w", padx=2)
lhb2 = tk.Label(frametwob, text = 'µmol/kg   ', background="white", wraplength=984, justify="left")
lhb2.config(font=("Times New Roman (serif)", 12))
lhb2.pack(side=tk.LEFT, anchor="w", padx=2)
frametwob.pack(expand=True)
ThrDens = tk.StringVar()
ThrDens.set('Density')
ThrDensBox = tk.Entry(frametwob,textvariable = ThrDens, width = 8)
ThrDensBox.config(font=("Times New Roman (serif)", 12))
ThrDensBox.pack(side=tk.LEFT, anchor="w", padx=2)
lhb3 = tk.Label(frametwob, text = 'kg/m\u00B3   ', background="white", wraplength=984, justify="left")
lhb3.config(font=("Times New Roman (serif)", 12))
lhb3.pack(side=tk.LEFT, anchor="w", padx=2)
frametwob.pack(expand=True)
framefour4a = tk.Frame(frametwoaa, bg='white',width = 1008)
output='                                                                                                                                                                                              '
loh = tk.Label(framefour4a, text = output, background="white", wraplength=984, justify="left")
loh.config(font=("Times New Roman (serif)", 12))
loh.pack(side=tk.LEFT, anchor="w", padx=12)
framefour4a.pack()
frametwoaa.pack(anchor="w", padx=12)

le2 = tk.Label(frame.scrollable_frame, text = "\nPlease Select ACTOM Tools and Settings:", background="white")
le2.pack(side=tk.TOP, anchor="w", padx=12)
le2.config(font=("Times New Roman (serif)", 18, 'bold'))
li3 = tk.Label(frame.scrollable_frame, text = '', background="white", wraplength=984, justify="left")
li3.config(font=("Times New Roman (serif)", 12))
li3.pack()

framefour5 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
ttmbut='False'
ttmbut_v1=tk.StringVar()
ttmbut = tk.Checkbutton(framefour5, text='The Tracer Transport Model', variable=ttmbut_v1,
	onvalue='True',offvalue='False',bg="white",highlightthickness=0,bd=0,command=TTMframe,fg=off_color)
ttmbut.config(font=("Times New Roman (serif)", 14, 'bold'))
ttmbut.deselect()
ttmbut.pack(side=tk.LEFT)
framefour5.pack(anchor="w", padx=12)
framefour51 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
framefour51.pack(anchor="w", padx=12)

framefour6 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
cseepbut='False'
cseepbut_v1=tk.StringVar()
cseepbut = tk.Checkbutton(framefour6, text="C\u209B\u2091\u2091\u209A", variable=cseepbut_v1,
	onvalue='True',offvalue='False',bg="white",highlightthickness=0,bd=0,command=cseepframe,fg=off_color)
cseepbut.config(font=("Times New Roman (serif)", 14, 'bold'))
cseepbut.deselect()
cseepbut.pack(side=tk.LEFT)
framefour6.pack(anchor="w", padx=12)
frametwofivea = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
frametwofivea.pack(anchor="w", padx=12)

framefour7 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
rocbut='False'
rocbut_v1=tk.StringVar()
rocbut = tk.Checkbutton(framefour7, text="Rate of Change Anomaly Criteria", variable=rocbut_v1,
	onvalue='True',offvalue='False',bg="white",highlightthickness=0,bd=0,command=rocframe,fg=off_color)
rocbut.config(font=("Times New Roman (serif)", 14, 'bold'))
rocbut.deselect()
rocbut.pack(side=tk.LEFT)
framefour7.pack(anchor="w", padx=12)
frametwofiveb = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
frametwofiveb.pack(anchor="w", padx=12)

framefour8 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
optbut='False'
optbut_v1=tk.StringVar()
optbut = tk.Checkbutton(framefour8, text="Optimal Cover Tool", variable=optbut_v1,
	onvalue='True',offvalue='False',bg="white",highlightthickness=0,bd=0,command=optframe,fg=off_color)
optbut.config(font=("Times New Roman (serif)", 14, 'bold'))
optbut.deselect()
optbut.pack(side=tk.LEFT)
framefour8.pack(anchor="w", padx=12)
frametwofivec = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
frametwofivec.pack(anchor="w", padx=12)

framefour9 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
carbbut='False'
carbbut_v1=tk.StringVar()
carbbut = tk.Checkbutton(framefour9, text="Carbonate System", variable=carbbut_v1,
	onvalue='True',offvalue='False',bg="white",highlightthickness=0,bd=0,command=carbframe,fg=off_color)
carbbut.config(font=("Times New Roman (serif)", 14, 'bold'))
carbbut.deselect()
carbbut.pack(side=tk.LEFT)
framefour9.pack(anchor="w", padx=12)
frametwofived = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
frametwofived.pack(anchor="w", padx=12)

framefour10 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
impbut='False'
impbut_v1=tk.StringVar()
impbut = tk.Checkbutton(framefour10, text="Impact Analysis", variable=impbut_v1,
	onvalue='True',offvalue='False',bg="white",highlightthickness=0,bd=0,command=impframe,fg=off_color)
impbut.config(font=("Times New Roman (serif)", 14, 'bold'))
impbut.deselect()
impbut.pack(side=tk.LEFT)
framefour10.pack(anchor="w", padx=12)
frametwofivee = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
frametwofivee.pack(anchor="w", padx=12)

framefour11 = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
repbut='False'
repbut_v1=tk.StringVar()
repbut = tk.Checkbutton(framefour11, text="Reporting", variable=repbut_v1,
	onvalue='True',offvalue='False',bg="white",highlightthickness=0,bd=0,command=repframe,fg=off_color)
repbut.config(font=("Times New Roman (serif)", 14, 'bold'))
repbut.deselect()
repbut.pack(side=tk.LEFT)
framefour11.pack(anchor="w", padx=12)
frametwofivef = tk.Frame(frame.scrollable_frame, bg='white',width = 1008)
frametwofivef.pack(anchor="w", padx=12)

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


# debugging inputs

if rateBox.get() == '':
   print('Leakage rate is missing, please restart.')
   quit()
try:
    float(rateBox.get())
except ValueError:
   print('Leakage rate value ('+rateBox.get()+') is not a number, please restart.')
   quit()
   
if rocbut_v1.get()=='True':
   if sampfreq.get() == '':
      print('Rate of Change sensor sampling frequency is missing, please restart.')
      quit()
   try:
      float(sampfreq.get())
   except ValueError:
      print('Rate of Change sensor sampling frequency ('+sampfreq.get()+') is not a number, please re-run.')
      quit()



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

with open('storage/Run-All.sh', 'r') as original1: data = original1.read()
with open('storage/Run-All.bat', 'r') as original2: data2 = original2.read()

OUTPUTS='The leakage rate of has been set as: '+rateBox.get()+' '+rateun+', and the user set pH detection thresholds have been set as: '+ThrBox.get()+'.'

TAoutput=re.findall(r"[-+]?(?:\d*\.*\d+)",ThrTABox.get())
DICoutput=re.findall(r"[-+]?(?:\d*\.*\d+)",ThrDICBox.get())
Densoutput=re.findall(r"[-+]?(?:\d*\.*\d+)",ThrDensBox.get())

if len(TAoutput)==1:
  config.set('General', 'ta', TAoutput[0])
  OUTPUTS=OUTPUTS+'The Total Alkalinity has been set as '+TAoutput+'µmol/kg.'
if len(DICoutput)==1:
  config.set('General', 'dic', DICoutput[0])
  OUTPUTS=OUTPUTS+'The Background Dissolved Inorganic Carbon has been set as '+DICoutput+'µmol/kg.'
if len(Densoutput)==1:
  config.set('General', 'dens', Densoutput[0])
  OUTPUTS=OUTPUTS+'The Background Seawater Density has been set as '+Densoutput+'kg/m\u00B3.'

#if ttm
if ttmbut_v1.get() == "True":
  if ThrdirBox22.get() !='':
    OUTPUTS=OUTPUTS+'\n\nThe directory for Tracer Transport Model input data has been set as '+ThrdirBox22.get()+'.'
    data = data.replace("#insertttminputhere", "--mount type=bind,source="+ThrdirBox22.get()+",target=/app/Indata \\")
    data2 = data2.replace("#insertttminputhere", "--mount type=bind,source="+ThrdirBox22.get()+",target=/app/Indata")
  else:
    data = data.replace("\n          #insertttminputhere", "")
    data2 = data2.replace("#insertttminputhere", "")
  OUTPUTS=OUTPUTS+''
else:
  data=re.sub('#1[^>]+#2', '', data)
  data2=re.sub('#1[^>]+#2', '', data2)
  OUTPUTS=OUTPUTS+'\n\nWARNING: TRACER TRANSPORT MODEL IS TURNED OFF...'
  OUTPUTS=OUTPUTS+'\nTHIS MAY AFFECT HOW OTHER TOOLS RUN!'
  
#if cseep
if cseepbut_v1.get() == "True":
  if ThrdirBox.get() !='':
    OUTPUTS=OUTPUTS+'\n\nThe directory for C\u209B\u2091\u2091\u209A input data has been set as '+ThrdirBox.get()+'.'
    data = data.replace("#insertcseepinputhere", "--mount type=bind,source="+ThrdirBox.get()+",target=/srv/actom-app/input/external \\")
    data2 = data2.replace("#insertcseepinputhere", "--mount type=bind,source="+ThrdirBox.get()+",target=/srv/actom-app/input/external")
  else:
    data = data.replace("\n          #insertcseepinputhere", "")
    data2 = data2.replace("#insertcseepinputhere", "")
else:
  data=re.sub('#2[^>]+#3', '', data)
  data2=re.sub('#2[^>]+#3', '', data2)
  OUTPUTS=OUTPUTS+'\n\nWARNING: C\u209B\u2091\u2091\u209A IS TURNED OFF...'
  OUTPUTS=OUTPUTS+'\nTHIS MAY AFFECT HOW OTHER TOOLS RUN!'
  
#if roc
if rocbut_v1.get() == "True":
  OUTPUTS=OUTPUTS+'\n\nThe Rate of Change Anomaly Criteria sampling frequency has been set as: '+sampfreq.get()+' '+sampfreqtimeoptions.get()
  if ThrdirBox2.get() !='':
    OUTPUTS=OUTPUTS+', with the directory for input data set as '+ThrdirBox2.get()
    data = data.replace("#insertrocinputhere", "--mount type=bind,source="+ThrdirBox2.get()+",target=/srv/actom-app/input/external \\")
    data2 = data2.replace("#insertrocinputhere", "--mount type=bind,source="+ThrdirBox2.get()+",target=/srv/actom-app/input/external")
  else:
    data = data.replace("\n          #insertrocinputhere", "")
    data2 = data2.replace("#insertrocinputhere", "")
  OUTPUTS=OUTPUTS+'.'
  config.add_section('RateOfChange')
  config.set('RateOfChange', 'Sampling_Frequency', sampfreq.get())
  config.set('RateOfChange', 'Sampling_Frequency_Units', sampfreqtimeoptions.get())
else:
  data=re.sub('#3[^>]+#4', '', data)
  data2=re.sub('#3[^>]+#4', '', data2)
  OUTPUTS=OUTPUTS+'\n\nWARNING: RATE OF CHANGE ANOMALY CRITERIA IS TURNED OFF!'

#if opt
if optbut_v1.get() == "True":
  OUTPUTS=OUTPUTS+'\n\nThe Optimal Cover tool has been set to '
  if c1_v1.get()=='True':
    OUTPUTS=OUTPUTS+'accept '
  else:
    OUTPUTS=OUTPUTS+'ignore '
  OUTPUTS=OUTPUTS+'cluster points, and '
  if c2_v1.get()=='1':
    OUTPUTS=OUTPUTS+'detect leaks with higher probability proportionally many more times'
  else:
    OUTPUTS=OUTPUTS+'analyse the maximum total number of leaks detected'
  if ThrdirBox3.get() !='':
    OUTPUTS=OUTPUTS+', with the directory for input data set as '+ThrdirBox3.get()
    data = data.replace("#insertoptinputhere", "--mount type=bind,source="+ThrdirBox3.get()+",target=/app/Input \\")
    data2 = data2.replace("#insertoptinputhere", "--mount type=bind,source="+ThrdirBox3.get()+",target=/app/Input")
  else:
    data = data.replace("#insertoptinputhere", "--mount type=bind,source=\"$(pwd)\"/../Advdiff/output,target=/app/Input \\")
    data2 = data2.replace("#insertoptinputhere", "--mount type=bind,source=\"%cd%\"/../Advdiff/output,target=/app/Input")
  OUTPUTS=OUTPUTS+'.'
  config.add_section('OptCover')
  config.set('OptCover', 'include_cluster_points', c1_v1.get())
  config.set('OptCover', 'cost_function', c2_v1.get())
  config.set('OptCover', 'time_limit', str(comptime))
else:
  data=re.sub('#4[^>]+#5', '', data)
  data2=re.sub('#4[^>]+#5', '', data2)
  OUTPUTS=OUTPUTS+'\n\nWARNING: OPTIMAL COVER TOOL IS TURNED OFF!'
   
#if carb
if carbbut_v1.get() == "True":
  if ThrdirBox4.get() !='':
    OUTPUTS=OUTPUTS+'\n\nThe directory for the carbonate system input data has been set as '+ThrdirBox4.get()+'.'
    data = data.replace("#insertcarbinputhere", "--mount type=bind,source="+ThrdirBox4.get()+",target=/external/input \\")
    data2 = data2.replace("#insertcarbinputhere", "--mount type=bind,source="+ThrdirBox4.get()+",target=/external/input")
  else:
    data = data.replace("#insertcarbinputhere", "--mount type=bind,source=\"$(pwd)\"/../Advdiff/output,target=/external/input \\")
    data2 = data2.replace("#insertcarbinputhere", "--mount type=bind,source=\"%cd%\"/../Advdiff/output,target=/external/input")
else:
  data=re.sub('#5[^>]+#6', '', data)
  data2=re.sub('#5[^>]+#6', '', data2)
  OUTPUTS=OUTPUTS+'\n\nWARNING: CARBONATE SYSTEM IS TURNED OFF...'
  OUTPUTS=OUTPUTS+'\nTHIS MAY AFFECT HOW OTHER TOOLS RUN!'
     
#if Imp
if impbut_v1.get() == "True":
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
  OUTPUTS=OUTPUTS+' animations'
  if ThrdirBox5.get() !='':
    OUTPUTS=OUTPUTS+', with the directory for input data set as '+ThrdirBox5.get()
    data = data.replace("#insertimpinputhere", "--mount type=bind,source="+ThrdirBox5.get()+",target=/external/input \\")
    data2 = data2.replace("#insertimpinputhere", "--mount type=bind,source="+ThrdirBox5.get()+",target=/external/input")
  else:
    data = data.replace("#insertimpinputhere", "--mount type=bind,source=\"$(pwd)\"/../carbon/output,target=/external/input \\")
    data2 = data2.replace("#insertimpinputhere", "--mount type=bind,source=\"%cd%\"/../carbon/output,target=/external/input")
  OUTPUTS=OUTPUTS+'.'
  config.add_section('Impacts')
  config.set('Impacts', 'animation', anioptions.get())
  config.set('Impacts', 'local', c3_v1.get())
  config.set('Impacts', 'output', 'pH')
else:
  data=re.sub('#6[^>]+#7', '', data)
  data2=re.sub('#6[^>]+#7', '', data2)
  OUTPUTS=OUTPUTS+'\n\nWARNING: IMPACT ANALYSIS IS TURNED OFF!'

#if report
if repbut_v1.get() == "True":
  if ThrdirBox6.get() !='':
    OUTPUTS=OUTPUTS+'\n\nThe directory for the reporting input data has been set as '+ThrdirBox6.get()+'.'
    data = data.replace("docker run -it $options --mount type=bind,source=\"$(pwd)\",target=/srv/actom-output/input", "docker run -it $options --mount type=bind,source="+ThrdirBox6.get()+",target=target=/srv/actom-output/input")
    data2 = data2.replace("docker run -it $options --mount type=bind,source=\"%cd%\",target=/srv/actom-output/input", "docker run -it $options --mount type=bind,source="+ThrdirBox6.get()+",target=target=/srv/actom-output/input")
else:
  data=re.sub('#7[^>]+#8', '', data)
  data2=re.sub('#7[^>]+#8', '', data2)
  OUTPUTS=OUTPUTS+'\n\nREPORTING IS TURNED OFF!'

print(OUTPUTS)

data=data.replace("mkdir -p ../","mkdir -p ",1)
data=data.replace("cd ../","cd ",1)

data2=data2.replace("mkdir ../","mkdir -p ",1)
data2=data2.replace("cd ../","cd ",1)

if c4=='Yes':
   data=data.replace("$options","")
   data2=data2.replace("$options","")
   config.set('General', 'options', '')
else:
   data=data.replace("$options","-a STDERR")
   data2=data2.replace("$options","-a STDERR")
   config.set('General', 'options', '-a STDERR')

with open('input/data.ini', 'w') as configfile:
    config.write(configfile)

data2=data2.replace("#","@ Rem ")
data2=data2.replace("-p","")
data2=data2.replace("$(pwd)","%cd%")
data2=data2.replace("-it","-i")

first_line="#!/bin/bash\n"
with open('input/Run-All.sh', 'w') as modified: modified.write(first_line+data)
with open('input/Run-All.bat', 'w') as modified2: modified2.write(data2)

print('\nThe run files have been generated...\n\n')

