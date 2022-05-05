#!/usr/bin/env python3.8
import scipy
import h5py
import sys
import numpy as np
from numpy import ndarray
import subprocess
import re

import pathlib
import time
from scipy import linalg
from scipy.linalg import eigh
import os

states = []
for i in range(2,len(sys.argv)):
    states.append(int(sys.argv[i]))

if len(sys.argv)<2:

    print("usage : t0 state0 state1 ... ")
    exit(1)


t0 = sys.argv[1]


chisqs = []
masses = []
masserrs = []
##### get masses and chisq ###############
for state in states:

    string = t0+"/prin_corr_fit_logs/prin_correlator_"+str(state)+".jack_prin_corr_exp_fit.log"
    tmp = subprocess.check_output("cat %(string)s | awk '{print $6}' | head -2 | tail -1" % locals(),shell=True).decode(sys.stdout.encoding)
    chisq = tmp[:-2]
    chisqs.append(chisq)
    mass = subprocess.check_output("cat %(string)s | awk '{print $10}' | head -2 | tail -1" % locals(),shell=True).decode(sys.stdout.encoding)
    masses.append(mass)
    tmp = subprocess.check_output("cat %(string)s | awk '{print $12}' | head -2 | tail -1" % locals(),shell=True).decode(sys.stdout.encoding)
    masserrs.append(tmp[0:6])



gnu = open("prin_tmp.gnu","w")




gnu.write("unset key \n")
#gnu.write("set multiplot \n")

num = len(states)
if num < 4:

   xsize = 1/num
   ysize = 1/num
else:
    xsize = 0.3
    ysize =0.3
cx=0
cy=0

for i in range(0,int(len(states))):
    ever1 = num * i
    ever2 = num * i + num-1
    
    
    yorigin = 1 - ysize
    if (cx*xsize + xsize) > 1:

        yorigin -= ysize
        cx=0
        cy-= ysize
    yorigin = 1 -ysize+cy
    xorigin = cx*xsize
    cx+=1
    gnu.write("set origin "+str(xorigin)+","+str(yorigin)+"\n")
    gnu.write("set size "+str(xsize)+","+str(ysize)+"\n")
    string = t0+"/prin_corr_fit_logs/prin_correlator_"+str(states[i])+".jack_prin_corr_exp_fit.ax"
    gnu.write("set label '{/Symbol c}^2/N_{dof} ="+str(chisqs[i])+"' at 20,1.21 font 'Times,8' \n")
    gnu.write("set label 'm ="+str(masses[i].strip())+" +/- "+str(masserrs[i])+"' at 20,1.15 font 'Times,8' \n")
    gnu.write("set term 'png' \n")
    gnu.write("set output 'prin.png'\n")
    tmp = "plot '"+string+"' index 3 u 1:2 w lines lt 6,\\\n"
    tmp = tmp +"     '"+string+"' index 4 u 1:2 w lines lt 6 ,\\\n"
    tmp = tmp +"     '"+string+"' index 5 u 1:2 w lines lt 6 ,\\\n"
    tmp = tmp +"     '"+string+"' index 0 u 1:2 w lines lt 4 ,\\\n"
    tmp = tmp +"     '"+string+"' index 1 u 1:2 w lines lt 4 ,\\\n"
    tmp = tmp +"     '"+string+"' index 2 u 1:2 w lines lt 4 ,\\\n"
    tmp = tmp +"     '"+string+"' index 6 u 1:2 w lines lt 4 ,\\\n"
    tmp = tmp +"     '"+string+"' index 7 u 1:2 w lines lt 4 ,\\\n"
    tmp = tmp +"     '"+string+"' index 8 u 1:2 w lines lt 4 ,\\\n"
    tmp = tmp +"     '"+string+"' index 9 u 1:2:3 w yerr ls 6 ,\\\n"
    tmp = tmp +"     '"+string+"' index 10 u 1:2:3 w yerr ls 4 ,\\\n"



    gnu.write(tmp+"\n")   
subprocess.Popen("gnuplot -persist prin_tmp.gnu", shell=True)
#time.sleep(2)
#os.system("rm mult_tmp*")
