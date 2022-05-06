#!/usr/bin/env python3.8
import scipy
import h5py
import sys
sys.path.append("../data_conversion")
import numpy as np
from numpy import ndarray
import subprocess
import re

import pathlib
import time
from scipy import linalg
from scipy.linalg import eigh
import os

t0s = []
for i in range(2,len(sys.argv)):
    t0s.append(int(sys.argv[i]))

if len(sys.argv)<2:

    print("usage : states_num t01 t02 ... ")
    exit(1)
masses = [ []  for i in range(len(t0s)) ]
masserrs = [ []  for i in range(len(t0s)) ]
state_nums = sys.argv[1]
for t0 in range(0,len(t0s)):

    tmp = "t0_"+str(t0s[t0])+"/mass_jacks/"
    
    
    for i in range(0,int(state_nums)):
        state = tmp+"/mass_state_"+str(i)+".jack"
        t = subprocess.check_output("calc %(state)s | awk '{print $2}'" % locals(),shell=True).decode(sys.stdout.encoding)
        m = t.strip()
        t = subprocess.check_output("calc %(state)s | awk '{print $3}'" % locals(),shell=True).decode(sys.stdout.encoding)
        merr = t.strip()
        masses[int(t0)].append(m)
        masserrs[t0].append(merr)

p = open("mult_tmp.dat","w")

for i in range(0,len(masses[0])):

    for t in range(len(t0s)):

        tmp = str(t0s[t])+" "+str(masses[t][i])+" "+str(masserrs[t][i])

        p.write(tmp+"\n")

p.close()

gnu = open("mult_tmp.gnu","w")
xmin = min (t0s)-0.5
xmax = max(t0s)+0.5
gnu.write("set xrange ["+str(xmin)+":"+str(xmax)+"]\n")
gnu.write("set xtics nomirror "+str(min(t0s))+",1,"+str(max(t0s))+"\n")
gnu.write("unset key \n")
gnu.write("set multiplot \n")
xsize = 0.3
ysize = 0.3
cx=0
cy=0
t0num = len(t0s)
for i in range(0,int(state_nums)):
    ever1 = t0num * i
    ever2 = t0num * i + t0num-1
    
    
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
    gnu.write("plot 'mult_tmp.dat' every ::"+str(ever1)+"::"+str(ever2)+" u 1:2:3 w yerr\n")
subprocess.Popen("gnuplot -persist mult_tmp.gnu", shell=True)
#time.sleep(2)
#os.system("rm mult_tmp*")
