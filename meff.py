#!/usr/bin/env python3.8

import h5py
import sys

import scipy
import numpy as np
from numpy import ndarray
import subprocess
import re



from scipy import linalg
from scipy.linalg import eigh
import os

if len(sys.argv)!=3:
    print("<corr> <name>")
filename = sys.argv[1]
name = sys.argv[2]

def read_file(fname):

    f = open(fname,"r")
    lines=f.readlines()
    Ncfgs = int(lines[0].split(' ')[0])
    Ntslices = int(lines[0].split(' ')[1])
    Corr = [ [] for t in range(Ntslices)]
    c=0
    for i in range(len(lines)):

        if i!=0:

            Corr[c%Ntslices].append(float(lines[i].split(' ')[1].strip()))
            c+=1
    f.close()
    return Corr,Ncfgs,Ntslices

def make_jack(Corr,Corravg,Ntslices,Ncfgs):
    
    Corr_jack = [[] for i in range(Ntslices)]
    
    for t in range(Ntslices):
        for cfg in range(len(Corr[t])):
            
            Corr_jack[t].append((Ncfgs/(Ncfgs-1))*Corravg[t]-(1/(Ncfgs-1))*Corr[t][cfg])

    return Corr_jack
def un_jack(Corrjack,Corravg,Nt,Ncfgs):

    Corr_unjack = [[] for j in range(Nt)]
    print(len(Corrjack[1]),len(Corravg))
    for t in range(0,Nt-1):
        for cfg in range(Ncfgs):
            if t ==0:
                Corr_unjack[t].append(0.0)
            else:
                Corr_unjack[t].append(Ncfgs*Corravg[t]-(Ncfgs-1)*Corrjack[t][cfg])
    return Corr_unjack

def main():

    Nt = 29
    Correlator = read_file(filename)[0]
    Ntslices = read_file(filename)[2]
    Ncfgs = read_file(filename)[1]
    Correlator_avg = []
    for t in range(Ntslices):
        print(Correlator[t])
        Correlator_avg.append(np.mean(Correlator[t]))

    Corr_jack = make_jack(Correlator,Correlator_avg,Ntslices,Ncfgs)
    meff_jack = [[] for i in range(Nt)]
    for cfg in range(Ncfgs):

        for t in range(1,Nt):

            meff_jack[t].append(np.log(Corr_jack[t+3][cfg]/Corr_jack[t][cfg])/(-3))
    meff_avg=[]
    for t in range(Nt-1):
        
        meff_avg.append(np.mean(meff_jack[t]))
    meff = un_jack(meff_jack,meff_avg,Nt,Ncfgs)

    f = open(name+".dat","w")

    f.write(str(Ncfgs)+" "+str(Nt)+" 0 0 1\n")

    for cfg in range(Ncfgs):
        for t in range(Nt-1):

            f.write(str(t)+" "+str(meff[t][cfg])+"\n")

    f.close()
if __name__ == "__main__":
  main()

