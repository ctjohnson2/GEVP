#!/usr/bin/env python3.8

import pandas as pd
import numpy as np
import pylab as plt
from fit_defs import eigen_solve, jack_utils, read_jacks
import sys
state = sys.argv[1]
Ncfgs,Nt,ops = 198,41,9
vecs = np.zeros((Ncfgs,Nt,ops))
for i in range(ops):
    vecs[:,:,i] = read_jacks.prin_corr("t0_2/eigenvectors/eigenvector_state_"+state+"_op_"+str(i)+".jack") 

t_test = 8
vecs_new =np.mean(vecs[t_test,:,:],axis=0)
vecs_err = np.var(vecs[t_test,:,:],axis=0)/Ncfgs


vecs_new = abs(vecs_new)
norm = np.amax(vecs_new)
vecs_new = abs(vecs_new) /norm

vecs_err  =  np.sqrt(vecs_err) / norm
  
labels=[]
for i in range(ops):
    labels.append(i)


x_pos = np.arange(len(labels))


fig, ax = plt.subplots()
ax.set_ylim([0,1.3])
ax.bar(x_pos, vecs_new, yerr=vecs_err)
plt.show()


