#!/usr/bin/env python3.8

import h5py
import sys

sys.path.append("../")
sys.path.append("./")
print(sys.path)
import scipy
import numpy as np
from numpy import ndarray
import subprocess
import re

import fit_defs 

from fit_defs import * 
from scipy import linalg
from scipy.linalg import eigh
import os


if len(sys.argv)!=3:

   print("usage: <h5 file> <ops.list>")

   exit(1)

h5_file = sys.argv[1]
ini_list = sys.argv[2]
ops_list, irrep = eigen_solve.read_input_ops(ini_list)

def main():
    
    Corr_real, Corr_imag = eigen_solve.get_corr_matrix(h5_file,irrep, ops_list) 
    
    ########## apply weights to Correlator ##############


    

    print("Checking Imaginary part")

    eigen_solve.check_imag_part(Corr_real,Corr_imag)
    Corr_jack = jack_utils.scale_up_corr(Corr_real)
   
    ########## Solve GEVP ###############################
    t0s=[2,3,4,5]
    for t0 in t0s:
      eigs_not_ordered , vecs_not_ordered = eigen_solve.gen_eigen_solver(Corr_jack,t0)
      eigs_ordered, vecs_ordered = eigen_solve.order_eigs_and_vecs(eigs_not_ordered, vecs_not_ordered)  


      ##### get true eigenvectors #########

    
      Ncfgs, Nt, op_num = eigs_ordered.shape[0], eigs_ordered.shape[1], eigs_ordered.shape[2]
      
      vecs_trans = eigen_solve.transform_vecs(Corr_jack,t0,vecs_ordered)
      #vecs_trans = vecs_ordered
      prin_corrs = np.zeros(eigs_ordered.shape)
      for op in range(op_num):
        for t in range(Nt):
            prin_corrs[:,t,op] = jack_utils.scale_down(eigs_ordered[:,t,op])
      vecs = np.zeros(vecs_trans.shape)
      
      for op1 in range(op_num):
        for op2 in range(op_num):
          for t in range(Nt):
            vecs[:,t,op1,op2] = jack_utils.scale_down(vecs_trans[:,t,op1,op2])
    #        print(vecs_trans[:,t,op1,op2].shape)
   #         vecs[:,t,op1,op2] = vecs_trans[:,t,op1,op2]

     
     # Zs = eigen_solve.Zs(vecs,t0,0.5)
      ######### test #####################
      cfg_test, t_test, state = 13, 10, 2
      #print(vecs_trans[cfg_test,t_test,:,state])
      #print(Corr_jack[cfg_test,:,:,t_test].dot(vecs_trans[cfg_test,t_test,:,state]))
      #print(eigs_ordered[cfg_test,t_test,state]*Corr_jack[cfg_test,:,:,t0].dot(vecs_trans[cfg_test,t_test,:,state]))
      



      ################ fit prin_corrs ####################
      subprocess.run(f"mkdir t0_{t0}", shell=True)
      
     
    
      
   
  
   #   for i in range(op_num):
   #       params = fit_defs.fit_(prin_corrs[:,:,i],"double_exp",t0,tmin,tmax,svd_cutoff)[0]
   #       write_jacks.param("t0_"+str(t0)+"/masses/mass_state_"+str(i)+".jack",params[1])
      
      
  ######### write results #############
   
      subprocess.run(f"mkdir t0_{t0}/principle_correlators", shell=True)

      for op in range(op_num):
        write_jacks.t_("t0_"+str(t0)+"/principle_correlators/principle_correlator_"+str(op)+".jack",prin_corrs[:,:,op])
      ############# run principal correlator fits and store fit logs ###################

      subprocess.run(f"mkdir t0_{t0}/prin_corr_fit_logs", shell=True)

      direct = "t0_"+str(t0)+"/prin_corr_fit_logs"
      os.chdir(direct)

      tmin = t0+1
      tmax = 26
      tslicenum = 5
    #  noise_cut_off = 0.3
      svd_cut_off = 1e-6
    #  fit_crit = "generic"
      for i in range(0,op_num):

        name = "../principle_correlators/principle_correlator_"+str(i)+".jack"
        subprocess.run(f"/Users/christopherjohnson/projects/GEVP/fitters/fitter {name} {t0} {tmin} {tmax} {tslicenum} {svd_cut_off}", shell=True)
        fitname = "principle_correlator_"+str(i)+".fit_info"
        subprocess.run(f"mv ../principle_correlators/{fitname}* .", shell=True)
    
      ############## copy fit results ################################
      os.chdir("../../")
      subprocess.run(f"mkdir t0_{t0}/mass_jacks", shell=True)
      direct = "t0_"+str(t0)+"/mass_jacks"
      os.chdir(direct)
      for i in range(0,op_num):
        mass_name = "mass_state_"+str(i)+".jack"
        fit_jack_name = "../principle_correlators/principle_correlator_"+str(i)+"_mass.jack"
        subprocess.run(f"mv {fit_jack_name} {mass_name}", shell=True)
  #    os.chdir("../../")
        #######################################################
   #   mass_avg = []
   #   for i in range(0,op_num):
   #       mass_avg.append(np.mean(read_jacks.param("t0_"+str(t0)+"/mass_jacks/mass_state_"+str(i)+".jack"),axis = 0))
      
   #   Zs = eigen_solve.Zs(vecs,t0,mass_avg)
      subprocess.run(f"mkdir t0_{t0}/eigenvectors", shell=True)
      for op in range(op_num):
        for row in range(op_num):
          write_jacks.t_("t0_"+str(t0)+"/eigenvectors/eigenvector_state_"+str(op)+"_op_"+str(row)+".jack",vecs[:,:,row,op])      

   #   subprocess.run(f"mkdir t0_{t0}/Zs", shell=True)
   #   for op in range(op_num):
   #     for row in range(op_num):
   #       write_jacks.t_("t0_"+str(t0)+"/Zs/Z_state_"+str(op)+"_op_"+str(row)+".jack",Zs[:,:,row,op])
if __name__=="__main__":
   main() 
