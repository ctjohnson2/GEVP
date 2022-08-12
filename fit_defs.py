#!/usr/bin/env python3.8

import numpy as np
import math
from scipy.optimize import minimize
import h5py
from scipy import linalg
from scipy.linalg import eigh

class eigen_solve():
    
    def read_input_ops(in_list):
      r = open(in_list,"r")
      lines = r.readlines()
      op_list=[]
      for line in lines:
        op_list.append((line.split(' || ')[1]).strip())


      irr=lines[0].split(' || ')[0]
      r.close()
      return op_list,irr

    def get_corr_matrix(h5_file,irrep, ops_list):

      with h5py.File(h5_file,"r") as f:
        found_irrep = False
        print("Available keys: ", f.keys())
        for key in f.keys():
          if key==irrep:
            print(irrep+" found key")
            print(f[irrep].keys())
            found_irrep = True

        if found_irrep == False:
          print("Can't find irrep "+irrep+" in the key list")
          print(f.keys())
          exit(1)
     
        total_ops = (f[irrep]).attrs['op_list']
        print("Total available operators: ",total_ops)
  

        # make ops map 
        ops_map=[] 
        op_num = len(ops_list)
        for i in range(op_num):
          c=False
          for j in range(len(total_ops)):
            if ops_list[i]==total_ops[j]:
              c=True
              ops_map.append([i,j])

          if c==False:
            print("Operator :",ops_list[i]," not in basis")
            exit(1)

     
        if(len(ops_list)==1 or len(total_ops)==1):
          print("Can't do GEVP on a 1x1 matrix")
          data = (f[irrep])['data']
          Ncfgs = data.shape[0]
        
          Ntslices = data.shape[-1]
          Creal = []
          p = open("corr.1x1.dat","w")
          p.write(str(Ncfgs)+" "+str(Ntslices)+" 0 0 1\n")
          print("Number Configs",Ncfgs)
          for cfg in range(0,Ncfgs):

            for t in range(0,Ntslices):
                 if len(total_ops)!=1:
                  
                  string = str(t)+" "+str(data[cfg][ops_map[0][1]][ops_map[0][1]][t].real)
                 else:
                  
                  string = str(t)+" "+str(data[cfg][t].real)
                 p.write(string+"\n")
         
          p.close()
          print("Wrote correlator into corr.1x1.dat")
          exit(1)


        data = (f[irrep])['data']
        Ncfgs = data.shape[0]
        Ntslices = data.shape[3]

        Creal=np.zeros(data.shape)
        Cimag=np.zeros(data.shape)
        print(f"Loading Correlator matrix",flush=True)

        #for cfg in range(0,Ncfgs):
        for i in range(op_num):
          for j in range(op_num):
               total = data[:,ops_map[i][1],ops_map[j][1],:].real
               totalim = data[:,ops_map[i][1],ops_map[j][1],:].imag
               zeros = np.zeros((Ncfgs,Ntslices))
               if np.array_equal(total, zeros)==True and np.array_equal(totalim,zeros)==True:
                  #### checking if symmetric part is nonzero ####
                  total = data[:,ops_map[j][1],ops_map[i][1],:].real
                  totalim = data[:,ops_map[j][1],ops_map[i][1],:].imag

               Creal[:,i,j,:] = total
               Cimag[:,i,j,:]  = totalim


        #Creal = Creal.reshape(Ncfgs,op_num,op_num,Ntslices)
        #Cimag = Cimag.reshape(Ncfgs,op_num,op_num,Ntslices)



      f.close()
      return Creal, Cimag 
  

    def gen_eigen_solver(Corr_jack,t0):
      print("Solving GEVP for t0="+str(t0))
      eigs_not_ordered_all,  vecs_not_ordered_all =[], []
      Ncfgs, op_num, Ntrange = Corr_jack.shape[0], Corr_jack.shape[1], Corr_jack.shape[3]
      ##### solve for each cfg

      for cfg in range(0,Ncfgs):

        Ct0 = Corr_jack[cfg,:,:,t0]
        
        L = np.linalg.cholesky(Ct0)
        Ldag = L.conj().T
        Linv = np.linalg.inv(L)
        Ldag_inv = np.linalg.inv(Ldag)
        ###### solve for each time slice
        eigs_not_ordered, vecs_not_ordered = [], []

        for t in range(0,Ntrange):

            B = np.matmul(Linv,np.matmul(Corr_jack[cfg,:,:,t],Ldag_inv))
            eigs, vecs = eigh(B)

            eigs_not_ordered = np.append(eigs_not_ordered,eigs)
            vecs_not_ordered = np.append(vecs_not_ordered,vecs)

        eigs_not_ordered_all = np.append(eigs_not_ordered_all,eigs_not_ordered)
        vecs_not_ordered_all = np.append(vecs_not_ordered_all,vecs_not_ordered)

      eigs_not_ordered_all = eigs_not_ordered_all.reshape(Ncfgs,Ntrange,op_num)
      vecs_not_ordered_all = vecs_not_ordered_all.reshape(Ncfgs,Ntrange,op_num,op_num)

      return eigs_not_ordered_all,vecs_not_ordered_all

    def check_imag_part(Corr_real,Corr_imag):
      Corravg_re=np.mean(Corr_real,axis=0)
      Corravg_im=np.mean(Corr_imag,axis=0)
      Ncfgs, op_num, Ntslices = Corr_real.shape[0], Corr_real.shape[2], Corr_real.shape[3]
      var_im=np.var(Corr_real,axis=0)/Ncfgs

      ############### check imaginary part of correlation matrix ################################     
#will want this at somepoint but it is taking too long
      for t in range(1,Ntslices):
        for i in range(0,op_num):
          for j in range(0,op_num):
            if abs(Corravg_im[i][j][t])> 4.0*(np.sqrt(var_im[i][j][t])):

             dis_from_zero = abs(Corravg_im[i][j][t])/np.sqrt(var_im[i][j][t])

            ### complain if matrix has significant imaginary part

             print("time = "+str(t)+" :op "+str(i)+" "+str(j)+" large imaginary part = "+str(dis_from_zero)+" x sigma")
             print("theta ="+str(np.arctan(Corravg_im[i][j][t]/Corravg_re[i][j][t])))


    def order_eigs_and_vecs(eigs_not_ordered_all,vecs_not_ordered_all):
      
      ############ order the eigenvectors and eigenvalues ##################################
      eigs_ordered_all = np.zeros(eigs_not_ordered_all.shape)
      vecs_ordered_all = np.zeros(vecs_not_ordered_all.shape)
      Ncfgs, op_num, Ntrange = eigs_not_ordered_all.shape[0], eigs_not_ordered_all.shape[2], eigs_not_ordered_all.shape[1]
      for cfg in range(0,Ncfgs):
        for t in range(0,Ntrange): 
          if t > -1:
            vstar_ref = vecs_not_ordered_all[0,t,:,:].conj()
            for i in range(0,op_num):
              norm = 0
              opnum = 0
              for j in range(0,op_num):
               
                dot = abs(vstar_ref[:,i].T.dot(vecs_not_ordered_all[cfg,t,:,j]))
            
             
                if dot > norm:
                  norm = dot
                  opnum = j
              
              if vecs_not_ordered_all[cfg,t,0,opnum]<0:
                vecs_not_ordered_all[cfg,t,:,opnum]*=-1
              vecs_ordered_all[cfg,t,:,opnum] = vecs_not_ordered_all[cfg,t,:,opnum]
              eigs_ordered_all[cfg,t,opnum] = eigs_not_ordered_all[cfg,t,opnum]
              
      return eigs_ordered_all, vecs_ordered_all
    def transform_vecs(Corr_jack,t0,vecs_ordered):
    
      Ncfgs, Nt, op_num = Corr_jack.shape[0], Corr_jack.shape[3], Corr_jack.shape[2]
      vecs_trans = np.zeros(vecs_ordered.shape)
 
      for cfg in range(Ncfgs):
        Ct0 = Corr_jack[cfg,:,:,t0]

        L = np.linalg.cholesky(Ct0)
        Ldag = L.conj().T
        Ldag_inv = np.linalg.inv(Ldag)
        for t in range(Nt):
          for op in range(op_num):
            vecs_trans[cfg,t,:,op] = Ldag_inv.dot(vecs_ordered[cfg,t,:,op])
        
            #vecs_trans[cfg,t,op,:] = vecs_trans[cfg,t,op,:].conj().dot(Ldag)
      return vecs_trans
    
    def Zs(vecs_trans,t0,m):
        Zs = np.zeros(vecs_trans.shape)
        Ncfgs, Nt , ops= vecs_trans.shape[0], vecs_trans.shape[1], vecs_trans.shape[2]
        for cfg in range(Ncfgs):
            for t in range(Nt):
                Zs[cfg,t,:,:] = np.linalg.inv(vecs_trans[cfg,t,:,:])
        for state in range(ops):
            Zs[:,:,state,:] = Zs[:,:,state,:]*np.sqrt(2*m[state])*np.exp(m[state]*t0/2)
        return Zs

class jack_utils():
###################### jackknife utility functions #####################
    def scale_up(unjack):
    ##### create jackknife enseble ##########
      Ncfgs = len(unjack)
      avg =np.mean(unjack,axis=0)
      jacked = []
      for cfg in range(Ncfgs):

        tmp = (Ncfgs/(Ncfgs-1))*avg-(1/(Ncfgs-1))*unjack[cfg]

        jacked=np.append(jacked,tmp)
      return jacked;

    def scale_up_corr(unjack):
    ##### create jackknife enseble ##########
      Ncfgs = unjack.shape[0]
      avg =np.mean(unjack,axis=0)
      jacked = []
      for cfg in range(Ncfgs):

        tmp = (Ncfgs/(Ncfgs-1))*avg-(1/(Ncfgs-1))*unjack[cfg]

        jacked=np.append(jacked,tmp)
      return jacked.reshape(unjack.shape);
    def scale_down(jack):
    #### undo jackknife bias ################

      Ncfgs = len(jack)
      avg = np.mean(jack)
      unjacked=[]
      for cfg in range(Ncfgs):

        tmp = Ncfgs*avg-(Ncfgs-1)*jack[cfg]
        unjacked = np.append(unjacked,tmp)
      return unjacked

    def my_cov(prin_corr, tmin, tmax):
        ############# inverse covariance matrix ########################
        Ncfgs = prin_corr.shape[0]
        
        Cov = np.zeros((tmax-tmin,tmax-tmin))
        avg,var =[],[]
        for t1 in range(tmin,tmax):
            avg = np.append(avg,np.mean(prin_corr[:,t1]))
       
        for t1 in range(tmax-tmin):
          for t2 in range(tmax-tmin):
            cov = 0
            for cfg in range(Ncfgs):
              
              cov+=(prin_corr[cfg,t1+tmin]-avg[t1])*(prin_corr[cfg,t2+tmin]-avg[t2])
     
            Cov[t1,t2] = cov/(Ncfgs*(Ncfgs-1))
        return Cov
    def mycinv(prin_corr, tmin, tmax, svd_cutoff):
        Cov = my_cov(prin_corr, tmin, tmax)
        cinv = np.linalg.pinv(Cov, rcond = svd_cut_off)
        
        return cinv
class read_jacks():

    def prin_corr(prin_corr_name):
      f = open(prin_corr_name,"r")
      lines = f.readlines()
      Ncfgs, Nt = int(lines[0].split(' ')[0]),int(lines[0].split(' ')[1])
      prin_corr = np.zeros((Ncfgs,Nt))
      t,c=0,0
      for i in range(len(lines)):
        if i!=0:

         prin_corr[c,t%Nt] = float(lines[i].split(' ')[1].strip())
         t+=1
         if t%Nt ==0 and t!=0:
          c+=1
      f.close()
      return prin_corr
    def param(mass_name):
      f = open(mass_name,"r")
      lines = f.readlines()
      Ncfgs = int(lines[0].split(' ')[0])
      masses = [] 
      for i in range(len(lines)):
        if i!=0:
          masses = np.append(masses, float(lines[i].split(' ')[2]))
      return masses

         
class write_jacks():

    def param(file_name, param_jack):

      Ncfgs = len(param_jack)
      f = open(file_name,"w")
      f.write(str(Ncfgs)+" 1 0 0 1\n")
      for cfg in range(Ncfgs):

       f.write(" 0 "+str(param_jack[cfg])+"\n")
      
      f.close()

    def t_(file_name, prin_corr_jack):

      Ncfgs, Nt = prin_corr_jack.shape[0], prin_corr_jack.shape[1]
      f = open(file_name,"w")
      f.write(str(Ncfgs)+" "+str(Nt)+" 0 0 1\n")
      for cfg in range(Ncfgs):
        for t in range(Nt):
          f.write(str(t)+" "+str(prin_corr_jack[cfg,t])+"\n")

      f.close()



class fit_functions():
  

    def lambda_(x,tmin,tmax,t0,bin_size,type_):
     ########## (1-A) * exp(-m ( t - t0) ) + A * exp(-m' (t - t0) )
      re_= np.zeros((bin_size,tmax-tmin))
      for t in range(tmin,tmax):
       for cfg in range(bin_size): 
        if type_=="single_exp":
          thing = x[0]*np.exp(-x[1]*(t-t0))
        elif type_=="double_exp":
          thing = (1-x[0])*np.exp(-x[1]*(t-t0))+x[0]*np.exp(-x[2]*(t-t0))
        else:
          print("Acceptable fit types are : single_exp, double_exp")
          exit(1)
        re_[cfg,t-tmin] = thing

      
      return re_


def fit_many_stochastic(prin_corr,t0,tmin,tmax,tslices_min,svd_cutoff):

  ########## first make smaller random sample ##########

  Ncfgs, Nt = prin_corr.shape[0], prin_corr.shape[1]
  
  N_st = math.ceil(Ncfgs*0.2)
 
  prin_corr_stoch = np.zeros((N_st,Nt))
  print(prin_corr.shape,prin_corr_stoch.shape)
  indexes = []
  for r in range(N_st):
    random_index = np.random.randint(Ncfgs)
    prin_corr_stoch[r,:] = prin_corr[random_index,:]
    indexes.append(random_index)
  chisq = 0
  tmax_w, tmin_w = 0,0
  for t_sl in range(tslices_min,tmax-tmin):
    for t in range(tmin,tmax):
      if t + t_sl < tmax:

        fit_(prin_corr_stoch,"double_exp",t0,t,t+t_sl,svd_cutoff)
        print(fit_(prin_corr_stoch,"double_exp",t0,t,t+t_sl,svd_cutoff))
        if fit_(prin_corr_stoch,"double_exp",t0,t,t+t_sl,svd_cutoff)[1] > chisq:
          chisq = fit_(prin_corr_stoch,"double_exp",t0,t,t+t_sl,svd_cutoff)
          tmin_w=t
          tmax_w= t+t_sl
  return fit_(prin_corr,"double_exp",t0,tmin_w,tmax_w,svd_cutoff)
def fit_(prin_corr,lambda_type,t0,tmin,tmax,svd_cutoff):
     
     
  debug_minimize=True
  Ncfgs, T_range = prin_corr.shape[0], prin_corr.shape[1]
  
  prin_corr_jack=np.zeros(prin_corr.shape)
  for t in range(T_range):
      prin_corr_jack[:,t] = jack_utils.scale_up(prin_corr[:,t])
  
  cinv = jack_utils.my_cinv(prin_corr, tmin, tmax, svd_cutoff)
  
  t_fit_range = tmax - tmin
  if lambda_type == "single_exp":
    x0 =[.1,1.]
  elif lambda_type == "double_exp":
    x0 =[.1,0.25,1.5]
  Ndof = tmax-tmin-len(x0)
  fit_params = [ [] for i in range(len(x0))]
  if len(x0)==3:
      def constraint(x):
        return  x[2] - x[1]
  for cfg in range(Ncfgs):
    
    
    chisq = lambda x:  (1/Ndof)*(((fit_functions.lambda_(x,tmin,tmax,t0,Ncfgs,lambda_type)[cfg,:] - prin_corr_jack[cfg,tmin:tmax]).T.dot(cinv.dot(fit_functions.lambda_(x,tmin,tmax,t0,Ncfgs,lambda_type)[cfg,:]-prin_corr_jack[cfg,tmin:tmax]))))
         
            
         
    res = minimize(chisq,x0,method='CG', tol=1e-2) 
    x0 = res.x
    if(debug_minimize):
      print(res,chisq(res.x))
    if res.success==False:

       print("WARNING MINIMIZE FAILED FOR CFG NUMBER : ",cfg)
       print(res)
                     
    for i in range(len(x0)):
      fit_params[i].append(res.x[i])
   
  centrals = [np.mean(fit_params[0]),np.mean(fit_params[1]),np.mean(fit_params[2])]
  print(centrals)
  chi = 0
  for cfg in range(Ncfgs):

      chi+=  ((fit_functions.lambda_(centrals,tmin,tmax,t0,Ncfgs,lambda_type)[cfg,:] - prin_corr[cfg,tmin:tmax]).T.dot(cinv.dot(fit_functions.lambda_(centrals,tmin,tmax,t0,Ncfgs,lambda_type)[cfg,:]-prin_corr[cfg,tmin:tmax])))/Ndof
  print("CHISQ =",chi/Ncfgs," for fit type: ",lambda_type," tmin: ",tmin," tmax: ",tmax)
  
  end_params = []
  for i in range(len(x0)):
    end_params.append(jack_utils.scale_down(fit_params[i]))
  return(end_params,chi/Ncfgs)
       

def main():
       
   prin_corr = read_jacks.prin_corr("principle_correlator_7.jack")
   print(prin_corr)
   result = fit_(prin_corr,"double_exp",3,4,20,1.0e-6) 

   write_jacks.param("m.jack",result[0][1])
if __name__== "__main__":
   main()
