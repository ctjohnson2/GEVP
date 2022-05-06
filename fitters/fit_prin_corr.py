#!/usr/bin/env python3.8


from ../fit_defs import *

if len(sys.argv)!=6:
    print("Usage: <file> <t0> <tmin> <tmax> <svd_cutoff>")
    exit(1)

prin_corr_file = sys.argv[1]
t0 = int(sys.argv[2])
tmin = int(sys.argv[3])
tmax = int(sys.argv[4])

svd_cutoff = float(sys.argv[5])


 
  
def main():

   prin_corr = read_jacks.prin_corr(prin_corr_file)
   print(prin_corr)
   result = fit_(prin_corr,"double_exp",t0,tmin,tmax,svd_cutoff)

   write_jacks.param("mass.jack",result[0][1])
if __name__== "__main__":
   main()
