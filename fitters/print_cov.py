#!/usr/bin/env python3.8

import sys
sys.path.append("../")
import fit_defs
from fit_defs import *


if len(sys.argv)!=4:
    print("<prin_corr_file> <tmin> <tmax>")
    exit(1)
prin_corr_file = sys.argv[1]
tmin = int(sys.argv[2])
tmax = int(sys.argv[3])




def main():

   prin_corr = read_jacks.prin_corr(prin_corr_file)
   print(jack_utils.my_cov(prin_corr,tmin,tmax) )
   

   
if __name__== "__main__":
   main()
