# script to copy data over to experiment data

import shutil, os

files = ['input.xml','out.txt']

dest_folder = '/ascldap/users/wmreese/Documents/Experiment_Data/ip_smoothed_T_4_nt_25_alpha_1e-2' 
 

for f in files:
    shutil.copy(f, dest_folder)
