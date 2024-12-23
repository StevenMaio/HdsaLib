# script to copy data over to experiment data

import shutil, os

files = ['input_data_generation.xml','out.txt']

dest_folder = '/ascldap/users/wmreese/Documents/Experiment_Data/fp_smoothed_T_4_nt_25' 
for f in files:
    shutil.copy(f, dest_folder)
