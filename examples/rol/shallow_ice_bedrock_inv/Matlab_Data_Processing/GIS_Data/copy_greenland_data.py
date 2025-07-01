# this is a script to move the data generated from a subset of Greenland to Data_Generation

# directions:
#1. make changes in Process.m, make sure that width, height, nx/ny are correct 

#2. copy all txt files to
# /scratch/Trilinos/build_opt/packages/rol/example/HdsaLib/examples/shallow_ice_bedrock_inv/Data_Generation
import shutil, os

files = ['Bedrock_Topography.txt', 'Forcing.txt','Log_Basal_Sliding.txt', 'Surface_Height.txt']

dest_folder = '/Users/joshart/software/Trilinos/build_opt/packages/rol/example/rol/shallow_ice_bedrock_inv/Data_Generation'
for f in files:
    shutil.copy(f, dest_folder)

