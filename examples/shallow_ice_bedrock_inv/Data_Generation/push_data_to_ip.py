# script to copy data over to experiment data

import shutil, os

files = ['Bedrock_Topography.txt','Forcing.txt','Log_Basal_Sliding.txt','Surface_Height.txt','nodes.txt','cell_to_node_quad.txt']

dest_folder = '/scratch/Trilinos/build_opt/packages/rol/example/HdsaLib/examples/shallow_ice_bedrock_inv'
for f in files:
	shutil.copy(f, dest_folder)

shutil.copy('clean_true_state.txt', '/scratch/Trilinos/build_opt/packages/rol/example/HdsaLib/examples/shallow_ice_bedrock_inv/clean_true_state.txt')

shutil.copy('coarse_data.txt', '/scratch/Trilinos/build_opt/packages/rol/example/HdsaLib/examples/shallow_ice_bedrock_inv/coarse_data.txt')
