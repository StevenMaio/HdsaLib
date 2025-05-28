# script to copy data over to experiment data

import shutil, os

files = ['Bedrock_Topography.txt','Forcing.txt','Log_Basal_Sliding.txt','Surface_Height.txt','nodes.txt','cell_to_node_quad.txt']

dest_folder = '../'
for f in files:
	shutil.copy(f, dest_folder)

shutil.copy('clean_true_state.txt', '../clean_true_state.txt')

shutil.copy('coarse_data.txt', '../coarse_data.txt')
