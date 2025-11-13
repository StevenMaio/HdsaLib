from mpi4py import MPI
import numpy as np
import exodusii
import subprocess

########## This block contains the variables/function specified by the user on a case-by-case basis #############################

lofi_opt = np.loadtxt('../optimization/final_params_.dat', skiprows=2, max_rows=441)
mask = np.abs(lofi_opt) < 1.e-15

# Specify the variable name assigned to the optimization variable in the exodus file
variable_name = 'thermal_boundary_control'

# Define file names
input_file = 'z_update_mean.exo'  # Base name for the input files
output_file = 'final_params_load_.txt'  # Base name for the output files

##################################################################################################################################

# Open the Exodus file
exo_in = exodusii.File(input_file, mode="r")

# Apply a node-wise transformation (example: square the values)
mean_opt = exo_in.get_node_variable_values(variable_name).squeeze()

mean_opt[mask] = 0.0

np.savetxt(output_file, mean_opt, fmt='%.2f') 
