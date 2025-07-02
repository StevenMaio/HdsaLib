from mpi4py import MPI
import numpy as np
import exodusii
import subprocess

########## This block contains the variables/function specified by the user on a case-by-case basis #############################

# Specify the variable name assigned to the optimization variable in the exodus file
variable_name = 'thermal_source'

# Define the desired z for high-fidelity simulation input
def Define_z(coords):
    return 1.0

# Set "add_to_lofi_opt" as True to add to the low-fidelity optimization solution to the output of "Define_z", 
# Set "add_to_lofi_opt" as False to use only the output of "Define_z"
add_to_lofi_opt = True

# Define file names
input_file_base = 'output_lofi_opt.exo'  # Base name for the input files
output_file_base = 'hifi_z2.exo'  # Base name for the output files

##################################################################################################################################

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Each process reads its corresponding parallel file
if(size > 1):
    input_file = f"{input_file_base}.{size}.{rank}"
    output_file = f"{output_file_base}.{size}.{rank}"
else:
    input_file = input_file_base
    output_file = output_file_base

# Open the Exodus file
exo_in = exodusii.File(input_file, mode="r")

# Apply a node-wise transformation (example: square the values)
coords = exo_in.get_coords()
transformed_data = Define_z(coords)

if add_to_lofi_opt:
    lofi_opt = exo_in.get_node_variable_values(variable_name).squeeze()
    transformed_data = transformed_data + lofi_opt

# Write the transformed data to a new Exodus file
exo_out = exodusii.File(output_file, mode='w')
exo_out.put_init(
        exo_in.title(),
        exo_in.num_dimensions(),
        exo_in.num_nodes(),
        exo_in.num_elems(),
        exo_in.num_blks(),
        exo_in.num_node_sets(),
        exo_in.num_side_sets(),
        num_edge=exo_in.num_edges(),
        num_edge_blk=exo_in.num_edge_blk(),
        num_face=exo_in.num_faces(),
        num_face_blk=exo_in.num_face_blk(),
    )
exodusii.copy.copy_mesh(exo_in, exo_out)

exo_out.put_node_variable_params(len(exo_in.get_node_variable_names()))
exo_out.put_node_variable_names(exo_in.get_node_variable_names())
exo_out.put_node_variable_values(1,variable_name,transformed_data)
