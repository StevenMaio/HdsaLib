import numpy as np
import exodusii
import pdb 
import matplotlib.pyplot as plt
import matplotlib.animation as animation

########## This block contains the variables/function specified by the user on a case-by-case basis #############################

# Specify the variable name assigned to the optimization variable in the exodus file
variable_name = 'thermal_source'


# Define file names
input_file = 'output_after_optimization.exo'  # Base name for the input files

# Open the Exodus file
exo_in = exodusii.File(input_file, mode="r")

# Apply a node-wise transformation (example: square the values)
x = exo_in.get_coords()
c1 = exo_in.get_node_variable_values('c').squeeze()
c2 = exo_in.get_node_variable_values('c2').squeeze()
source2 = exo_in.get_node_variable_values('source2').squeeze()

m = c1.shape[0]

# Create a figure and axis
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

# Initialize the lines for both c1 and c2
line1, = ax1.plot(x, c1[0, :], color='blue', label='c1')
line2, = ax2.plot(x, c2[0, :], color='orange', label='c2')

# Set the limits and labels for c1
ax1.set_xlim(x.min(), x.max())
ax1.set_ylim(c1.min(), c1.max())
ax1.set_ylabel('Field Value (c1)')
ax1.set_title('Time Evolution of c1')
ax1.legend()

# Set the limits and labels for c2
ax2.set_xlim(x.min(), x.max())
ax2.set_ylim(c2.min(), c2.max())
ax2.set_xlabel('Spatial Coordinate (x)')
ax2.set_ylabel('Field Value (c2)')
ax2.set_title('Time Evolution of c2')
ax2.legend()

# Animation update function
def update(frame):
    line1.set_ydata(c1[frame, :])  # Update the line for c1
    line2.set_ydata(c2[frame, :])  # Update the line for c2
    ax1.set_title(f'Time Evolution of c1')  # Update title for c1
    ax2.set_title(f'Time Evolution of c2')  # Update title for c2
    return line1, line2

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=m, blit=True, repeat=True)

# Show the animation
plt.tight_layout()
plt.show()