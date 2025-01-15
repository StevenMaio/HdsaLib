import xml.etree.ElementTree as ET
import subprocess

# Define the command and arguments
command = ["./ROL_example_rol_darcy_inv_example_inversion.exe"] 
# Open a file to write the output
with open('optimization.txt', 'w') as output_file:
    # Start the process
    process = subprocess.Popen(command, stdout=output_file, stderr=subprocess.STDOUT)
    # Wait for the process to complete
    process.wait()
    # Print a success message
    print("Output written to optimization.txt")

# Define the command and arguments
command = ["./ROL_example_rol_darcy_inv_example_perturbed_inversion.exe"]
# Open a file to write the output
with open('perturbed_optimization.txt', 'w') as output_file:
    # Start the process
    process = subprocess.Popen(command, stdout=output_file, stderr=subprocess.STDOUT)
    # Wait for the process to complete
    process.wait()
    # Print a success message
    print("Output written to perturbed_optimization.txt")

# Load the XML file
tree = ET.parse('Sensitivity_input.xml')  
root = tree.getroot()

# Loop from 2 to 9
for i in range(0, 25, 5):
    for j in range(0, 25, 5):
        rank_param = root.find(".//Parameter[@name='rank']")  # Use XPath to find the element
        if rank_param is not None:
            rank_param.set('value', str(i))  # Update the value attribute with the current loop index
        update_param = root.find(".//Parameter[@name='Maximum Block Update Storage']")  # Use XPath to find the element
        if update_param is not None:
            update_param.set('value', str(j))  # Update the value attribute with the current loop index
        # Save the changes back to the XML file
        tree.write('Sensitivity_input.xml')  # Save the changes to the same file
        
        command = ["./ROL_example_rol_darcy_inv_example_sensitivity.exe"]  
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print("Completed sensitivity computation with initial rank "+str(i)+" and update rank "+str(j))

        current_file_name = 'Modified_Euler_Cost_Report.txt'
        new_file_name = 'Modified_Euler_Cost_Report_'+str(i)+'_'+str(j)+'.txt'
        subprocess.run(['mv', current_file_name, new_file_name], check=True)

