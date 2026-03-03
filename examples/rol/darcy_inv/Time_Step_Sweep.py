import xml.etree.ElementTree as ET
import subprocess

# Define the command and arguments
command = ["./darcy_inv_example_inversion.exe"] 
# Open a file to write the output
with open('optimization.txt', 'w') as output_file:
    # Start the process
    process = subprocess.Popen(command, stdout=output_file, stderr=subprocess.STDOUT)
    # Wait for the process to complete
    process.wait()
    # Print a success message
    print("Output written to optimization.txt")

# Define the command and arguments
command = ["./darcy_inv_example_perturbed_inversion.exe"]
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

root.find(".//Parameter[@name='Gradient Tolerance']").set('value','1.e-8') 
root.find(".//Parameter[@name='CG Tolerance']").set('value','1.e-4') 
root.find(".//Parameter[@name='use_qn_prec']").set('value','false') 
root.find(".//Parameter[@name='Use Block Update']").set('value','false') 
tree.write('Sensitivity_input.xml')  

# Loop from 2 to 9
for i in range(2, 10):
    
    root.find(".//Parameter[@name='N_fe']").set('value', str(i)) 
    root.find(".//Parameter[@name='N_me']").set('value', str(i))  
    tree.write('Sensitivity_input.xml')  # Save the changes to the same file
    
    command = ["./darcy_inv_example_sensitivity.exe"]  
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    print("Completed sensitivity computation with "+str(i)+" time steps")

    current_file_name = 'Forward_Euler_Cost_Report.txt'
    new_file_name = 'Forward_Euler_Cost_Report_'+str(i)+'.txt'
    subprocess.run(['mv', current_file_name, new_file_name], check=True)
    current_file_name = 'Modified_Euler_Cost_Report.txt'
    new_file_name = 'Modified_Euler_Cost_Report_'+str(i)+'.txt'
    subprocess.run(['mv', current_file_name, new_file_name], check=True)

root.find(".//Parameter[@name='use_qn_prec']").set('value','true') 
root.find(".//Parameter[@name='Use Block Update']").set('value','true') 
tree.write('Sensitivity_input.xml')  

# Loop from 2 to 9
for i in range(2, 10):
    
    root.find(".//Parameter[@name='N_fe']").set('value', str(i)) 
    root.find(".//Parameter[@name='N_me']").set('value', str(i))  
    tree.write('Sensitivity_input.xml')  # Save the changes to the same file
    
    command = ["./darcy_inv_example_sensitivity.exe"]  
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    print("Completed sensitivity computation with "+str(i)+" time steps")

    current_file_name = 'Forward_Euler_Cost_Report.txt'
    new_file_name = 'Preconditioned_Forward_Euler_Cost_Report_'+str(i)+'.txt'
    subprocess.run(['mv', current_file_name, new_file_name], check=True)
    current_file_name = 'Modified_Euler_Cost_Report.txt'
    new_file_name = 'Preconditioned_Modified_Euler_Cost_Report_'+str(i)+'.txt'
    subprocess.run(['mv', current_file_name, new_file_name], check=True)
