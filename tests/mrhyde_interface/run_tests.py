import os
import subprocess

# ANSI escape codes for colored output
RED = "\033[91m"
GREEN = "\033[92m"
RESET = "\033[0m"

# Set the file path to the "mrhyde" executable
mrhyde_path = os.path.expanduser('~/software/MrHyDE/build_opt/src/mrhyde')  

# Set the base directory to the current working directory
base_directory = os.getcwd()  

# Function to execute commands and handle errors
def execute_command(command, cwd=None):
    try:
        subprocess.run(command, check=True, shell=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}\n{e}")

# Walk through all subdirectories of the base directory
for subdir, _, _ in os.walk(base_directory):
    # Skip the base directory itself
    if subdir == base_directory:
        continue

    # Print only the folder name
    folder_name = os.path.basename(subdir)
    print(f"Processing: {folder_name}")

    # Create a symbolic link to "mrhyde" in the current subdirectory
    link_command = f"ln -s {mrhyde_path} ."
    execute_command(link_command, cwd=subdir)

    # Execute the mpi command from the current subdirectory
    mpi_command = f"mpiexec -n 2 ./mrhyde > output.log"
    execute_command(mpi_command, cwd=subdir)

    # Compare the log files
    compare_command = f"diff output.log output.gold"
    comparison_result = subprocess.run(compare_command, shell=True, cwd=subdir)

    if comparison_result.returncode != 0:
        print(f"{RED}Test Failed: output.log differs from output.gold{RESET}")
    else:
        print(f"{GREEN}Test Passed: output.log matches output.gold{RESET}")

    # Clean up the output files
    cleanup_command = f"rm -rf hdsa_output output.log mrhyde"
    execute_command(cleanup_command, cwd=subdir)
