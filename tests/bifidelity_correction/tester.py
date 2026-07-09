import os
import subprocess
import filecmp
import shutil
from pathlib import Path
import numpy as np

def execute_executable(executable_path):
    """Execute the specified executable file."""
    try:
        subprocess.run(executable_path, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        #print(f"Executed {executable_path} successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing {executable_path}: {e}")

def find_all_files(directory):
    """Find all files in the specified directory and its subdirectories."""
    all_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            all_files.append(os.path.join(root, file))
    return all_files

def load_vector(filename):
    with open(filename, "r") as f:
        first_line = f.readline().strip()

    if first_line.startswith("%%MatrixMarket"):
        return np.loadtxt(filename, skiprows=2)
    else:
        return np.loadtxt(filename)

def compare_files(file1, file2):
    """Compare two files and return True if they are the same, False otherwise."""
    threshold = 1.e-7
    v1 = load_vector(file1)
    v2 = load_vector(file2)
    if v1.shape != v2.shape:
        raise ValueError(f"Shape mismatch: {v1.shape} vs {v2.shape}")

    denom = np.linalg.norm(v1)
    if denom < 1.e-8:
        rel_diff = np.linalg.norm(v1 - v2)
    else:
        rel_diff = np.linalg.norm(v1 - v2) / denom
    return rel_diff < threshold

def main():
    # Define paths
    exe_files = list(Path(".").glob("*.exe"))
    if len(exe_files) != 1:
        raise ValueError(f"Expected exactly one .exe file, found {len(exe_files)}")
    exe_file = exe_files[0].name
    executable_path = './'+exe_file
    output_directory = '.'  # Assuming the .exe outputs files in the current directory and subdirectories
    reference_directory = 'reference_output'

    # Execute the .exe file
    execute_executable(executable_path)

    # Find all files in the output directory and its subdirectories
    reference_files = find_all_files(reference_directory)

    # Compare each output file with the corresponding reference file
    for reference_file in reference_files:
        # Create a relative path for the reference file
        output_file = reference_file[17:]
        RED = "\033[31m"
        GREEN = "\033[32m"
        if os.path.exists(output_file):
            if not compare_files(output_file, reference_file):
                print(RED + f"{output_file} does NOT match the reference output.")
            else:
                print(GREEN + f"{output_file} matches the reference output.")
        else:
            print(f"Reference file {reference_file} does not exist.")

    # Define the command to execute
    command = 'rm -rf z_update.txt'

    # Execute the command
    try:
        subprocess.run(command, check=True, shell=True)
        #print("Files removed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error removing files: {e}")

if __name__ == "__main__":
    main()
