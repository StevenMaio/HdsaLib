import os
import subprocess
import filecmp
import shutil

def execute_executable(executable_path):
    """Execute the specified executable file."""
    try:
        subprocess.run(executable_path, check=True)
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

def compare_files(file1, file2):
    """Compare two files and return True if they are the same, False otherwise."""
    return filecmp.cmp(file1, file2, shallow=False)

def main():
    # Define paths
    executable_path = './synthetic_test_transient_example_sensitivity.exe'
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
    command = 'rm -rf prior_discrepancy_* posterior_* hessian_evals.txt'

    # Execute the command
    try:
        subprocess.run(command, check=True, shell=True)
        #print("Files removed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error removing files: {e}")

if __name__ == "__main__":
    main()
