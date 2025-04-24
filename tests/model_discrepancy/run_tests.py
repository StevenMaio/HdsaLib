import os
import subprocess

def execute_script_in_subdirectories(base_directory, subdirectories):
    """Execute 'python3 run_test.py' in each specified subdirectory."""
    for subdirectory in subdirectories:
        # Construct the full path to the subdirectory
        full_path = os.path.join(base_directory, subdirectory)
        
        # Check if the subdirectory exists
        if os.path.isdir(full_path):
            print(f"Executing in: {full_path}")
            try:
                # Execute the command in the subdirectory
                subprocess.run(['python3', 'run_test.py'], cwd=full_path, check=True)
                #print(f"Successfully executed 'python3 run_test.py' in {full_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error executing 'python3 run_test.py' in {full_path}: {e}")
        else:
            print(f"Subdirectory {full_path} does not exist.")

def main():
    # Define the base directory where the subdirectories are located
    base_directory = '.'  # Change this to the appropriate base directory if needed
    subdirectories = [
        'PDE_Test_Problem',
        'synthetic_test',
        'synthetic_test_elliptic_prior',
        'synthetic_test_hessian_gevp',
        'synthetic_test_hyperparam_1D'
    ]

    # Execute the script in each subdirectory
    execute_script_in_subdirectories(base_directory, subdirectories)

if __name__ == "__main__":
    main()
