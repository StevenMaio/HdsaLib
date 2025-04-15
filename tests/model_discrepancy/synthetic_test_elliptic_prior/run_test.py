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

def find_txt_files(directory, exclude_directories):
    """Find all .txt files in the specified directory and its subdirectories, excluding specific directories."""
    txt_files = []
    for root, _, files in os.walk(directory):
        # Skip any directories listed in exclude_directories
        if any(exclude in root for exclude in exclude_directories):
            continue
        for file in files:
            if file.endswith('.txt'):
                txt_files.append(os.path.join(root, file))
    return txt_files

def compare_files(file1, file2):
    """Compare two files and return True if they are the same, False otherwise."""
    return filecmp.cmp(file1, file2, shallow=False)

def delete_txt_files(directory, exclude_directories):
    """Delete all .txt files in the specified directory and its subdirectories, excluding specific directories."""
    for root, _, files in os.walk(directory):
        # Skip any directories listed in exclude_directories
        if any(exclude in root for exclude in exclude_directories):
            continue
        for file in files:
            if file.endswith('.txt'):
                os.remove(os.path.join(root, file))
                #print(f"Deleted: {os.path.join(root, file)}")

def delete_empty_directories(directory, exclude_directories):
    """Delete empty directories in the specified directory, excluding specific directories."""
    for root, dirs, _ in os.walk(directory, topdown=False):
        # Skip any directories listed in exclude_directories
        if any(exclude in root for exclude in exclude_directories):
            continue
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            try:
                os.rmdir(dir_path)
                #print(f"Deleted empty directory: {dir_path}")
            except OSError:
                # Directory is not empty, so we skip it
                continue

def copy_reference_file(reference_file_path, destination_directory):
    """Copy the specified reference file to the destination directory."""
    try:
        shutil.copy(reference_file_path, destination_directory)
        #print(f"Copied {reference_file_path} to {destination_directory}.")
    except Exception as e:
        print(f"Error copying file: {e}")

def main():
    # Define paths
    executable_path = './ROL_example_tests_model_discrepancy_synthetic_test_elliptic_prior_example_sensitivity.exe'
    output_directory = '.'  # Assuming the .exe outputs files in the current directory and subdirectories
    reference_directory = 'reference_output'
    exclude_directories = ['reference_output', 'CMakeFiles']  # Directories to exclude

    # Execute the .exe file
    execute_executable(executable_path)

    # Find all .txt files in the output directory and its subdirectories, excluding specified directories
    output_files = find_txt_files(reference_directory, exclude_directories)

    # Compare each output file with the corresponding reference file
    for output_file in output_files:
        # Create a relative path for the reference file
        relative_path = os.path.relpath(output_file, output_directory)
        reference_file = os.path.join(reference_directory, relative_path)

        if os.path.exists(reference_file):
            if not compare_files(output_file, reference_file):
                print(f"{output_file} does NOT match the reference output.")
        else:
            print(f"Reference file {reference_file} does not exist.")

    # Delete all .txt files generated from the executable
    delete_txt_files(output_directory, exclude_directories)

    # Delete empty directories after removing .txt files
    delete_empty_directories(output_directory, exclude_directories)

    # Copy the reference file to the output directory
    reference_file_to_copy = os.path.join(reference_directory, 'random_numbers.txt')
    copy_reference_file(reference_file_to_copy, output_directory)

if __name__ == "__main__":
    main()
