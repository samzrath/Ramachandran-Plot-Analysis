import os

# Set the root directory for your project
root_dir = '/home/samith/projects/ctb-rmansbac/samith/proj/Benchmark/rama/experimental/'

# Function to rename files
def rename_files():
    for subdir, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith('_pdbfixed_angles.csv'):
                old_file_path = os.path.join(subdir, file)
                # Creating new file name
                base_name = file.replace('_pdbfixed_angles.csv', '')
                new_file_name = f'{base_name}_experimental_prediction_pdbfixed_angles.csv'
                new_file_path = os.path.join(subdir, new_file_name)
                # Rename the file
                os.rename(old_file_path, new_file_path)
                print(f'Renamed: {old_file_path} to {new_file_path}')

# Call the function
rename_files()
