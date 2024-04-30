"""
@author: samith
"""

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

def get_trajectory_files(root_directory):
    """Retrieve PDB files that contain all frames for trajectory analysis from subdirectories."""
    trajectory_files = []
    for subdir, dirs, files in os.walk(root_directory):
        for file in files:
            if file.endswith('_coordinates_peptideOnly_pbcFixed.pdb'):
                trajectory_files.append(os.path.join(subdir, file))
    return trajectory_files

def analyze_trajectory(trajectory_path):
    """Analyze a PDB trajectory and save Ramachandran plot and angles to CSV."""
    u = mda.Universe(trajectory_path, trajectory_path)  # Load file as both topology and trajectory
    rama = Ramachandran(u.select_atoms('protein'))
    rama.run()

    file_base = os.path.splitext(trajectory_path)[0]
    
    plt.figure(figsize=(8, 6))
    rama.plot()
    plt.title(f"Ramachandran Plot for {os.path.basename(trajectory_path)}")
    plot_path = f"{file_base}_ramachandran.png"
    plt.savefig(plot_path)
    plt.close()

    # Check data dimension and handle accordingly
    if hasattr(rama.results, 'angles') and rama.results.angles.ndim == 3:
        # Taking the mean of angles across all frames for simplicity
        mean_angles = np.mean(rama.results.angles, axis=0)
        phi_angles = mean_angles[:, 0]
        psi_angles = mean_angles[:, 1]
        data = pd.DataFrame({'phi': phi_angles, 'psi': psi_angles})
        csv_path = f"{file_base}_angles.csv"
        data.to_csv(csv_path, index=False)
        return file_base, plot_path, csv_path
    else:
        print(f"Skipping file {trajectory_path}, insufficient or invalid angle data")
        return file_base, plot_path, None

def process_files(files):
    """Process a list of PDB files."""
    results = {}
    for file in files:
        base, plot, csv = analyze_trajectory(file)
        results[base] = {'plot': plot, 'csv': csv}
    return results

if __name__ == "__main__":
    root_directory = '/home/samith/projects/ctb-rmansbac/samith/proj/Benchmark/omegafold/'
    trajectory_files = get_trajectory_files(root_directory)
    results = process_files(trajectory_files)

    for base, files in results.items():
        if files['csv']:
            print(f"Processed {base}: Plot: {files['plot']}, CSV: {files['csv']}")
        else:
            print(f"Skipped processing for {base}, check warnings.")
            
    print("Analysis complete. Results are stored alongside the input files.")

