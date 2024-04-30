"""
@author: samith
"""

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import matplotlib.pyplot as plt
import pandas as pd
import os
import warnings

# Suppress specific common warnings that are non-critical
warnings.filterwarnings("ignore", category=UserWarning, message="Cannot determine phi and psi angles for the first or last residues")
warnings.filterwarnings("ignore", category=UserWarning, message="Reader has no dt information, set to 1.0 ps")

def get_pdb_files(root_directory):
    """Retrieve PDB files that contain a single frame for analysis from subdirectories."""
    pdb_files = []
    for subdir, _, files in os.walk(root_directory):
        for file in files:
            if file.endswith('_prediction_pdbfixed.pdb'):
                pdb_files.append(os.path.join(subdir, file))
    return pdb_files

def analyze_pdb(pdb_path):
    """Analyze a single frame PDB and save Ramachandran plot and angles to CSV."""
    try:
        u = mda.Universe(pdb_path)
        protein_atoms = u.select_atoms('protein')
        if len(protein_atoms) == 0:
            print(f"No protein atoms found in {pdb_path}. Skipping.")
            return None

        rama = Ramachandran(protein_atoms)
        rama.run()

        if not hasattr(rama.results, 'angles') or not rama.results.angles.any():
            print(f"No angles calculated for {pdb_path}. Skipping.")
            return None

        # Flatten the angles array if it's not 2D
        if rama.results.angles.ndim > 2:
            angles_array = rama.results.angles.reshape(-1, 2)
        else:
            angles_array = rama.results.angles

        # Define base filename for outputs
        file_base = os.path.splitext(pdb_path)[0]
        
        # Generate and save Ramachandran plot
        plt.figure(figsize=(8, 6))
        rama.plot()
        plt.title(f"Ramachandran Plot for {os.path.basename(pdb_path)}")
        plot_path = f"{file_base}_ramachandran.png"
        plt.savefig(plot_path)
        plt.close()

        # Save phi and psi angles to CSV
        data = pd.DataFrame(angles_array, columns=['phi', 'psi'])
        csv_path = f"{file_base}_angles.csv"
        data.to_csv(csv_path, index=False)

        return file_base, plot_path, csv_path

    except Exception as e:
        print(f"An error occurred while processing {pdb_path}: {e}")
        return None

def process_files(files):
    """Process a list of PDB files."""
    results = {}
    for pdb_file in files:
        output = analyze_pdb(pdb_file)
        if output:
            base, plot, csv = output
            results[base] = {'plot': plot, 'csv': csv}
        else:
            results[os.path.splitext(pdb_file)[0]] = {'plot': None, 'csv': None}
    return results

if __name__ == "__main__":
    root_directory = '/home/samith/projects/ctb-rmansbac/samith/proj/Benchmark/rama/omegafold/'  # Adjust this path to your root directory
    pdb_files = get_pdb_files(root_directory)
    results = process_files(pdb_files)
    
    for base, files in results.items():
        if files['plot'] and files['csv']:
            print(f"Processed {base}:")
            print(f"  Plot: {files['plot']}")
            print(f"  CSV: {files['csv']}")
        else:
            print(f"Skipped processing for {base}, check warnings.")
            
    print("Analysis complete. Results are stored alongside the input files.")

