#!/usr/bin/env python3
# =============================================================================
# PLOT COVARIANCE MATRICES FROM FORTRAN PROGRAM RUN_EXAMPLE
# 
# Author: Miroslav HALLO, Kyoto University
# E-mail: hallo.miroslav.2a@kyoto-u.ac.jp
# Revision 2026/04: First version
# Tested with: Python 3.12.3, Matplotlib 3.10.8, NumPy 2.4.4
# Method:
# Hallo, M., Gallovic, F. (2016): Fast and cheap approximation of Green 
#      functions uncertainty for waveform-based earthquake source inversions,
#      Geophys. J. Int., 207 1012-1029. https://doi.org/10.1093/gji/ggw320
# 
# Copyright (C) 2026 Miroslav Hallo
# 
# This program is published under the GNU General Public License (GNU GPL).
# 
# This program is free software: you can modify it and/or redistribute it
# or any derivative version under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
# 
# This code is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
# and don't remove their names from the code.
# 
# You should have received copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.
# 
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# INPUT PARAMETERS
# =============================================================================

# Length of data [samples]
N = 200

# Sampling interval [s]
DT = 0.1

# List of input files with covariance matrices (from Fortran run_example)
INPUTS = ['example_ACF.txt', 'example_AXCF.txt', 
         'example_SACF.txt', 'example_SAXCF.txt']

# List of titles for covariance matrices
TITLES = ['Cov. matrix by ACF', 'Cross-cov. matrix by AXCF', 
         'Cov. matrix by SACF', 'Cross-cov. matrix by SAXCF']
         
# Output filename
OUTPUT = 'example_Fortran.png'
               
# =============================================================================
# FUNCTIONS
# =============================================================================

# Plot covariance matrices
def plot_covariance_matrices(n, dt, files, titles, output='output.png'):
    """
    Plot covariance matrices.

    Args:
        n (int): Length of data [samples].
        dt (float): Sampling interval [s].
        files (list): List of input files with covariance matrices.
        titles (list): List of titles for covariance matrices.
        output (str): Output filename.
    Returns:
        None
    """
    # Prepare time vector
    t = np.arange(n) * dt
    # Start figure
    plt.rcParams.update({
        'font.size': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 16,
    })
    fig, axs = plt.subplots(2, 2, figsize=(12, 10), facecolor='white')
    # Plot subplots
    for i, (ax, file, title) in enumerate(zip(axs.flat, files, titles)):
        try:
            data = np.loadtxt(file)
            im = ax.imshow(data, extent=[t[0], t[-1], t[-1], t[0]], cmap='RdBu_r')
            ax.set_title(title, fontweight='normal')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Time (s)')
            fig.colorbar(im, ax=ax)
            ax.set_aspect('equal')
        except Exception as e:
            print(f"Error during reading {file}: {e}")

    plt.tight_layout()
    
    # Save figure
    plt.savefig(output, dpi=300)
    plt.close()
    print(f"[+] SUCCESS: Saved to '{OUTPUT}'")


# -----------------------------------------------------------------------------
# Main function
def main():
    """
    Main function for standalone execution (plot covariance matrices).

    - Read covariance matrices from INPUTS.
    - Plot covariance matrices.
    - Save results into the OUTPUT.
    """
    print("-" * 50)
    
    print("[*] Plot covariance matrices")
    plot_covariance_matrices(N, DT, INPUTS, TITLES, OUTPUT)


# -----------------------------------------------------------------------------
# Entry point
if __name__ == "__main__":
    main()