#!/usr/bin/env python3
"""
Short script to retrieve data from calculations in autodE.
It generates a CSV file summarizing results. 
All autodE generated folders muts be on the current directory.

@author: Dr. Freddy A. Bernal
"""

import os

from lib.utils import retrieve_reaction_profile


if __name__ == '__main__':
    # Get results from all folders in current directory
    summary = retrieve_reaction_profile()
    # Define output file
    dirname = os.path.split(os.getcwd())[-1]
    results_file = f"Results_{dirname}.csv"
    # Export dataframe to csv
    summary.to_csv(results_file, index=False)

