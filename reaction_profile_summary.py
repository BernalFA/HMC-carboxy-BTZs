#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Short script to retrieve data from calculations in autodE.
It generates a summary csv file.

@author: Dr. Freddy A. Bernal
"""

import os
import sys

FOLDER_MODULE = '/home/tga_user/Documents/Scripts/Reactivity_autodE'

sys.path.append(FOLDER_MODULE)

from autodE_utils import retrieve_data_rxn_profile


if __name__ == '__main__':
    # Get results from all folders in current directory
    summary = retrieve_data_rxn_profile()
    # Define output file
    dirname = os.path.split(os.getcwd())[-1]
    results_file = f"Results_{dirname}.csv"
    # Export dataframe to csv
    summary.to_csv(results_file, index=False)

