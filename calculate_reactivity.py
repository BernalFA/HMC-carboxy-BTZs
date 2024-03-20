#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script was created to run DFT calculations for a set of reactions of 
interest. A new folder is created to store the results for all the reactions.

It requires autodE (conda env 'autode').
    
Author: Dr. Freddy Bernal
"""


# Import essential
import autode as ade
import pandas as pd
# import numpy as np
import threading
import argparse
import sys
import os


# Define parser
def arg_parser():
    """
    Parses command-line arguments and returns parsed arguments.

    Returns
    -------
    Parsed arguments.

    """
    
    script_usage = """python {} -i infile -j jobname 
    -basis basis_set -func functional -solv solvent -ncores ncores 
    """.format(sys.argv[0][2:])
    
    parser = argparse.ArgumentParser(usage=script_usage, 
                                     description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-i',
                        dest='infile', 
                        help='file with SMILES for reactants and products (TXT, CSV)')
    parser.add_argument('-obasis',
                        dest='basis_set_opt', 
                        default='ma-def2-SVP',
                        help='Basis set for optimizations (default: ma-def2-SVP)')
    parser.add_argument('-sbasis',
                        dest='basis_set_sp', 
                        default='ma-def2-TZVP',
                        help='Basis set for single point calculations (default: ma-def2-TZVP)')
    parser.add_argument('-func',
                        dest='functional', 
                        default=None,
                        help='Functional (default: PBE0-D3BJ)')
    parser.add_argument('-solv',
                        dest='solvent', 
                        default='water',
                        help='Solvent (default: water)')
    parser.add_argument('-ncores',
                        dest='ncores', 
                        default=14,
                        help='Number of cores (default: 14)')
    parser.add_argument('-j',
                        dest='jobname',
                        help='jobname')
    
    args = parser.parse_args()
    
    return args



# Create helper function to run the simulations for each compound
# from pandas row
def run_simulation(row, solvent, jobname): 
    # Define reactants and products
    reac1 = ade.Reactant(name='R1',
                         smiles=row['R1'])
    if not pd.isnull(row['R2']):
        reac2 = ade.Reactant(name='R2', 
                             smiles=row['R2'])
    else:
        reac2 = None
    prod1 = ade.Product(name='P1',
                        smiles=row['P1'])
    if not pd.isnull(row['P2']):
        prod2 = ade.Product(name='P2',
                            smiles=row['P2'])
    else:
        prod2 = None
    
    # Define reaction
    rxn = ade.Reaction(reac1, reac2, prod1, prod2, 
                       name=jobname, 
                       solvent_name=solvent)
    
    # Calculate reaction profile
    rxn.calculate_reaction_profile(free_energy=True)
    
    return None

    

def main():
    # Define arguments
    args = arg_parser()
    
    # define initial path
    path = os.getcwd()
    
    # New folder 
    new_folder = os.path.join(path, args.jobname)
    # Create new folder for job
    if not os.path.exists(new_folder):
        os.mkdir(new_folder)
    
    # Customize autodE
    # Define number of cores to use
    # 16 seem to be the max allowed for the available RAM
    ade.Config.n_cores = args.ncores 
    # Add diffuse functions because of anions
    # (ma-def2-SVP is recommendation from autodE troubleshooting)
    ade.Config.ORCA.keywords.set_opt_basis_set(args.basis_set_opt)
    ade.Config.ORCA.keywords.sp.basis_set = args.basis_set_sp
    
    # Read structures
    if args.infile.endswith('txt'):
        df = pd.read_csv(args.infile, sep='\t')
    elif args.infile.endswith('csv'):
        df = pd.read_csv(args.infile)
    else:
        raise ValueError('File format not supported. Please provide csv or txt file')
    
    # Change dir
    os.chdir(new_folder)
    
    # Iterate over dataframe to run reaction for each compound
    for _, row in df.iterrows():
        name = row["name"]
        print(f'Running reaction for {name}...')
        
        jobname=f'{args.jobname}_{name}'
        # run_simulation(row)
        thread = threading.Thread(target=run_simulation(row, 
                                                        args.solvent, 
                                                        jobname))
        thread.start()
        
        resulting_file = os.path.join(new_folder, 
                                      jobname + '_reaction_profile.pdf')
        while not os.path.exists(resulting_file):
            pass
        
        print(f'Simulation for {name} finished.\n')
        
        

if __name__ == '__main__':
    main()
