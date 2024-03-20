#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script will run optimization and frequency calculations for a series
of compounds (in SMILES strings). Based on the Reaction class from autodE.

It requires autodE (conda env 'autode').
    
Author: Dr. Freddy Bernal
"""


# Import essential
import autode as ade
import argparse
import sys

sys.path.append("/home/tga_user/Documents/Scripts/Reactivity_autodE/")

from autode_module import Redox


# Define parser
def arg_parser():
    """
    Parses command-line arguments and returns parsed arguments.

    Returns
    -------
    Parsed arguments.

    """
    
    script_usage = """python {} smiles -j jobname 
    -basis basis_set -func functional -solv solvent -ncores ncores 
    """.format(sys.argv[0])
    
    parser = argparse.ArgumentParser(usage=script_usage, 
                                     description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('smiles',
                        help='file with list of SMILES with compound ID (.smi)')
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
    parser.add_argument('-temp',
                        dest='temp',
                        default=298.15,
                        help='Temperature (default: 298.15 K)')
    parser.add_argument('-ncores',
                        dest='ncores', 
                        default=14,
                        help='Number of cores (default: 14)')
    parser.add_argument('-j',
                        dest='jobname',
                        help='jobname')
    
    args = parser.parse_args()
    
    return args



def read_file(filename):
    """
    Read contents of the specified file.

    Parameters
    ----------
    filename : str
        Name of the file to be read.

    Returns
    -------
    lines_out : list of str
        Contents of the file split by line.

    """

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    lines_out = []
    for line in lines:
        if not (line.startswith('#')) | (line.isspace()):
            lines_out.append(line)
           
    return lines_out


# Create helper function to run the simulations for each compound
# from pandas row
def calculate(smiles, name, solvent, temp=298.15): 
    # Define structures to submit for calculation
    compounds = read_file(smiles)

    # Define calculation
    calc = Redox(compounds,
                 name=name,
                 solvent_name=solvent,
                 temp=temp)
    
    # Calculate reaction profile
    calc.run_calculation()
    
    return None

        

if __name__ == '__main__':
    # Define arguments
    args = arg_parser()
    
    # Customize autodE
    # Define number of cores to use
    # 16 seem to be the max allowed for the available RAM
    ade.Config.n_cores = args.ncores 
    # Add diffuse functions because of anions
    # (ma-def2-SVP is recommendation from autodE troubleshooting)
    ade.Config.ORCA.keywords.set_opt_basis_set(args.basis_set_opt)
    ade.Config.ORCA.keywords.sp.basis_set = args.basis_set_sp
    # Modify input for sp to overcome SCF convergence problems for radical anions
    # according to ORCA directions
    ade.Config.ORCA.keywords.sp.append(
        "\n"
        "%scf\n"
        "soscfmaxit 12\n"
        "directresetfreq 1\n"
        "maxiter 250\n"
        "end\n"
    )

    # run_simulation
    if args.temp != float(298.15):
        calculate(args.smiles, args.jobname,
                  args.solvent, float(args.temp))
    else: 
        calculate(args.smiles, args.jobname,
                  args.solvent)
