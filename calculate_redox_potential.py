#!/usr/bin/env python3
"""
This script will run optimization and frequency calculations for a series
of compounds (SMILES strings). It automatically run the calculation for the 
given structures (defined as the oxidized forms) and for the corresponding 
radical anions (reduced forms). The results are to be used for redox potential
calculations.

Based on the `Reaction` class from autodE.

@author: Dr. Freddy Bernal
"""


# Import essential
import autode as ade
import argparse

from lib.redox import Redox


# Define parser
def arg_parser() -> argparse.Namespace:
    """Parses command-line arguments and returns parsed arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """    
    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('smiles',
                        help='file listing compounds as SMILES and ID (.smi)')
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
                        help='Temperature (default: 298.15 °C)')
    parser.add_argument('-ncores',
                        dest='ncores', 
                        default=14,
                        help='Number of cores (default: 14)')
    parser.add_argument('-j',
                        dest='jobname',
                        help='jobname')
    
    args = parser.parse_args()
    
    return args



def read_file(filename: str) -> list:
    """Reads text file as a list of rows.

    Args:
        filename (str): path to file.

    Returns:
        list: content of file split by line.
    """    

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    lines_out = []
    for line in lines:
        if not (line.startswith('#')) | (line.isspace()):
            lines_out.append(line)
           
    return lines_out


# Helper function to run the simulations for each compound
def calculate(filename: str, name: str, solvent: str, temp: float=298.15):
    """Performs optimization and single point energy calculations for oxidized
    and reduced forms of the given compound. It uses the `Redox` custom class.

    Args:
        filename (str): path to SMI file containing SMILES strings and ID for 
                        the compounds (oxidized, neutral form).
        name (str): jobname or ID.
        solvent (str): solvent name (autodE valid).
        temp (float, optional): temperature (Celsius). Defaults to 298.15.

    """     
    # Define structures to submit for calculation
    compounds = read_file(filename)

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
    # 14 seem to be the max allowed for the available RAM
    ade.Config.n_cores = args.ncores 
    # Add diffuse functions because of anions
    # (ma-def2-SVP is recommendation from autodE troubleshooting)
    ade.Config.ORCA.keywords.set_opt_basis_set(args.basis_set_opt)
    ade.Config.ORCA.keywords.sp.basis_set = args.basis_set_sp
    # Modify input for sp to overcome SCF convergence problems for 
    # radical anions according to ORCA directions
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
