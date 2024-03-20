"""
Module containing useful functions for data processing and analysis.
It includes combination of calculated and experimental data, calculation
of correlation statistics, and retrieval of reaction energy data calculated 
by autodE. 

@author: Dr. Freddy Bernal
"""

import os
import numpy as np
import pandas as pd

from scipy.stats import pearsonr, spearmanr, kendalltau
from scipy.stats import shapiro

ha2kcalmol = 627.509 
ha2Volt = 2625000.5 / 96500


def combine_data(prop: pd.DataFrame, exp: pd.DataFrame) -> pd.DataFrame:
    """Create a concatenated dataframe between predicted and experimental data.

    Args:
        prop (pd.DataFrame): calculated properties/descriptors.
        exp (pd.DataFrame): experimentally determined property.

    Returns:
        pd.DataFrame: concatenated dataframe.
    """    
    # Define a filtering mask
    mask = exp.Name.isin(prop.Name.unique().tolist())
    # Filter exp data frame
    cols = [x for x in exp.columns if ('M.' in x) or ('MH' in x)]
    exp_red = exp[mask].loc[:, ['Name'] + cols].copy()
    # Rename columns
    exp_red.columns = ['Name', 'MIC_Mv', 'MIC_Ms', 'MHC']
    exp_red['pMIC_Mv'] = -np.log10(exp_red.MIC_Mv / 1e6)
    exp_red['pMHC'] = -np.log10(exp_red.MHC)

    # Now sort exp_red such same order as all_atomic_neutral (sequence of 3 steps)
    exp_red.set_index('Name', inplace=True)
    exp_red = exp_red.reindex(index=prop.Name)
    exp_red.reset_index(inplace=True)
    
    # Combine the two dataframes
    out = pd.concat((exp_red, prop.reset_index().iloc[:,2:]), axis=1)
        
    return out 


def format_compound_name(name: str) -> str:
    """transform compound names by removing unnecessary suffixes 
    (linked by hyphen).
    """    
    name_list = name.split('-')
    if len(name_list) >= 2:
        newname = ''.join(x for x in name_list)
    elif len(name_list) == 1:
        newname = name_list[0]
    return newname


def correlation_statistics(df: pd.DataFrame, 
                           variables: list[str], 
                           reference: str) -> None:
    """Calculate Pearson, Spearman, and Kendall correlations for the listed 
    variables in a dataframe (against a reference variable). It first check
    for normality

    Args:
        df (pd.DataFrame): dataframe.   
        variables (list[str]): list of columns to compare with reference.
        reference (str): reference column for calculation.
    """    
    labels = ['Pearson correlation', 'Spearman correlation', 'Kendall correlation']

    # Check normality
    print('Shapiro Test:')
    for var in variables + [reference]:
        shap = shapiro(df[var])
        print(f'  {var:20}: {shap.pvalue:.2}')
    print('\n')

    # Calculate correlations (individually)
    for var in variables:
        pearson = pearsonr(df[reference], df[var])
        spearman = spearmanr(df[reference], df[var])
        kendall = kendalltau(df[reference], df[var])
        print(var + ':')
        for lab, test in zip(labels, [pearson, spearman, kendall]):
            try:
                print(f'  {lab:20}: {test.statistic:.2} ({test.pvalue:.2})')
            except AttributeError:
                print(f'  {lab:20}: {test.correlation:.2} ({test.pvalue:.2})')


############################
# autodE related functions #
############################

def delta(df: pd.DataFrame, units:str='Ha') -> np.array:
    """Calculate reaction energy and activation energy after automatic autodE
    reaction profile calculation.

    Args:
        df (pd.DataFrame): dataframe with energy values from autodE (energies.csv). 
        units (str, optional): Energy units ('Ha', 'kcal mol-1'). Defaults to 'Ha'.

    Returns:
        np.array: array containing reaction energy and activation energy values. 
                  The values follow the order dE, dEts, dG, dGts 
                  (ts=Transition State).
    """

    # Define free energy
    df['G'] = df[' E_sp'] + df[' G_cont']
    # Separate energy and free energy for reactants, products, and TS
    reacts = df[df.Species.str.startswith("R")]
    prods = df[df.Species.str.startswith("P")]
    ts = df[df.Species.str.startswith("TS")]

    # Calculate deltas for reaction and barrier
    dE = prods[' E_sp'].sum() - reacts[' E_sp'].sum()
    dG = prods['G'].sum() - reacts['G'].sum()
    # Verify TS existence
    if len(ts) != 0:
        dEts = ts[' E_sp'].values[0] - reacts[' E_sp'].sum()
        dGts = ts['G'].values[0] - reacts['G'].sum()
    else:
        if dG <= 0:
            dEts = 0.00694 # Effective free energy barrier of 4.35 kcal/mol (autodE)
            dGts = 0.00694 # Effective free energy barrier of 4.35 kcal/mol (autodE)
        else:
            dEts = dE + 0.00694 
            dGts = dG + 0.00694 

    
    res = [dE, dEts, dG, dGts]
    
    if units == 'Ha':
        return np.asarray(res)
    
    elif units == 'kcal mol-1':
        res = [val * ha2kcalmol for val in res]
        
        return np.asarray(res)
    


def get_name(folder):
    name = folder.split('_')[-1]
    return name


def retrieve_reaction_profile() -> pd.DataFrame:
    """Search for results from autodE run in the current directory and compile
    the information in a single dataframe.

    Returns:
        pd.DataFrame: reaction and activation energy values in kcal/mol for
                      all the reaction calculations found.
    """    
    ref_path = os.getcwd()
    # Define list of folders
    folders = [folder for folder in os.listdir(ref_path)]
    # Create empty dict to store results
    results = {}
    # Iterate over folders
    for folder in folders:
        # Create a path to look for
        path = os.path.join(ref_path, folder)
        # Verify it is dir
        if os.path.isdir(path):
            # Define file with energy values
            file = os.path.join(ref_path, folder, 'output/energies.csv')
            # read file with Pandas
            try:
                df = pd.read_csv(file, skiprows=1)
            except FileNotFoundError:
                print(f'No output file in {folder}')
                continue
            # Get delta values
            res = delta(df, units='kcal mol-1')
            name = get_name(folder)
            # Store results in dict
            results[name] = res
    # Create empty dataframe            
    summary = pd.DataFrame(columns=['Cmpd', 'dE', 'dE‡', 'dG', 'dG‡'])
    # Convert dict to dataframe
    for key, val in results.items():
        summary.loc[len(summary)] = [key] + val.tolist()
    
    return summary


def reformat_species_col(col):
    """Conveniently separate species/state information contained on 
    given string.
    """    
    parts = col.split('_')
    if len(parts) == 2:
        _, name = parts
        s = 'gas'
    elif len(parts) == 3:
        _, name, s = parts
    if col.startswith('ox'):
        b = 'ox'
    elif col.startswith('red'):
        b = 'red'
    return [b, name, s]


def calculate_redox_potential(filename: str) -> pd.DataFrame:
    """Perform redox potential calculation for a series of results from 
    calculate_redox_potential.py script.

    Args:
        filename (str): path to file with autodE energies (energies.csv)

    Returns:
        pd.DataFrame: Gibbs free energy (kcal/mol) and redox potential (V) 
                      values for each compound.
    """

    # Load file as DataFrame
    df = pd.read_csv(filename, skiprows=1)
    # Sort df by species name
    df.sort_values(by='Species', inplace=True, ignore_index=True)    
    # Convert species to name, type, and id
    df['Species'] = df['Species'].apply(lambda x: reformat_species_col(x))
    df[['species', 'name', 'state']] = df['Species'].tolist()
    df.drop(columns='Species', inplace=True)
    # Check E_sp was really calculated (if ORCA failed, second E = first E)
    df['diff'] = df['E_opt'] - df['E_sp']
    exclude = df.query('diff == 0')['name'].unique()
    newdf = df.query('name not in @exclude').copy()
    # Group by compound and calculate dGs
    result = pd.DataFrame(columns=['Name', 'dG_aq (Ha)', 'Ered (V)'])
    for name, group in newdf.groupby(by='name'):
        # define properties by species and state
        ox_g = group.query('state == "gas" & species == "ox"').copy()
        red_g = group.query('state == "gas" & species == "red"').copy()
        ox_s = group.query('state == "solv" & species == "ox"').copy()
        red_s = group.query('state == "solv" & species == "red"').copy()
        for tmp in [ox_g, red_g, ox_s, red_s]:
            # reindex to allow operations
            tmp.reset_index(inplace=True, drop=True)
        # Calculate results for those compounds with full ORCA runs
        dG_rxn_gas = (red_g['E_sp'] + red_g['G_cont']) - (ox_g['E_sp'] + ox_g['G_cont'])
        dG_solv_ox = ox_s['E_sp'] - ox_g['E_sp']
        dG_solv_red = red_s['E_sp'] - red_g['E_sp']
        dG_rxn_sln = dG_rxn_gas + dG_solv_red - dG_solv_ox
        Ered = -dG_rxn_sln * ha2Volt
        result.loc[len(result)] = [name, dG_rxn_sln.values[0], Ered.values[0]]
            
    return result
