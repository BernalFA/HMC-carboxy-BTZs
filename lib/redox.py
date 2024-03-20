"""
Class to run Standard Redox potential calculations using 
ORCA and autodE.
Based on the `Reaction` class from autodE. 

Author: Dr. Freddy Bernal

"""

from typing import Optional
from autode.methods import get_hmethod
from autode.config import Config
from autode.utils import work_in
from autode.species.molecule import Molecule
from autode.solvent.solvents import get_solvent
from autode.exceptions import CouldNotGetProperty
from datetime import date
import tqdm


BAR_FORMAT = '{l_bar}{bar:20}{r_bar}'


class Redox:
    def __init__(
            self, 
            smiles: list, 
            name: str,
            solvent_name: Optional[str] = None, 
            temp: float = 298.15
    ):
        
        self.name = name
        self.oxidized, self.reduced = [], []
        self.solvent = get_solvent(solvent_name, kind="implicit")
        self.temp = float(temp)
        
        if smiles:
            self._init_from_smiles(smiles)
       

    def run_calculation(self):

        @work_in(self.name)
        def calculate(redox):
            desc = f"Calculation {self.name} using {type(self).__name__}.\n"
            print(desc + f"{len(self.oxidized)} structures.")
            print("Calculation started.")
            redox.find_lowest_energy_conformers()
            redox.optimise_structures()
            redox.calculate_thermochemical_cont()
            redox._solvate()    
            redox.calculate_single_points()
            redox.print_output()
            return None
        
        calculate(self)

        return None
    

    def _init_from_smiles(self, smiles) -> None:
        """
        Initialise from a SMILES string 
        -----------------------------------------------------------------------
        Arguments:
            list of smiles with names (str):
        """
        # Add all the oxidized and reduced species
        for i, smi in enumerate(smiles):
            ox = Molecule(smiles=smi.split()[0])
            ox.name = f"ox{i}_{smi.split()[1]}"
            self.oxidized.append(ox)

        return None


    def find_lowest_energy_conformers(self) -> None:
        print("Conformational search started.")
        h_method = get_hmethod() if Config.hmethod_conformers else None
        for mol in tqdm.tqdm(self.oxidized, 
                             desc="ConfSearch",
                             bar_format=BAR_FORMAT):
            mol.find_lowest_energy_conformer(hmethod=h_method)
        print("Conformational search completed.")

        return None
    

    @work_in("optimization")
    def optimise_structures(self) -> None:
        """Perform a geometry optimisation on all the reactants and products
        using the method"""

        print("Optimizations started.")
        h_method = get_hmethod()

        for mol in tqdm.tqdm(self.oxidized, 
                             desc="OptOx",
                             bar_format=BAR_FORMAT):
            mol.optimise(h_method)
            red = mol.copy()
            red.charge = -1 
            red.mult = 2
            red.name = red.name.replace("ox", "red")
            self.reduced.append(red)

        for mol in tqdm.tqdm(self.reduced, 
                             desc="OptRed",
                             bar_format=BAR_FORMAT):
            mol.optimise(h_method)

        print("Optimizations completed.")

        return None

    
    @work_in("thermo")
    def calculate_thermochemical_cont(self, free_energy=True, enthalpy=True):
        
        if not (free_energy or enthalpy):
            return None
        print("Thermochemical contributions started.")
        for mol in tqdm.tqdm(self.oxidized + self.reduced, 
                             desc="g_count",
                             bar_format=BAR_FORMAT):
            mol.calc_thermo(temp=self.temp)
            
        print("Thermochemical contributions completed.")

        return None


    def _solvate(self) -> None:

        if self.solvent is not None:
            for mol in self.oxidized + self.reduced:
                mol_solv = mol.copy()
                mol_solv.name = mol_solv.name + "_solv"
                mol_solv.solvent = self.solvent
                if 'ox' in mol_solv.name:
                    self.oxidized.append(mol_solv)
                elif 'red' in mol_solv.name:
                    self.reduced.append(mol_solv)
        
        print(f"Structures solvated ({self.solvent}).")

        return None
    

    @work_in("single_points")
    def calculate_single_points(self) -> None:
        print("Single point calculations started.")
        h_method = get_hmethod()
        for mol in tqdm.tqdm(self.oxidized + self.reduced, 
                             desc="SP",
                             bar_format=BAR_FORMAT):
            try:
                mol.single_point(h_method)
            except CouldNotGetProperty:
                print(f"Could not get energy for {mol.name}")
            
        print("Single point calculations completed.")

        return None
        
    
    @work_in("output")
    def print_output(self):
        print("Printing results.")
        from autode.log.methods import methods
        
        with open("methods.txt", "w") as out_file:
            print(methods, file=out_file)
            
        csv_file = open("energies.csv", "w")
        method = get_hmethod()
        print(
            f"Energies generated by autodE on: {date.today()}. Single point "
            f"energies at {method.keywords.sp.bstring} and optimisations at "
            f"{method.keywords.opt.bstring}",
            "Species,E_opt,G_cont,H_cont,E_sp",
            sep="\n",
            file=csv_file,
        )
        
        def print_energies_to_csv(_mol):
            print(
                f"{_mol.name}",
                f"{_mol.energies.first_potential}",
                f"{_mol.g_cont}",
                f"{_mol.h_cont}",
                f"{_mol.energies.last_potential}",
                sep=",",
                file=csv_file,
            )

        frequencies = []
        names = []            
        for mol in self.oxidized + self.reduced:
            mol.print_xyz_file()
            print_energies_to_csv(mol)
            if "_solv" not in mol.name:
                names.append(mol.name)
                frequencies.append(mol.imaginary_frequencies)
            
        with open("imaginary_frequencies.csv", "w") as imfreq:
            # for n, f in zip(names, frequencies):
            #     imfreq.write(f"{n}: {f}")
                imfreq.writelines([f"{n}: {freq}" for n, freq in zip(names, frequencies)])
    
        return None
