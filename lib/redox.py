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


class Redox:
    """Help to automatically perform optimization and single point energy
    calculations for neutral compounds (oxidized) and their corresponding 
    reduced forms (radical anion), using autodE. Create a dedicated folder
    to stored results, including conformers, optimization, single_points, 
    thermo, and an overall output. 
    """    

    def __init__(
            self, 
            smiles: list, 
            name: str,
            solvent_name: Optional[str] = None, 
            temp: float = 298.15
    ):
        """
        Args:
            smiles (list): SMILES strings for compounds to process.
            name (str): jobname or ID.
            solvent_name (Optional[str], optional): Solvent name (autodE valid). Defaults to None.
            temp (float, optional): Temperature. Defaults to 298.15 °C.
        """        
        self.name = name
        self.oxidized, self.reduced = [], []
        self.solvent = get_solvent(solvent_name, kind="implicit")
        self.temp = float(temp)
        
        if smiles:
            self._init_from_smiles(smiles)
       

    def run_calculation(self) -> None:
        """Perform whole calculation as iterative process."""        

        @work_in(self.name)
        def calculate(redox):
            desc = f"Calculation {self.name} using {type(self).__name__}.\n"
            print(desc + f"{len(self.oxidized)} structures.")
            print("Calculation started.")
            print("CONFORMATIONAL SEARCH")
            redox.find_lowest_energy_conformers()
            print("GEOMETRY OPTIMIZATION")
            redox.optimise_structures()
            print("THERMOCHEMICAL CONTRIBUTIONS")
            redox.calculate_thermochemical_cont()
            redox._solvate()    
            print("SINGLE POINT CALCULATIONS")
            redox.calculate_single_points()
            redox.print_output()
            return None
        
        calculate(self)

        return None
    

    def _init_from_smiles(self, smiles: list) -> None:
        """Initialize calculation from SMILES string

        Args:
            smiles (list): list of SMILES

        """        
        # Add all the oxidized and reduced species
        for i, smi in enumerate(smiles):
            ox = Molecule(smiles=smi.split()[0])
            ox.name = f"ox{i}_{smi.split()[1]}"
            self.oxidized.append(ox)

        return None


    def find_lowest_energy_conformers(self) -> None:
        """Perform conformational search for each compound at DFT level."""        
        
        print("Conformational search started.")
        h_method = get_hmethod() if Config.hmethod_conformers else None
        for mol in self.oxidized:
            print("  " + mol.name)
            mol.find_lowest_energy_conformer(hmethod=h_method)
        print("Conformational search completed.")

        return None
    

    @work_in("optimization")
    def optimise_structures(self) -> None:
        """Perform geometry optimization for every compound in oxidized 
        and reduced form, using the indicated DFT method"""

        print("Optimizations started.")
        h_method = get_hmethod()

        for mol in self.oxidized:
            print("  " + mol.name)
            mol.optimise(h_method)
            red = mol.copy()
            red.charge = -1 
            red.mult = 2
            red.name = red.name.replace("ox", "red")
            self.reduced.append(red)

        for mol in self.reduced:
            print("  " + mol.name)
            mol.optimise(h_method)

        print("Optimizations completed.")

        return None

    
    @work_in("thermo")
    def calculate_thermochemical_cont(self, free_energy=True, enthalpy=True):
        """Perform calculation of thermochemical contributions for each 
        compound in both forms.

        Args:
            free_energy (bool, optional): if True, calculate free energy 
                                          component. Defaults to True.
            enthalpy (bool, optional): if True, calculate enthalpy. 
                                       Defaults to True.

        """        
        
        if not (free_energy or enthalpy):
            return None
        print("Thermochemical contributions started.")
        for mol in self.oxidized + self.reduced:
            print("  " + mol.name)
            mol.calc_thermo(temp=self.temp)
            
        print("Thermochemical contributions completed.")

        return None


    def _solvate(self) -> None:
        """add solvent to each species (if provided)."""

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
        """Perform single point energy calculation for each compound
        in both forms.
        """

        print("Single point calculations started.")
        h_method = get_hmethod()
        for mol in self.oxidized + self.reduced:
            try:
                print("  " + mol.name)
                mol.single_point(h_method)
            except CouldNotGetProperty:
                print(f"Could not get energy for {mol.name}")
            
        print("Single point calculations completed.")

        return None
        
    
    @work_in("output")
    def print_output(self):
        """Help storing XYZ files of each optimized structure. Create 
        a CSV file with energies and another with imaginary frequency checks. 
        """

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
            for n, f in zip(names, frequencies):
                imfreq.write(f"{n}: {f}\n")
    
        return None
