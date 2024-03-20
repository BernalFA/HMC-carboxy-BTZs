"""
Module containing a useful class to preprocess QM descriptors calculated 
with Jaguar.

@author: Dr. Freddy Bernal
"""

import ast
import pandas as pd


class QMDescriptors:
    """Class to process CDFT descriptors from Jaguar calculations."""    

    def __init__(self, file_properties: str, file_atoms: str):
        """Basic init

        Args:
            file_properties (str): path to CSV file with QM descriptors 
                                   (e.g. obtained with qm_descriptors.py 
                                   from Schrödinger)
            file_atoms (str): path to TXT file with atom numbers (core under 
                              study) per compound
        """        
        self.file_props = file_properties
        self.file_atoms = file_atoms
        self._glob_props = [
            'HOMO_', 
            'LUMO_', 
            'Mean_ESP', 
            'Dipole_Moment'
        ]
        self._prop_str2name = {
            'NN_HOMO' : 'fhomo',
            'NN_LUMO' : 'flumo', 
            'Mulliken': 'mul',
            'from_ESP': 'cesp', 
            'Lowdin'  : 'lowdin',
            'Potential_at': 'espn'
        }


    def get_global_descriptors(self) -> pd.DataFrame:
        """Calculate chemical potential (u) and electrophilicity index (w) as 
        global (molecular) CDFT descriptors, as follows:

            u = (Eh + El) / 2
            w = u**2 / (2 * n)

        where n is the Hardness (represented by the HOMO-LUMO gap), El and Eh 
        are LUMO and HOMO energies, respectively.

        Returns:
            pd.DataFrame: molecular CDFT descriptors
        """        

        self._get_data()
        
        # First take only columns with HOMO, LUMO, gap and charge 
        sel = []
        for col in self._rawdata.columns:
            # if ('HOMO_' in col) | ('LUMO_' in col) | ('Mol_Ch' in col):
            if any(x in col for x in self._glob_props):
                sel.append(col)
                
        # Create a new reduced dataframe with those columns
        desc = self._rawdata[['Name'] + sel].copy()
        
        # Define columns for calculation
        for col in sel:
            if '_LUMO_' in col:
                El = col
            elif '_HOMO_' in col:
                Eh = col
            elif 'Gap' in col:
                gap = col
        # Calculate descriptors and add them to df
        desc['Chemical_potential'] = (desc[El] + desc[Eh]) / 2
        desc['Electrophilicity_index'] = desc['Chemical_potential']**2 / (2 * desc[gap])
        desc.set_index("Name", inplace=True, drop=True)

        return desc


    def get_local_descriptors(self, core: str = 'nitrobtz') -> pd.DataFrame:
        """Retrieve and organize local (atomic) QM descriptors for a specified
        core/scaffold.

        Args:
            core (str): name of the core (benzene, btz, nitrobtz). The number
                        of atoms in the selected core must match the number 
                        of atoms given in the file_atoms parameter.
                        Defaults to 'nitrobtz'.

        Returns:
            pd.DataFrame: atomic QM descriptors on selected core.
        """        
        self._get_data()

        data = self.retrieve_atomic_descriptors()

        cols_explode = [x for x in data.columns if x != "Name"]
        data = data.explode(cols_explode)
        data["Atom_Number"] = data["Atom_Name"].apply(self._atoms2num)

        indexes = self._read_atom_idx()
        # check consistency 
        self._check_index_consistency(indexes, core)
        
        results = []
        for name in data.Name.unique():
            if name not in indexes.keys():
                print(f"Warning: {name} not found in {self.file_atoms}")
                continue
            
            idxs = indexes[name]
            cmpd = data[data.Name == name]
            data4atoms = self._sel_core_atomic_desc(cmpd, idxs)

            data4atoms["Atom"] = self._colname_by_core(core=core)

            results.append(data4atoms)

        results = pd.concat(results)

        cols_pivot = [x for x in results.columns if x not in ["Name", "Atom"]]
        descriptors = results.pivot_table(
            index="Name",
            columns="Atom",
            values=cols_pivot
        )
        descriptors.columns = ["_".join(col) for col in descriptors.columns.values]
        return descriptors


    def normalized_fukui(self, core: str = 'nitrobtz', orb: str = 'lumo') -> pd.DataFrame:
        """Normalize Fukui indexes over the whole molecule and return values
        for specified core/scaffold atoms.

        Args:
            core (str, optional): name of core (benzene, btz, nitrobtz). The 
                                  number of atoms in the selected core must 
                                  match those in the given file_atoms. 
                                  Defaults to 'nitrobtz'.
            orb (str, optional): orbital of interest (lumo, homo). 
                                 Defaults to 'lumo'.

        Returns:
            pd.DataFrame: normalized Fukui indexes.
        """        
        self._get_data()

        data = self.retrieve_atomic_descriptors()
        # # Add compound charge
        # results = results.join(self._rawdata.i_j_Mol_Charge)
        fukui = data[['f' + orb, 'Atom_Name', 'Name']]

        cols_explode = [x for x in fukui.columns if x != "Name"]
        fukui = fukui.explode(cols_explode)
        fukui["Atom_Number"] = fukui["Atom_Name"].apply(self._atoms2num)

        indexes = self._read_atom_idx()
        # check consistency 
        self._check_index_consistency(indexes, core)
        
        results = []
        for name, group in fukui.groupby('Name'):
            if name not in indexes.keys():
                print(f"Warning: {name} not found ")
                continue

            idxs = indexes[name]
            
            norm_data = group.flumo.sub(group.flumo.min()).div(group.flumo.max() - group.flumo.min())

            norm_df = pd.concat((group[['Atom_Number']], norm_data), axis=1)

            data4atoms = self._sel_core_atomic_desc(norm_df, idxs)
            data4atoms.columns = [name]
            data4atoms["Atom"] = self._colname_by_core(core=core)
            data4atoms.set_index('Atom', inplace=True)

            results.append(data4atoms.T)

        results = pd.concat(results)
        return results


    def retrieve_atomic_descriptors(self) -> pd.DataFrame:
        """Extract atomic properties according to _prop_str2name attribute.

        Returns:
            pd.DataFrame: full set of atomic descriptors.
        """        
        # atomic properties are stored as string representing a list
        # they need to be converted to individual objects 
        results = {}
        for string in self._prop_str2name.keys():
            col = self._column_selector(self._rawdata.columns, string)
            if col is not None:
                # numeric data can be obtained from literal_eval
                data = self._rawdata[col].apply(lambda x: ast.literal_eval(x))
                newcol = self._prop_str2name[string]
                results[newcol] = data
            else:
                continue
            # print(f"{string}\n")

        atom_col = self._column_selector(self._rawdata.columns, "_Name")
        atoms = self._rawdata[atom_col].apply(self._str2list)
        results["Atom_Name"] = atoms

        results = pd.DataFrame(results)
        # Add compound names
        results = results.join(self._rawdata.Name)
        # Add compound charge
        # results = results.join(self._rawdata.i_j_Mol_Charge)
        return results


    def _atoms2num(self, atom : str) -> int:
        """Transform atom numbering to integer (e.g. C12 -> 12).

        Args:
            atom (str): atom numbering given in Maestro.

        Returns:
            int: number of the specified atom.
        """        
        try:
            num = int(atom[1:])
        except ValueError:
            num = int(atom[2:])
        return num


    def _check_index_consistency(self, indexes: dict, core: str) -> None:
        """Verify specified core atoms match given atoms.

        Args:
            indexes (dict): atom indexes per compound from _read_atom_idx.
            core (str): name of core (benzene, btz, nitrobtz).

        Raises:
            ValueError: in case there is no number match between core atoms 
                        and given atoms
        """        
        index_len = {'nitrobtz': 13, 'btz': 10, 'benzene': 6}
        for key, val in indexes.items():
            if len(val) != index_len[core]:
                raise ValueError(f"Length of indexes for {key} does not match {core=}")


    def _colname_by_core(self, core: str) -> list[str]:
        """Create column names based on selected core.

        Args:
            core (str): name of core (benzene, btz, nitrobtz).

        Raises:
            NotImplementedError: in case core name does not match those implemented.

        Returns:
            list[str]: column names representing the core atoms.
        """        
        if core == 'nitrobtz':
            newcols = ['C4', 'C4a', 'C5', 'C6', 'C7', 
                       'C8', 'N', 'O1', 'O2', 'C8a',
                       'S1', 'C2', 'N3']
        elif core == 'btz':
            newcols = ['S1', 'C2', 'N3', 'C4', 'C4a', 
                       'C5', 'C6', 'C7', 'C8', 'C8a']
        elif core == 'benzene':
            newcols = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']
        else:
            raise NotImplementedError(f"{core=} not recognized.\
                                    Please use 'btz', 'nitrobtz', or 'benzene'")
        return newcols


    def _column_selector(self, columns, string):
        """Select dataframe columns containing a string.

        Args:
            columns (pd columns): group of columns to analyze.
            string (str): string to search for in columns.

        Returns:
            str: selected column containing the string. Returns None if
                 string not found in columns. 
        """        
        try:
            selected_col = [x for x in columns if string in x][0]
        except IndexError:
            selected_col = None
        return selected_col


    def _get_data(self):
        """Read CSV file and prepare data as internal attribute."""        
        df = pd.read_csv(self.file_props)
        df.rename(columns={"s_m_title": "Name"}, inplace=True)
        df['Name'] = df.Name.apply(self._rename)
        self._rawdata = df


    def _read_atom_idx(self) -> dict[str, list[int]]:
        """Read atom indexes file (file_atoms).

        Returns:
            dict[str, list[int]]: atom indexes per compound.
        """        
        # Read file with chosen atom numbers (i.e. BTZ core)
        lines = []
        with open(self.file_atoms, 'r') as f:
            for line in f:
                # remove blank lines and commented lines
                if line.strip() and not line.startswith('#'):
                    lines.append(line)
        # Define empty dictionary to store results
        indexes = {}
        # Iterate over compound
        for line in lines:
            # Transform string into list of numbers
            name, atom_str = line.split()
            atoms = []
            for string in atom_str.split(','):
                if '-' in string:
                    numrange = [int(i) for i in string.split('-')]
                    atoms.extend(list(range(numrange[0], numrange[1] + 1)))
                else:
                    atoms.append(int(string))
                    
            indexes[name] = atoms
                
        return indexes
    

    def _rename(self, name: str) -> str:
        return name.split("_")[0]


    def _sel_core_atomic_desc(self, data: pd.DataFrame, idxs: list[int]) -> pd.DataFrame:
        """Select atomic descriptors for specified atoms.

        Args:
            data (pd.DataFrame): full set of atomic descriptors.
            idxs (list[int]): atomic indexes.

        Returns:
            pd.DataFrame: reduced dataframe with descriptors for selected atoms.
        """        
        to_remove = data.columns[data.columns.str.contains("Atom_N")]

        data4atoms = data[data.Atom_Number.isin(idxs)].copy()
        data4atoms.sort_values(by="Atom_Number", 
                            key=lambda col: col.map(lambda e: idxs.index(e)),
                            inplace=True)
        data4atoms.drop(columns=to_remove, inplace=True)

        return data4atoms


    def _str2list(self, string):
        """Transforms str of list of data to actual list."""        
        return [s.strip() for s in string[1:-1].split(',')]