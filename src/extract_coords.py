"""
PDB Atom Information Extractor

This program extracts atom information from a Protein Data Bank (PDB) file using Biopython.
It processes the PDB file, excluding hydrogen atoms, and creates a Pandas DataFrame
containing information about each atom, including the residue, residue ID, atom name, atom type,
and Cartesian coordinates (x, y, z).

Functions:
    - extract_pdb(filename):
        Extract atom information from a PDB file and return it as a Pandas DataFrame.
Usage:
    - Load a protein structure into a Pandas DataFrame using the 'extract_pdb' function.

Example:
    - To extract atom information from a PDB file 'protein.pdb':
        df = extract_pdb('protein.pdb')
        print(df)

Requirements:
    - Biopython
    - Pandas
"""

from Bio import PDB
from Bio.PDB import *
import pandas as pd

def extract_pdb(pdb_name, file_path):
    """
    Extract atom information from a PDB file.

    Parameters:
        pdb_name (str): The name of the PDB file.
        file_path (str): The path of the PDB file.

    Returns:
        atom_data: A Pandas DataFrame containing atom information.
    """
    # Initialize the PDB parser
    parser = PDBParser(QUIET=True)

    # Load the structure from the PDB file
    structure = parser.get_structure(pdb_name, file_path)

    # List to store atom information
    atom_data = []

    # Iterate through the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):  # Check if the residue is an amino acid
                    for atom in residue:
                        if atom.element != 'H':  # Exclude hydrogen atoms
                            atom_info = {
                                'residue': f"{residue.get_resname()}",
                                'id_residue': residue.id[1],
                                'atom_name': atom.get_name(),
                                'atom_type': atom.element,
                                'x': atom.get_coord()[0],
                                'y': atom.get_coord()[1],
                                'z': atom.get_coord()[2],
                            }
                            atom_data.append(atom_info)

    # Create a Pandas DataFrame from the atom data
    return pd.DataFrame(atom_data)

if __name__ == "__main__":
    import extract_coords
