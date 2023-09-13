"""
Solvent Accessibility Comparison

This program calculates and compares solvent accessibility values using the Biopython's program.

Functions:
    - accessibility_Shrake(filename):
        Calculate solvent accessibility using the ShrakeRupley algorithm.

    - comparison(df_access, df_access_shrake):
        Compare the solvent accessibility values between two DataFrames.

    - comparison_plot(df_access, df_access_shrake):
        Create a comparison plot of solvent accessibility values.

     
    - comparison_identity(df_access, df_access_shrake):
        Calculate the percentage identity of ASA values between two programs.


Usage:
1. Load a protein structure into a Pandas DataFrame using the 'extract_pdb' function.
2. Call the 'solvent_accessibility' function to obtain the dataframe of asa informations.
3. Call the 'accessibility_Shrake' to calculate solvent accessibility using the Biopython's program using Shrake & Rupley method.
4. Call the "comparison_plot" and "comparison_plot" functions to obtain respectively a bar diagram of comparison of the asa for each residue 
and the percentage of identity of this value between the two algorithms


Example:
    import extract_coords
    import accessibility

    # Load the protein structure into a DataFrame.
    filename = "protein"
    atom_df = extract_coords.extract_pdb("filename")

    # Calculate solvent accessibility and group the results by residue.
    result_df = accessibility.solvent_accessibility(atom_df)

    # Calculate solvent accessibility by residue with Shrake and Rupley algorithm.
    result_shrake_df = accessibility_Shrake(filename)

    # Create a comparison plot of solvent accessibility values.
    comparison_plot(result_df, result_shrake_df)


Requirements:
- Biopython
- Pandas
- Matplotlib
- extract_coords
- accessibility
"""

import pandas as pd
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import matplotlib.pyplot as plt


def accessibility_Shrake(pdb_name, file_path):
    """
    Calculate solvent accessibility using the ShrakeRupley algorithm.

    Parameters:
    pdb_name (str): The name of the PDB file.
    file_path (str): The path of the PDB file.

    Returns:
    DataFrame: A DataFrame containing residue ID, residue name, and surface accessibility.
    """
    VDW_RADIUS = {'C': 2, 'N': 1.5, 'O': 1.4, 'S': 1.85}

    p = PDBParser(QUIET=True)


    structure = p.get_structure(pdb_name, file_path)


    sr = ShrakeRupley(radii_dict=VDW_RADIUS)


    sr.compute(structure, level="R")
    
    atom_data=[]
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    atom_info = {
                                'id_residue': residue.id[1],
                                'residue': f"{residue.get_resname()}",
                                'surface_accessible_shrake': round(residue.sasa,2)
                     }
                    atom_data.append(atom_info)
    dataf=pd.DataFrame(atom_data)
    dataf = dataf.set_index('id_residue')

    sr.compute(structure, level = "S")
    structure.sasa
    return dataf




def comparison(df_access, df_access_shrake):
    """
    Compare the solvent accessibility values between two DataFrames.

    Parameters:
    df_access (DataFrame): DataFrame containing accessibility values.
    df_access_shrake (DataFrame): DataFrame containing accessibility values calculated using Biopython's program.

    Returns:
    DataFrame: DataFrame with absolute differences added.
    """
    
    absolute_diff = abs(df_access['surface_accessible'] - df_access_shrake['surface_accessible_shrake'])
    absolute_diff = pd.DataFrame(absolute_diff)
    absolute_diff.columns = ["absolute_difference"]
    diff_df = df_access_shrake.copy()
    diff_df["surface_accessible"] = tuple(df_access['surface_accessible'].to_list())
    diff_df["absolute_difference"] = tuple(absolute_diff['absolute_difference'].to_list())
    return diff_df



def comparison_plot(df_access, df_access_shrake):
    """
    Create a comparison plot of solvent accessibility values.

    Parameters:
    df_access (DataFrame): DataFrame containing accessibility values.
    df_access_shrake (DataFrame): DataFrame containing accessibility values calculated using Biopython's program.
    Returns:
    fig (figure): Figure object containing the bar chart of ASA (Accessible Surface Area) values given by the two programs for each residue
    """
    
    diff_df = comparison(df_access, df_access_shrake)

    plt.figure(figsize=(10, 6)) 
    bar_width = 0.35
    index = diff_df.index

    plt.bar(index - bar_width/2, diff_df['surface_accessible'], bar_width, label='ASA Value with our program', color='cornflowerblue')
    plt.bar(index + bar_width/2, diff_df['surface_accessible_shrake'], bar_width, label="ASA Value with Biopython's program", color='indianred')

    plt.xlabel("Residue ID")
    plt.ylabel("ASA Value (Å)")
    plt.title("Comparison of ASA Values by Residue ID")
    plt.xticks(index)
    plt.legend()
    plt.grid(True)

 



def comparison_identity(df_access, df_access_shrake):
    """
    Calculate the percentage identity of ASA values between two programs.

    Parameters:
    df_access (DataFrame): DataFrame containing accessibility values.
    df_access_shrake (DataFrame): DataFrame containing accessibility values calculated using ShrakeRupley.

    Returns:
    str: A string containing the percentage identity of ASA values.
    """
    
    diff_df = comparison(df_access, df_access_shrake)
    diff_tot = diff_df['absolute_difference'].sum()
    asa_tot_shrake =  diff_df['surface_accessible_shrake'].sum()
    identity_pourcent = (1-diff_tot/asa_tot_shrake)*100
    
    return f"The percentage of identity of the asa between these two programs is {round(identity_pourcent,2)} %"

if __name__ == "__main__":
    import comp
