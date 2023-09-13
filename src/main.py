"""
Main script to perform solvent accessibility calculations and comparisons.

This script reads a PDB file, calculates solvent accessibility using Shrake & Rupley method 
with our program and Biopython's program, and generates a comparison plot while also calculating 
the percentage identity between the methods.
"""

import argparse
import extract_coords as coor
import accessibility as access
import pandas as pd
import comp 
import time
import matplotlib.pyplot as plt



if __name__ == "__main__":
    # Start measuring execution time
    start_time = time.time()
    
    # Define command-line argument parser
    parser = argparse.ArgumentParser()

    # Define command-line arguments   
    parser.add_argument("pdb_name",
                       help="the pdb name of the protein",
                       type = str)

    parser.add_argument("file_path",
                        help="the path of the pdb file",
                        type=str)
   
    parser.add_argument("n_points",
                        help="the number of points on the spheres",
                        type=int)
                        
    # Parse command-line arguments
    args = parser.parse_args()
    
    # Extract parsed arguments
    pdb_name = args.pdb_name
    file_path = args.file_path
    n_points = args.n_points
    if n_points < 2:
        raise ValueError("Number of points must be 2 or greater. Please choose a larger value for n_points.")

   
    print("Step 1 : Reading pdb file and extracting atoms coordinates\n")

    df_atoms = coor.extract_pdb(pdb_name, file_path)
   
    print("Step 2 : Calculating accessible surface area for each residue and for all protein\n")

    asa = access.solvent_accessibility(df_atoms, n_points)
    print(asa)
    
    print(f"\nWith our program, the surface accessible to the solvent of the protein is {round(asa['Surface_accessible'].sum(),2)} Å.\n")
    
    print("Step 3 : Calculating accessible surface of residus and protein using the ShrakeRupley method\n")
   
    shrake = comp.accessibility_Shrake(pdb_name, file_path)

    print("Step 4 : Comparison of accessible surface between our method and Shrake & Rupley method\n")      
    
    diff = comp.comparison(asa, shrake)
    print(diff)
    
    print(f"\nWith Biopython's program, the surface accessible to the solvent of the protein is {round(shrake['surface_accessible_shrake'].sum(),2)} Å.\n")

    # Generate and save a comparison plot
    comp.comparison_plot(asa, shrake)
    plt.savefig('results/plot.png')

    print(f"\nYou will find the plot in your repertory!\n")
    
    # Calculate and print the percentage identity of ASA values between the two methods
    identity = comp.comparison_identity(asa, shrake)
    
    print(identity)

    # End measuring execution time
    end_time = time.time()

    execution_time = end_time - start_time

    print(f"\nExecution time : {round(execution_time)} seconds")
