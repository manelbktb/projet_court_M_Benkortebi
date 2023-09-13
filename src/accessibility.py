"""
Solvent Accessibility Calculation

This program calculates solvent accessibility for atoms in a protein structure and groups the results by residue.

Functions:
    - calculate_distance(x1, y1, z1, x2, y2, z2)
        Calculates the Euclidean distance between two points in 3D space.
    - find_neighbors(atom, atom_df)
        Finds neighboring atoms within a threshold distance of a given atom.
    - distribute_points_on_sphere(num_points, radius)
        Distributes points evenly on the surface of a sphere.
    - place_sphere_around_atom(atom_x, atom_y, atom_z, num_points, sphere_radius)
        Places points on a sphere centered around an atom's coordinates.
    - radius_sphere(atom_type)
        Calculates the effective radius of a sphere based on atom type.
    - solvent_accessibility(atom_df, points_nb=100)
        Calculates solvent accessible surface and groups the results by residue.

Usage:
1. Load a protein structure into a Pandas DataFrame using the 'extract_pdb' function.
2. Call the 'solvent_accessibility' function to calculate solvent accessibility.

Example:
    import extract_coords
    import accessibility

    # Load the protein structure into a DataFrame
    atom_df = extract_coords.extract_pdb("protein.pdb")

    # Calculate solvent accessibility and group the results by residue
    result_df = accessibility.solvent_accessibility(atom_df)

    print(result_df)

Requirements:
- math
- Pandas
- NumPy
- extract_coords
"""



import math
import numpy as np
import pandas as pd


def calculate_distance(x1, y1, z1, x2, y2, z2):
    """
    Calculate the Euclidean distance between two points in 3D space.

    Parameters:
        x1 (float): X-coordinate of the first point.
        y1 (float): Y-coordinate of the first point.
        z1 (float): Z-coordinate of the first point.
        x2 (float): X-coordinate of the second point.
        y2 (float): Y-coordinate of the second point.
        z2 (float): Z-coordinate of the second point.

    Returns:
        float: The Euclidean distance between the two points.
    """
    
    distance = ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5
    
    return distance


def find_neighbors(atom, atom_df):
    """
    Find neighboring atoms within a threshold distance of the given atom.

    Parameters:
        atom (pd.Series): A Pandas Series representing the atom of interest with 'x', 'y', and 'z' coordinates.
        atom_df (pd.DataFrame): A Pandas DataFrame containing atom information with 'x', 'y', 'z' coordinates.

    Returns:
        pd.DataFrame: A subset of the input DataFrame containing neighboring atoms.
    """
    # Van der Waals radii for common elements and solvent radius
    VDW_RADIUS = {'C': 2, 'N': 1.5, 'O': 1.4, 'S': 1.85}
    SOLVENT_RADIUS = 1.4

    # Calculate the threshold distance for finding neighbors
    Thrsd = math.ceil(2 * VDW_RADIUS['C'] + SOLVENT_RADIUS)

    # Extract the coordinates of the atom
    x, y, z = atom['x'], atom['y'], atom['z']

    # Calculate the distance from the atom to all other atoms in atom_df
    neighbors = atom_df.copy()  # Make a copy of the DataFrame to avoid modifying the original
    neighbors['Distance'] = neighbors.apply(lambda row: calculate_distance(x, y, z, row['x'], row['y'], row['z']), axis=1)

    # Filter neighbors based on the threshold distance and exclude the atom itself
    neighbors = neighbors[(neighbors['Distance'] <= Thrsd) & (neighbors.index != atom.name)]

    return neighbors




def distribute_points_on_sphere(num_points, radius):
    """
    Distribute points evenly on the surface of a sphere.

    Parameters:
        num_points (int): Number of points to distribute.
        radius (float): Radius of the sphere.

    Returns:
        np.array: An array of (x, y, z) coordinates representing the points on the sphere.
    """
    points = []
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))  # Golden angle
    
    for i in range(num_points):
        z = 1.0 - (i / float(num_points - 1)) * 2  # Varying height
        radius_at_height = np.sqrt(1.0 - z * z) * radius
        theta = golden_angle * i  # Angle around the z-axis
        x = np.cos(theta) * radius_at_height
        y = np.sin(theta) * radius_at_height
        points.append((x, y, z * radius))

    return np.array(points)

def place_sphere_around_atom(atom_x, atom_y, atom_z, num_points, sphere_radius):
    """
    Place points on a sphere centered around an atom's coordinates.

    Parameters:
        atom_x (float): X-coordinate of the atom.
        atom_y (float): Y-coordinate of the atom.
        atom_z (float): Z-coordinate of the atom.
        num_points (int): Number of points to distribute on the sphere.
        sphere_radius (float): Radius of the sphere.

    Returns:
        list: A list of (x, y, z) coordinates representing the points on the sphere centered around the atom.
    """

    sphere_points = distribute_points_on_sphere(num_points, sphere_radius)
    centered_points = [(x + atom_x, y + atom_y, z + atom_z) for x, y, z in sphere_points]
    return centered_points

def radius_sphere(atom_type):
    """
    Calculate the effective radius of a sphere based on the atom type.

    Parameters:
        atom_type (str): The type of the atom within the residue ('C', 'N', 'O' or 'S').

    Returns:
        float: The effective radius of the sphere for the given atom type.
    """
    VDW_RADIUS = {'C': 2, 'N': 1.5, 'O': 1.4, 'S': 1.85}
    SOLVENT_RADIUS = 1.4

    # Calculate the effective radius by adding the Van der Waals radius of carbonate and solvent radius
    radius = VDW_RADIUS.get(atom_type, 0) + SOLVENT_RADIUS

    return radius



def solvent_accessibility(atom_df, points_nb=100):
    """
    Calculate solvent accessibility for atoms in a DataFrame and group the results by residue.

    Parameters:
        atom_df (pd.DataFrame): DataFrame containing atom information.
        points_nb (int, optional): Number of points to distribute on the sphere. Default is 100.

    Returns:
        pd.DataFrame: DataFrame containing information on the solvent accessible surface grouped by residue.
    """
    
    ratio_accessibility = []  # List to store percentage accessibility for each atom
    accessible_surface_atoms = []  # List to store accessible surface area for each atom
    total_surface_atoms = []  # List to store total surface area for each atom
    
    # Create a column for sphere radius based on atom type in atom_df
    atom_df['sphere_radius'] = atom_df['atom_type'].apply(radius_sphere)

    for atom_index, atom_row in atom_df.iterrows():
        ref_coords = atom_df.loc[atom_index, ['x', 'y', 'z']]  # Reference atom coordinates
        neighbors = find_neighbors(atom_row, atom_df)  # Find neighboring atoms

        # Create an array of points on the sphere centered around the atom
        sphere_points = place_sphere_around_atom(*ref_coords, points_nb, atom_row['sphere_radius'])

        exposed_points = 0  # Counter for exposed points on the sphere

        # Check if each point on the sphere is exposed or buried
        for i in range(len(sphere_points)):
            is_exposed = True
            for _, neighb_row in neighbors.iterrows():
                threshold = radius_sphere(neighb_row['atom_type'])
                distance = calculate_distance(*sphere_points[i], *neighb_row[['x', 'y', 'z']])
                if distance <= threshold:
                    is_exposed = False
                    break
            if is_exposed:
                exposed_points += 1

        ratio = exposed_points / points_nb  # Calculate ratio accessibility for each atom
        ratio_accessibility.append(ratio)

        # Calculate accessible and total surface area for the atom
        surface_accessible = ratio * 4 * np.pi * (atom_row['sphere_radius'])**2
        surface_total = 4 * np.pi * (atom_row['sphere_radius'])**2
        accessible_surface_atoms.append(surface_accessible)
        total_surface_atoms.append(surface_total)

    # Add columns for accessible and total surface area to the DataFrame
    atom_df["surface_accessible"] = accessible_surface_atoms
    atom_df["surface_total"] = total_surface_atoms

    # Create a copy of the DataFrame for residue-level calculations
    access_residu = atom_df.copy()
    
    # Calculate the accessible and total surface area for each residue
    access_residu = round(access_residu.groupby(['id_residue', 'residue'])[['surface_accessible', 'surface_total']].sum(), 2)
    
    # Calculate the percentage of surface accessibility for each residue
    access_residu['pourcent_surface_accessible'] = round(access_residu['surface_accessible'] / access_residu['surface_total'] * 100, 2)
    
    # Drop the 'Surface_total' column from the resulting DataFrame
    access_residu = access_residu.drop('surface_total', axis=1)
    
    
    access_total = access_residu['surface_accessible'].sum()
    access_tot = f"L'accessibilité totale de la protéine est de {access_total} Angstroms"
    return access_residu    


if __name__ == "__main__":
    import accessibility
