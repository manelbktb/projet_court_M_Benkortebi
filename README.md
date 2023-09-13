
# SOLVENT-ACCESSIBLE SURFACE AREA CALCULATION ON PROTEINS
This project aims to calculate the solvent-accessible surface area (SASA) of residues and the entire protein from a PDB file using the Shrake and Rupley method.  The project consists of several Python programs that perform different steps in the calculation and comparison of results.

## Included Programs
The following programs are located in the src directory:
### extract_coords.py
This program extracts information and coordinates of atoms from the PDB file.

### accessibility.py
This program calculates solvent accessibility for each residue of the protein using a specific method.

### comp.py
This program uses the Biopython module (Bio.PDB.SASA) to calculate solvent accessibility and compares the results to values obtained with the accessibility.py program. It also generates a bar chart for visualizing the comparison of values and calculates the percentage identity between the two programs.

### main.py
This main script executes the entire process. It outputs the solvent-accessible surface area in angstroms for each residue using the accessibility.py program, as well as this value for the entire protein. It also executes the calculation of solvent-accessible surface area in angstroms for each residue using the Biopython program, along with the value for the entire protein. The script saves the comparison chart of values in a folder named "results" and displays the percentage identity between the results of both methods.

## Requirements

***Python Version:*** Python 3

***Operating System:*** Unix

The following Python modules are used in this project:

### Biopython
Used for calculating solvent accessibility using the SASA ShrakeRupley method. You can install it in the Conda environment using the command:

```bash
conda install -c conda-forge biopython
```
### pandas
Used for data manipulation and analysis. 

```bash
conda install -c conda-forge pandas
```

### numpy
Used for efficient numerical calculations.

```bash
conda install -c conda-forge numpy
```

### matplotlib
Used for generating graphs, including the bar chart for comparing values. 
```bash
conda install -c conda-forge matplotlib
```


## Execution Instructions
To run this program, follow these steps:
Clond the repository :
```bash
git clone https://github.com/hocinebib/Projet_Court_H-H.git
```

Import the Conda environment with the required dependencies using the project.yml file:

```bash
conda env create -f env_project.yml
```

Activate the newly created Conda environment:
```bash
conda activate env_project
```

Execute the main program main.py by specifying the PDB file name, the path to the PDB file, and the desired number of points on the sphere :
```bash
python src/main.py pdb_name pdb_path nb_points
```

Make sure to replace pdb_name, pdb_path, and nb_points with the appropriate values.

The program will generate detailed results, including the solvent-accessible surface area in angstroms for each residue and the entire protein, as well as a comparison chart of values in the "results" folder. The percentage identity between the results will also be displayed.

## Author

***Manel Benkortebi***
