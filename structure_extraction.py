
# Structure Data Extraction from PDB/AlphaFold Protein Models

import numpy as np
import pandas as pd
import pymol
# Opens a gateway for PyMOL to interact with Python scripts
from pymol import cmd
# Itertools is a library of functions that creates iterators
# for more efficient looping, product is a function that produces
# all possible ordered combinations of elements (aka a Cartesian product) 
from itertools import product

# Configuration, variable initialization
pdb_file = r"C:\Users\safaa\Downloads\4OLU.pdb"
antibody_chains = ["H", "L"] # Creates sections for both heavy and light chains
antigen_chains = ["G"]
output_distance_file = "distance_matrix.csv"
output_interface_contacts = "interface_contacts.csv"

# Initialize PyMOL
pymol.finish_launching(['pymol', '-cq'])
cmd.load(pdb_file, "complex")

# Basic approach: calculate pairwise Euclidean distances
# between Cα in a protein structure

def get_alpha_carbons(chains):
    atoms = []
    # Joins each chain together to later use as a string input for
    # the chains of interest
    section = " or ".join([f"chain {c}" for c in chains])
    # Iterates over the chain to look for the Cα and adds its position, its residue information
    # to the dictionary
    cmd.iterate_state(
        1,
        f"({section}) and name CA",
        "atoms.append({'chain': chain, 'resindex': resi, 'resiname': resn, 'x': x, 'y': y, 'z': z})",
        # Maps variables from the Python script to atom properties in the iteration
        space = {'atoms': atoms}
    )
    return atoms

antibody_atoms = get_alpha_carbons(antibody_chains)
antigen_atoms = get_alpha_carbons(antigen_chains)
all_atoms = antibody_atoms + antigen_atoms

# Checks for the length of the antibody chains and antigen chains
print(f"Antibody residues: {len(antibody_atoms)}")
print(f"Antigen residues: {len(antigen_atoms)}")
print(f"Total residues: {len(all_atoms)}")

# Creates the distance matrix
n = len(all_atoms)
coordinates = np.zeros((n, 3))
labels = []

for i, a in enumerate(all_atoms):
    coordinates[a, 0] = a['x']
    coordinates[a, 1] = a['y']
    coordinates[a, 2] = a['z']

# Uses dictionary comprehension to create a dictionary from an iterable
# object in a single line
labels = [f"{a['chain']}_{a['resindex']}_{a['resiname']}" for a in all_atoms]

# Calculates the distance matrix
distance_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        dx = coordinates[i, 0] - coordinates[j, 0]
        dy = coordinates[i, 1] - coordinates[j, 1]
        dz = coordinates[i, 2] - coordinates[j, 2]
        distance_matrix[i, j] = np.sqrt((dx**2) + (dy**2) + (dz**2))

df_distances = pd.DataFrame(distance_matrix, index = labels, columns = labels)
df_distances.to_csv(output_distance_file)
print(f"Saved distance matrix! Name of the file is: {output_distance_file}")

# Identify interface contacts
contacts = []
for antibody, antigen in product(antibody_atoms, antigen_atoms):
    dx = antibody['x'] - antigen['x']
    dy = antibody['y'] - antigen['y']
    dz = antibody['z'] - antigen['z']
    distance = np.sqrt((dx**2) + (dy**2) + (dz**2))
    if distance <= 12:
        contacts.append({
            'antibody_chain': antibody['chain'],
            'antibody_resindex': antibody['resindex'],
            'antibody_resiname': antibody['resiname'],
            'antigen_chain': antigen['chain'],
            'antigen_resindex': antigen['resindex'],
            'antigen_resiname': antigen['resiname'],
            'distance': round(distance, 3)
        })

df_interface = pd.DataFrame(contacts).sort_values("distance")
df_interface.to_csv(output_interface_contacts, index = False)
print(f"Saved interface contacts file! File name is: {output_interface_contacts}")

cmd.quit()