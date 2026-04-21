"""
Structure Data Extraction from PDB/AlphaFold Protein Models
"""

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
# Some form of AlphaFold file initialization once I figure that out
antibody_chains = ["H", "L"] # Creates sections for both heavy and light chains
antigen_chains = ["G"]
output_distance_file = "distance_matrix.csv"
output_interfact_contacts = "interface_contacts.csv"

# Initialize PyMOL
pymol.finish_launching(['pymol', '-cq'])


# Basic approach: calculate pairwise Euclidean distances
# between Cα in a protein structure
