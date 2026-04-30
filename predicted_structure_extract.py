# Import necessary libraries 
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqUtils import seq1 

# Parse mmCIF file
parser = MMCIFParser(QUIET = True)
# mmCIF file with the mature site model will provide the mature structure of each antibody-antigen complex
SF162_mature = parser.get_structure('protein', 'fold_sf162_mature_site_model_0_.cif')
CH505TF_mature = parser.get_structure('protein', 'fold_ch505tf_mature_site_model_0.cif')

# Getting structures for germline model of antibody-antigen complex
# SF162 germline sites
SF162_g0 = parser.get_structure('protein', 'fold_sf162_germline_site_0_model_0.cif')
SF162_g1 = parser.get_structure('protein', 'fold_sf162_germline_site_1_model_0.cif')
SF162_g2 = parser.get_structure('protein', 'fold_sf162_germline_site_2_model_0.cif')
SF162_g3 = parser.get_structure('protein', 'fold_sf162_germline_site_3_model_0.cif')
SF162_g4 = parser.get_structure('protein', 'fold_sf162_germline_site_4_model_0.cif')
SF162_g5 = parser.get_structure('protein', 'fold_sf162_germline_site_5_model_0.cif')
SF162_g6 = parser.get_structure('protein', 'fold_sf162_germline_site6_model_0.cif')
SF162_g7 = parser.get_structure('protein', 'fold_sf162_germline_site_7_model_0.cif')
SF162_g8 = parser.get_structure('protein', 'fold_sf162_germline_site_8_model_0.cif')
SF162_g9 = parser.get_structure('protein', 'fold_sf162_germline_site_9_model_0.cif')
# CH505TF germline sites
CH505TF_g0 = parser.get_structure('protein', 'fold_ch505tf_germline_site_0_model_0.cif')
CH505TF_g1 = parser.get_structure('protein', 'fold_ch505tf_germline_site_1_model_0.cif')
CH505TF_g2 = parser.get_structure('protein', 'fold_ch505tf_germline_site_2_model_0.cif')
CH505TF_g3 = parser.get_structure('protein', 'fold_ch505tf_germline_site_3_model_0.cif')
CH505TF_g4 = parser.get_structure('protein', 'fold_ch505tf_germline_site_4_model_0.cif')
CH505TF_g5 = parser.get_structure('protein', 'fold_ch505tf_germline_site_5_model_0.cif')
CH505TF_g6 = parser.get_structure('protein', 'fold_ch505tf_germline_site6_model_0.cif')
CH505TF_g7 = parser.get_structure('protein', 'fold_ch505tf_germline_site_7_model_0.cif')
CH505TF_g8 = parser.get_structure('protein', 'fold_ch505tf_germline_site_8_model_0.cif')
CH505TF_g9 = parser.get_structure('protein', 'fold_ch505tf_germline_site_9_model_0.cif')

# Defines a function to get the minimum distance between sites on the antibody and the closest residue(s) on the antigen of interest
def find_min_distance(structure, residue_interest):
    # Returns a generator over all chains in the structure, saves as a list stuck to chains
    chains = list(structure.get_chains())

    # Want to compare the site we changed to all the residues on the antigen
    # First chain is the residue source --> heavy chain
    H_chain = chains[0]

    # Chain C = target chain
    C_chain = None
    for chain in chains: 
        if chain.id == "C":
            C_chain = chain
            break

    # Find the residue of interest in H chain
    res = None
    for residue in H_chain:
        if residue.id[0] == " " and residue.id[0] == residue_interest:
            res = residue
            break

    # Use alpha-carbon atom for distance comparison
    # Find distance between residue of interest in heavy chain and antigen residues, 
    # then calculate minimum
    atom1 = res["CA"]

    all_residue = []
    for residue in C_chain:
        if residue.id[0] == " " and "CA" in residue:
            atom2 = residue["CA"]
            # Biopython calculates Euclidean distance
            distance = atom1 - atom2
            all_residue.append(distance)
    min_distance = min(all_residue)

    return min_distance

# Finds minimum distance for each site in the mature model, using the amture structure for each antigen
# Store each result in a variable
# SF162
SF162_m0 = find_min_distance(SF162_mature, 3)
SF162_m1 = find_min_distance (SF162_mature, 18)
SF162_m2 = find_min_distance(SF162_mature, 35)
SF162_m3 = find_min_distance(SF162_mature, 45)
SF162_m4 = find_min_distance(SF162_mature, 51)
SF162_m5 = find_min_distance(SF162_mature, 52)
SF162_m6 = find_min_distance(SF162_mature, 55)
SF162_m7 = find_min_distance(SF162_mature, 57)
SF162_m8 = find_min_distance(SF162_mature, 59)
SF162_m9 = find_min_distance(SF162_mature, 127)
# CH505TF
CH505TF_m0 = find_min_distance(CH505TF_mature, 3)
CH505TF_m1 = find_min_distance(CH505TF_mature, 18)
CH505TF_m2 = find_min_distance(CH505TF_mature, 35)
CH505TF_m3 = find_min_distance(CH505TF_mature, 45)
CH505TF_m4 = find_min_distance(CH505TF_mature, 51)
CH505TF_m5 = find_min_distance(CH505TF_mature, 52)
CH505TF_m6 = find_min_distance(CH505TF_mature, 55)
CH505TF_m7 = find_min_distance(CH505TF_mature, 57)
CH505TF_m8 = find_min_distance(CH505TF_mature, 59)
CH505TF_m9 = find_min_distance(CH505TF_mature, 127)

