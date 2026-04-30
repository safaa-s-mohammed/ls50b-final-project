# Import necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats

# Provides fit_tobit_model()
import run_tobit_model_w7 as run_tobit_model

df = pd.read_csv('VRC07_kds.csv')

# Define variables reused throughout this analysis program
site_cols = [f'site_{i}' for i in range(10)]
antigen_col = 'antigen'
kd_column = 'nlog10_Kd'
left_censor = 6

# Define antigens
antigen1 = 'SF162'
antigen2 = 'CH505TF'

# Relation to Week 7 assignment: Found order 3 fit  data the best --> means there is epistatic effect --> if I fit order 1, it wont model the epistatic effects
# Therefore, it wont give us realistic data --> won't capture true structure of the data

# Looking at additive effects --> pull out on contribution to the binding affinity (ignoring epistatic effects)
# Find out what the effect on any given mutation is

# Full-data fit (no held-out fold)
_, _, _, df_params_antigen1 = run_tobit_model.fit_tobit_model(
    df=df,
    kd_column=kd_column,
    order=3,
    antigen=antigen1,
    antigen_column=antigen_col,
    site_columns=site_cols,
    genotype_id_column='id',
    left_censor=left_censor,
    seed=42,
    verbose=False,
    alpha_l1 = 0.1
)

# Full-data fit (no held-out fold)
_, _, _, df_params_antigen2 = run_tobit_model.fit_tobit_model(
    df=df,
    kd_column=kd_column,
    order=4,
    antigen=antigen2,
    antigen_column=antigen_col,
    site_columns=site_cols,
    genotype_id_column='id',
    left_censor=left_censor,
    seed=42,
    verbose=False,
    alpha_l1 = 0.1
)

# Add new column with distances matching up with sites
df_params_antigen1 = df_params_antigen1[df_params_antigen1['order'] == 1]

# Load the saved distances
df_distances = pd.read_csv('min_distances.csv')

df_params_antigen1['SF162_g_values'] = df_distances['SF162_germline'].values
df_params_antigen1['SF162_m_values'] = df_distances['SF162_mature'].values

df_params_antigen2 = df_params_antigen2[df_params_antigen2['order'] == 1]

df_params_antigen2['CH505TF_g_values'] = df_distances['CH505TF_germline'].values
df_params_antigen2['CH505TF_m_values'] = df_distances['CH505TF_mature'].values

plt.scatter(df_params_antigen1['SF162_g_values'],df_params_antigen1['coefficient'], label='germline structure')
plt.title('additive effects (SF162)')
plt.ylabel('coefficients')
plt.xlabel('min dist from a given residue to the antigen (angstrom (Å))')
plt.scatter(df_params_antigen1['SF162_m_values'],df_params_antigen1['coefficient'], label='mature structure')
plt.legend()

plt.show()
plt.scatter(df_params_antigen2['CH505TF_g_values'],df_params_antigen2['coefficient'], label='germline structure')
plt.title('additive effects (CH505TF)')
plt.ylabel('coefficients')
plt.xlabel('min dist from a given residue to the antigen (angstrom (Å))')
plt.scatter(df_params_antigen2['CH505TF_m_values'],df_params_antigen2['coefficient'], label='mature structure')
plt.legend()