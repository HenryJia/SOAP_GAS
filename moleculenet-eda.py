# Exploratory data analysis of the MoleculeNet dataset
# with the Smooth Overlap of Atomic Positions (SOAP) GA code

import pandas as pd

from quippy import descriptors

from refactoring import GeneParameters
from utils import load_xyz

from tqdm import tqdm

descDict1 = {
    'lower': 1, 'upper': 50, 'centres': '{8, 7, 6, 1, 16, 17, 9}',
    'neighbours': '{8, 7, 6, 1, 16, 17, 9}', 'mu': 0, 
    'mu_hat': 0, 'nu': 2, 'nu_hat': 0, 'mutation_chance': 0.50, 
    'min_cutoff': 1, 'max_cutoff': 50, 'min_sigma': 0.1, 
    'max_sigma': 0.9,
    'message_steps': 0}

num_gens = 100
best_sample, lucky_few, population_size, number_of_children = 4, 2, 12, 4
early_stop = 2
early_number = 3 
min_generations = 5

params1 = GeneParameters(**descDict1)

example_gene_set = params1.make_gene_set()

df = pd.read_csv('BBBP/BBBP_clean.csv')

print(params1)
print(example_gene_set)
print(example_gene_set.gene_parameters)
print(example_gene_set.get_soap_string())

xyz = []
for i, row in tqdm(df.iterrows(), total=df.shape[0]):
    xyz.append(load_xyz('BBBP/xyz/' + str(row['num']) + '.xyz'))

df = df.assign(xyz=xyz)
print(df.head())

#soaps = comp_soaps([params1], df['xyz'])
soaps = []
for i, row in tqdm(df.iterrows(), total=df.shape[0]):
    soaps.append(descriptors.Descriptor(example_gene_set.get_soap_string()).calc(row['xyz']))

df = df.assign(soaps=soaps)

print(df.head())
