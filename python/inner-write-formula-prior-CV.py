import argparse
import numpy as np
import os
import pandas as pd
from rdkit.Chem import Descriptors, rdMolDescriptors
from tqdm import tqdm

# set working directory
git_dir = os.path.expanduser("~/git/invalid-smiles-analysis")
python_dir = git_dir + "/python"
os.chdir(python_dir)

# import functions
from functions import clean_mol, clean_mols, read_smiles

# suppress rdkit errors
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

### dynamically build CLI
parser = argparse.ArgumentParser()
## build the CLI
grid_file = git_dir + '/sh/grids/write-formula-prior.txt'
grid = pd.read_csv(grid_file, sep='\t')
for arg_name in list(grid):
    param_name = '--' + arg_name
    param_dtype = str(grid[arg_name].dtype)
    # convert to pandas
    param_type = {'object': str,
                  'int64': int,
                  'float64': float,
                  'bool': str
                  }[param_dtype]
    parser.add_argument(param_name, type=param_type)

# parse all arguments
args = parser.parse_args()
print(args)

if not os.path.isdir(args.output_dir):
    os.makedirs(args.output_dir)

# read training and test sets
train_smiles = read_smiles(args.train_file)
test_smiles = read_smiles(args.test_file)

# calculate the masses and molecular formulas of the train and test sets
train_mols = clean_mols(train_smiles)
train_masses = [round(Descriptors.ExactMolWt(mol), 4) for mol in train_mols]
train_fmlas = [rdMolDescriptors.CalcMolFormula(mol) for mol in train_mols]
test_mols = clean_mols(test_smiles)
test_masses = [round(Descriptors.ExactMolWt(mol), 4) for mol in test_mols]
test_fmlas = [rdMolDescriptors.CalcMolFormula(mol) for mol in test_mols]
test = pd.DataFrame({'smiles': test_smiles})
test = test.assign(mass=test_masses, formula=test_fmlas)
test = test.assign(mass_known=test['mass'].isin(train_masses))
test = test.assign(formula_known=test['formula'].isin(train_fmlas))

# create training set
train = pd.DataFrame({'smiles': train_smiles, 
                      'mass': train_masses,
                      'formula': train_fmlas})
# assign frequencies as NAs
train = train.assign(size=np.nan)

# read PubChem file
pubchem = pd.read_csv(args.pubchem_file, delimiter='\t', header=None,
                      names=['smiles', 'mass', 'formula'])
# assign frequencies as NAs
pubchem = pubchem.assign(size=np.nan)

# read sample file from the generative model
gen = pd.read_csv(args.sample_file)

# set up outputs
rank_df = pd.DataFrame()

# iterate through PubChem vs. generative model
inputs = {'model': gen.assign(source='model'),
          'PubChem': pubchem.assign(source='PubChem'),
          'train': train.assign(source='train')}

# match on formula and mass
for key, query in inputs.items():
    for row in tqdm(test.itertuples(), total=test.shape[0]):
        # get formula and exact mass
        query_mol = clean_mol(row.smiles)
        query_mass = Descriptors.ExactMolWt(query_mol)
        
        # compute 10 ppm range
        min_mass = (-args.err_ppm / 1e6 * query_mass) + query_mass
        max_mass = (args.err_ppm / 1e6 * query_mass) + query_mass
        
        # get matches to mass
        if key != "model":
            # get matches to mass and shuffle them to random order
            matches = query[query['mass'].between(min_mass, max_mass)].\
                sample(frac=1)
            matches = matches.assign(rank=np.arange(matches.shape[0]))
        else:
            # generative model: sort descending by mass
            matches = query[query['mass'].between(min_mass, max_mass)].\
                sort_values('size', ascending=False)
            matches = matches.assign(rank=np.arange(matches.shape[0]))
        
        # add 'target_' to all columns in matches
        matches.columns = 'target_' + matches.columns
        
        # get the rank of the correct formula
        rank = matches[matches['target_formula'] ==
                       row.formula][['target_size', 'target_rank',
                                    'target_source']]
        if rank.shape[0] > 1:
            rank = rank.head(1)
        elif rank.shape[0] == 0:
            rank = pd.DataFrame({'target_size': np.nan,
                                 'target_rank': np.nan,
                                 'target_source': key}, index=[0])
        
        # assign number of candidates
        n_formulas = len(matches.target_formula.unique())
        rank = rank.assign(n_candidates=n_formulas)

        # now, append rank to the growing data frame
        rank_row = pd.concat([pd.DataFrame([row]).reset_index(drop=True),
                              rank.reset_index(drop=True)], axis=1)
        rank_df = rank_df.append(rank_row)

# write to output files
rank_df.to_csv(args.output_dir + '/CV-ranks.csv.gz', index=False,
               compression='gzip')