"""Functions and stuff related to the model."""

import re
import subprocess
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem import MolFromSmiles, MolToSmiles

def canonicalize_smiles(smiles):
  return MolToSmiles(MolFromSmiles(smiles))

def smiles_tokenizer(smiles):
    '''Tokenize SMILES based on regex pattern'''

    SMI_REGEX_PATTERN =  r"(\%\([0-9]{3}\)|\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\||\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"
    smiles_regex = re.compile(SMI_REGEX_PATTERN)
    tokens = [token for token in smiles_regex.findall(smiles)]
    return ' '.join(tokens)

def smiles_to_txt(smiles):
  '''Create .txt file from tokenized smiles'''

  with open('input.txt', 'w') as f:
    f.write(smiles)


def smiles_from_txt():
  #read predicted SMILES from .txt file and return as a list
  with open('predictions.txt', 'r') as f:
    s = f.read().replace(' ', '').split('\n')
    return s

def translate():
    subprocess.call(["onmt_translate -model models/USPTO480k_model_step_400000.pt -src input.txt -output predictions.txt  -n_best 1 -beam_size 5 -max_length 300 -batch_size 64"],
                        shell=True)


def react_multiproducts(mol, smiles):
  '''React a single molecule with a list of SMILES using the Transformer
     and return a list of products.
     
     Args:
       mol: a SMILES string
       smiles: a list of SMILES strings
     Returns:
       final_reactions: a list of reaction SMILES strings'''
  
  #canonicalize the input SMILES 
  pairs = [canonicalize_smiles(mol) + '.' + canonicalize_smiles(s) for s in smiles]
  #tokenize the input SMILES
  tokenized = [smiles_tokenizer(p) for p in pairs]
  #join the tokenized SMILES with a newline character
  tokenized = '\n'.join(tokenized)
  #write the tokenized SMILES to a .txt file
  smiles_to_txt(tokenized)
  #translate the tokenized SMILES to product SMILES
  translate()
  #read the product SMILES from the .txt file
  predicted = smiles_from_txt()[:-1]
  #combine pairs and predicted SMILES
  final_reactions = [pairs[i] + '>>' + predicted[i] for i in range(len(predicted))]

  return final_reactions

if __name__ == '__main__':
    
    smi = 'CCCCC.CCC'

    print(react_multiproducts('CC', ['SCCC', 'CCCCCCC', 'c1ccccc1']))