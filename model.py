"""Functions and stuff related to the model."""

import re
import subprocess
from rdkit import Chem
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem.QED import qed
#import logp from rdkit
from rdkit.Chem.Crippen import MolLogP

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


def score_reactions(predictions, top_n,  criterion='QED'):
  '''Score reactions according to selected criterion
  
  '''
  #get product SMILES of predictions and get indices of the top-n based on criterion
  products = [Chem.MolFromSmiles(p.split('>>')[1]) for p in predictions]
  
  if criterion == 'QED':
    scores = [round(qed(p),3) for p in products if p is not None]
  elif criterion == 'LogP':
    scores = [round(MolLogP(p), 3) for p in products if p is not None]
  
  top_n_idx = sorted(range(len(scores)), key=lambda i: scores[i])[-top_n:]
  #invert top_n_idx to get descending order
  top_n_idx = top_n_idx[::-1]
  top_n_scores = [scores[i] for i in top_n_idx]
  top_n_predictions = [predictions[i] for i in top_n_idx]

  return top_n_predictions, top_n_scores


if __name__ == '__main__':
    
    smi = 'CCCCC.CCC'

    p = react_multiproducts('CC', ['SCCC', 'CCCCCCC', 'c1ccccc1', 'n1cccnc1'])

    print(score_reactions(p, 3))