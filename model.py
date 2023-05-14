"""Functions and stuff related to the model."""

import re
import subprocess
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem import Draw

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
  with open('predictions.txt', 'r') as f:
    s = f.read().replace(' ', '').replace('\n', '')
    return s

def translate():
    subprocess.call(["onmt_translate -model models/USPTO480k_model_step_400000.pt -src input.txt -output predictions.txt  -n_best 1 -beam_size 5 -max_length 300 -batch_size 64"],
                        shell=True)