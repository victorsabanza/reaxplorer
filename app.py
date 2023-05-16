'''Streamlit app to deploy the Reaxplorer'''
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from model import *
from utils import download_model_uspto480k, download_mcule_molecules, take_random_subset_mols
from streamlit_ketcher import st_ketcher


#start the app
st.title('Reaxplorer')

st.write('''#### Explore the chemical space of purchasable molecules using the molecular transformer.''')

n_mols = st.sidebar.number_input('Number of molecules from catalogue', min_value=1, 
                        max_value=100, value=10, step=1, 
                        help='Number of molecules to select from Mcule database')
random_seed = st.sidebar.number_input('Random seed', min_value=1, 
                        max_value=100, value=33, step=1,
                        help='Random seed to select molecules from Mcule database')

filtering_criteria = st.sidebar.selectbox('Filtering criteria', 
                        ['QED', 'MW', 'LogP'], 
                        help='Scoring functions for product molecules')

#download the model
download_model_uspto480k()

#download the mcule molecules
download_mcule_molecules()

st.success('Downloaded model and molecules')

tab1, tab2 = st.tabs(['Input', 'Output'])

with tab1:
    st.write('''### Draw your molecule of interest''')
    molecule = st_ketcher(value='CCCCC', key='molecule')

    st.write('''### Draw substructures to match in products''')
    substructure = st_ketcher(value='', key='substructure')

    #read only a random subset of n_mols molecules from the .smi file
    mols = take_random_subset_mols(n_mols, random_seed)

    #display molecules as rdkit mol objects
    mols = [Chem.MolFromSmiles(mol) for mol in mols]

    #display images of molecules
    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200, 200))
    st.image(img, use_column_width=True)

    if molecule:
        tokenized = smiles_tokenizer(canonicalize_smiles(molecule))
        smiles_to_txt(tokenized)
        translate()
        predicted = smiles_from_txt()


with tab2:
    st.write(f'Top {n_mols} reactions''')
    
    if predicted:
        rxn = molecule + '>>' + predicted
        rxn = ReactionFromSmarts(rxn)
        img = Draw.ReactionToImage(rxn)
        st.image(img, use_column_width=True)
    
    else:
        pass

# #read the molecules
# molecules = Chem.SmilesMolSupplier('data/molecules.smi', delimiter='\t', titleLine=False,
#                                    nameColumn=-1, smilesColumn=0,)

# #select the molecules
# selected_molecules = st.sidebar.selectbox('Select a molecule', molecules)
