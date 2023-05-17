'''Streamlit app to deploy the Reaxplorer'''
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from model import *
from utils import download_model_uspto480k, download_mcule_molecules, take_random_subset_mols
from streamlit_ketcher import st_ketcher
from rdkit.Chem.QED import qed

#start the app
st.title('Reaxplorer')

st.write('''#### Explore the chemical space of purchasable molecules using the molecular transformer.''')

n_mols = st.sidebar.number_input('Number of molecules from catalogue', min_value=1, 
                        max_value=500, value=10, step=1, 
                        help='Number of molecules to select from Mcule database')

random_seed = st.sidebar.number_input('Random seed', min_value=1, 
                        max_value=100, value=33, step=1,
                        help='Random seed to select molecules from Mcule database')

filtering_criteria = st.sidebar.selectbox('Filtering criteria', 
                        ['QED', 'LogP'], 
                        help='Scoring functions for product molecules')

n_products = st.sidebar.number_input('Number of products', min_value=1, 
                        max_value=100, value=5, step=1, 
                        help='Number of products to display')

#download the model
download_model_uspto480k()

#download the mcule molecules
download_mcule_molecules()

st.success('Downloaded model and molecules')

tab1, tab2 = st.tabs(['Input', 'Output'])

with tab1:
    st.write('''### Draw your molecule of interest''')
    st.write('''Draw the molecule you want to react with the molecules from the Mcule database
    and click **Apply**''')
    molecule = st_ketcher(value='', key='molecule')

    if filtering_criteria == 'substructure':
        st.write('''### Draw substructures to match in products''')
        substructure = st_ketcher(value='', key='substructure')

    #read only a random subset of n_mols molecules from the .smi file
    mols = take_random_subset_mols(n_mols, random_seed)

    #display molecules as rdkit mol objects
    mols_img = [Chem.MolFromSmiles(mol) for mol in mols]

    #display images of molecules
    st.write('''#### Selected molecules from Mcule database''')
    img = Draw.MolsToGridImage(mols_img, molsPerRow=5, subImgSize=(200, 200))
    st.image(img, use_column_width=True)

with tab2:

    st.write('''Click to predict the reactions''')
    start = st.button('Predict!')

    if start:
        predicted = react_multiproducts(molecule, mols)
        top_reactions, top_scores = score_reactions(predicted, n_products, filtering_criteria)
        
        #return top n products
        st.write('''Top reactions''')

        for i, rxn in enumerate(top_reactions):
            st.write(f'**Top {i+1} reaction**')
            rxn = ReactionFromSmarts(rxn)
            img = Draw.ReactionToImage(rxn)
            st.image(img, use_column_width=True)
            st.write(f'{filtering_criteria}: {top_scores[i]}')


