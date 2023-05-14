#imports 
import os 
import gdown
import requests

#download mcule file
def download_mcule_molecules():

    link = 'https://mcule.s3.amazonaws.com/database/mcule_purchasable_in_stock_230514.smi.gz'

    #if molecules file already exists, do not download it again
    if os.path.exists('data/molecules.smi'):
        print('Molecules file already exists')
        return

    print('Downloading mcule purchasable molecules from mcule database')

    #download file from the internet using previous link
    r = requests.get(link, allow_redirects=True)

    #save file in the data folder
    open('data/mcule_purchasable_in_stock_230514.smi.gz', 'wb').write(r.content)

    #unzip file
    os.system('gunzip data/mcule_purchasable_in_stock_230514.smi.gz')

    # #remove zip file
    # os.system('rm data/mcule_purchasable_in_stock_230514.smi.gz')

    #rename file
    os.system('mv data/mcule_purchasable_in_stock_230514.smi data/molecules.smi')


#download molecular transformer file
def download_model_uspto480k():

    # Model trained on USPTO-480k
    trained_model_url = 'https://drive.google.com/uc?id=1ywJCJHunoPTB5wr6KdZ8aLv7tMFMBHNy'

    models_dir = os.path.join(os.path.dirname(__file__), 'models')
    model_name = "USPTO480k_model_step_400000.pt"   

    target_path = os.path.join(models_dir, model_name)

    if not os.path.exists(target_path):
        os.makedirs(models_dir, exist_ok=True)
        gdown.download(trained_model_url, target_path, quiet=False)
    else:
        print(f"{target_path} already exists")



if __name__ == "__main__":
    download_mcule_molecules()
    download_model_uspto480k()