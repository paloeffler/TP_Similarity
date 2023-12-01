### TP similarity
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, MACCSkeys, rdMolDescriptors, rdFingerprintGenerator, rdMolHash, Fingerprints, rdDescriptors, Lipinski
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.AtomPairs import Pairs

### import function 
def calculate_similarity(x, y, fptype='Morgan', fpsize=2048):
    if fptype == 'RDKit':
        fp_generator = lambda mol: Chem.RDKFingerprint(mol, maxPath=7, fpSize=fpsize) # Use the default RDKit fingerprint
    elif fptype == 'Morgan':
        radius = 2 # Specify the radius for the Morgan fingerprint
        fp_generator = lambda mol: AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=fpsize)
    elif fptype == 'MACCSKeys':
        fp_generator = MACCSkeys.GenMACCSKeys
    elif fptype == 'Daylight-like':
        fp_generator = lambda mol: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
    elif fptype == 'Topological':
        targetSize = fpsize  # Number of torsions in the fingerprint
        fp_generator = lambda mol: rdMolDescriptors.GetTopologicalTorsionFingerprint(mol, targetSize=targetSize)
    elif fptype == 'Avalon':
        fp_generator = lambda mol: pyAvalonTools.GetAvalonFP(mol)
    else:
        raise ValueError(f"Fingerprint type {fptype} is not supported.")
    
    parent_fp = [fp_generator(mol) for mol in x]
    tp_fp = [fp_generator(mol) for mol in y]
    similar_fp = []
    
    for parent_fp_single, tp_fp_single in zip(parent_fp, tp_fp):
        similarity = DataStructs.TanimotoSimilarity(parent_fp_single, tp_fp_single)
        similar_fp.append(similarity)
    
    similar_count = len([sim for sim in similar_fp if sim > 0.95])
    percentage = round((similar_count / len(x)) * 100, 2)
    print(f'Percentage of parent-TPs that are similar using {fptype} is: {percentage}%')
    
### prepare data
data = pd.read_csv('/mnt/c/Users/pllo0001/OneDrive - Sveriges lantbruksuniversitet/Skrivbordet/PhD/Review/Viewpoint/similarity_compounds.csv')
data.head()
parent = data["SMILES1"]
tp = data["SMILES2"]
canon_parent_smiles = [Chem.CanonSmiles(i) for i in parent]
canon_tp_smiles = [Chem.CanonSmiles(i) for i in tp]
print(f'There are {len(set(canon_parent_smiles))} unique parents, {len(set(canon_tp_smiles))} unique tps (RDKit canonical SMILES).')
print(f'There are {len(set(data.Orcid1))} unique parents, {len(set(data.Orcid2))} unique tps(CIDs)')
parent_mols = [Chem.MolFromSmiles(s) for s in parent]
tp_mols = [Chem.MolFromSmiles(s) for s in tp]

calculate_similarity(parent_mols, tp_mols, fptype='RDKit', fpsize=2048)    
calculate_similarity(parent_mols, tp_mols, fptype='Morgan', fpsize=2048)
calculate_similarity(parent_mols, tp_mols, fptype='MACCSKeys', fpsize=2048)
calculate_similarity(parent_mols, tp_mols, fptype='Daylight-like', fpsize=2048)
calculate_similarity(parent_mols, tp_mols, fptype='Topological', fpsize=7)
calculate_similarity(parent_mols, tp_mols, fptype='Avalon', fpsize=2048)