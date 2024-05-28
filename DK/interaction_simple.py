import ray
import numpy as np

from pathlib import Path
from tqdm import tqdm
from scipy.spatial.distance import cdist

ray.init()

@ray.remote
def process_vina_interaction(file, output, hotspot, cutoff=3.5, type='H', save=False):    
    with open(file, 'r') as f:
        lines = f.readlines()

    rank = {}
    pdbqt = []
    for idx, line in enumerate(lines):
        pdbqt.append(line)
        
        if line.startswith('MODEL '):
            name = line.split()[-1] + '_' + lines[idx + 1].split()[3]
        
        if line.startswith('ENDMDL'):
            
            atom_infos = []
            for atom in pdbqt:
                if atom.startswith('ATOM'):
                    atom_infos.append(np.array(atom.split())[[2, 5, 6, 7]])

            atom_infos_np = np.array(atom_infos)
            HBD = atom_infos_np[np.isin(atom_infos_np[:, 0], ['H'])][:, 1:]
            HBA = atom_infos_np[np.isin(atom_infos_np[:, 0], ['O', 'N'])][:, 1:]
            
            if type == 'H':
                matrix = cdist(HBD.astype(float), np.array(hotspot).astype(float), metric='euclidean')
            elif type == 'A':
                matrix = cdist(HBA.astype(float), np.array(hotspot).astype(float), metric='euclidean')
            
            if np.any(matrix <= cutoff):
                rank[name] = pdbqt
            
            pdbqt = []

    if save and rank:
        output_path = Path(output)
        output_path.mkdir(parents=True, exist_ok=True)
        for k, v in rank.items():
            with open(output_path / f'{Path(file).stem}_{k}_.pdbqt', 'w') as f:
                f.write(''.join(v))

hotspot = [[1, 2, 3]]
cutoff = 350
type = 'H'
save = True

results = []
for i in tqdm(range(10000)):
    result = process_vina_interaction.remote(file, output, hotspot, cutoff, type, save)

# ! Ori code
'''
file = Path(r'D:\HM\drugbank1_out.pdbqt')
hotspot = [['0.0', '0.0', '0.0'], ['5.0', '5.0', '5.0']]
output = file.parent / 'job'

for i in tqdm(range(10000)):
    vina_interaction(file, output, hotspot, type='H', save=True)

def vina_interaction(file, output, hotspot, cutoff=3.5, type='H', save=False):    
    with open(file, 'r') as f:
        lines = f.readlines()

        rank = {}
        pdbqt = []
        for idx, line in enumerate(lines):
            pdbqt.append(line)
            
            if line.startswith('MODEL '):
                name = line.split()[-1] + '_' + lines[idx + 1].split()[3]
            
            if line.startswith('ENDMDL'):
                
                atom_infos = []
                for atom in pdbqt:
                    
                    if atom.startswith('ATOM'):
                        atom_infos.append(np.array(atom.split())[[2,5,6,7]])

                HBD = np.array(atom_infos)[np.isin(np.array(atom_infos)[:, 0], ['H'])][:, 1:]
                HBA = np.array(atom_infos)[np.isin(np.array(atom_infos)[:, 0], ['O', 'N'])][:, 1:]
                
                if type == 'H':
                    matrix = cdist(HBD.astype(float), np.array(hotspot).astype(float), metric='euclidean')
                elif type == 'A':
                    matrix = cdist(HBA.astype(float), np.array(hotspot).astype(float), metric='euclidean')
                
                if np.any(matrix <= cutoff):
                    rank[name] = pdbqt
                
                pdbqt = []
    
    if save:
        
        if rank:
            for k, v in rank.items():
                
                with open(output / f'{file.stem}_{k}.pdbqt', 'w') as f:
                    f.write(''.join(v))
'''