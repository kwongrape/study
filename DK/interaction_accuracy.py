import pandas as pd

from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport
from pathlib import Path
from tqdm import tqdm


def profiling(pdb):
    interaction_types = ("hydrophobic", "hbond", "waterbridge", "saltbridge", "pistacking" "pication", "halogen", "metal",)
    structure = PDBComplex()
    structure.load_pdb(pdb)
    
    for ligand in structure.ligands:
        structure.characterize_complex(ligand)
        
    sites = {}
    for key, site in sorted(structure.interaction_sets.items()):
        binding_site = BindingSiteReport(site)
        
    interactions = {k: [getattr(binding_site, k + "_features")] + getattr(binding_site, k + "_info")for k in interaction_types}
    sites[key] = interactions
    
    selected_site_interactions = list(sites.key())[0]
    ligand = sites[selected_site_interactions]
    
    interaction_lst = pd.DataFrame()
    for inter in interaction_types:
        interactions = pd.DataFrame.from_records(ligand[inter][1:], columns=ligand[inter][0],)
        interactions['Type'] = inter
        interaction_lst = pd.concat([interaction_lst, interactions])
    
    return interaction_lst

pdb = Path(r'D:\HM\selected_compound_PLK4_23.pdb')

for i in tqdm(range(1000)):
    interactions_by_site = prepare(str(pdb))

stack.to_excel('D:/HM/test.xlsx')

