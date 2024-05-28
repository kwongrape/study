#!/bin/bash

read -p "input file name : "  input
read -p "tell me name of ID col : " id
read -p "number of core : "  core


path=`pwd`

prepare_2D=$(python << end
import pandas as pd
import numpy as np
import multiprocessing as mp
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools, SDMolSupplier

def canonicalize_smiles(smiles):
    return Chem.CanonSmiles(max(smiles.split('.'), key=len))

def parallel_apply(df, func, num_processes):
    df_split = np.array_split(df, num_processes)
    pool = mp.Pool(num_processes)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

def canonicalization_wrapper(df):
    df['Smiles'] = df['Smiles'].apply(canonicalize_smiles)
    return df

DB = pd.read_csv(f'$path/$input', sep=',')
DB = DB[:3000]
num_cores = mp.cpu_count()
DB = parallel_apply(DB, canonicalization_wrapper, num_cores)

cols = [col for col in DB.columns]
PandasTools.AddMoleculeColumnToFrame(DB)
PandasTools.WriteSDF(DB, out=str('$path/2D.sdf'), properties=cols)
end
)

chunk=$(python << chunk_end
import time
with open(f'$path/2D.sdf', 'r') as f:
    lines = f.readlines()
    
    r, rs = [], {}
    for i, l in enumerate(lines):
        r.append(l)
        if l.find('<$id>') != -1:
            ID = lines[i+1].replace('\n', '')
        if l.startswith('\$\$\$\$'):
            rs[ID] = r
            r, ID = [], []
    core = 18
    chunk_point = len(rs) // int(core)
    end_point = len(rs) // chunk_point
    end_point2 = len(rs) % chunk_point
	
    n, m = 0, 0
    rows = []    
    for k, v in rs.items():
        n += 1
        rows.append(''.join(v))
        if n == int(chunk_point):
            n = 0
            with open(f'$path/chunk{m}.sdf', 'w') as f:
                f.write(''.join(rows))
            rows = []
            m += 1
        elif (m == end_point) & (n == end_point2):
            with open(f'$path/chunk{m}.sdf', 'w') as f:
                f.write(''.join(rows)) 
chunk_end
)

file_lst=`ls | grep chunk`
for file in $file_lst; do
	if [ ! -d "${file%.sdf}" ]; then
		mkdir "${file%.sdf}"
	fi
done

split=$(python <<split_end
import os
chunk_lst = [f for f in os.listdir('$path') if f.startswith('chunk') and os.path.isdir(os.path.join('$path', f))]
for cid in chunk_lst:
    with open(f'$path/{cid}.sdf', 'r') as f:
        lines = f.readlines()

        n = 0
        r, rs = [], {}
        for i, l in enumerate(lines):
            r.append(l)
            if l.find('<$id>') != -1:
                ID = lines[i+1].replace('\n', '')
            if l.startswith('\$\$\$\$'):
                rs[ID] = r
                r, ID = [], []
        
        for k, v in rs.items():
            with open(f'$path/{cid}/{k}.sdf', 'w') as f:
                f.write(''.join(v))
split_end
)    

for file in $file_lst; do
	obabel -isdf ./"${file%.sdf}"/*.sdf -opdb -O ./"${file%.sdf}"/3D_*.pdb --gen3D -pH 7.4 -P --partialcharge gasteiger &
done

wait

if [ ! -d "fusion" ]; then
	mkdir "fusion"
fi

file_lst=`ls | grep chunk`
for file in $file_lst; do
	rm ./"${file%.sdf}"/*.sdf
	mv ./"${file%.sdf}"/*.pdb ./fusion/ &
done

#TODO will be fix
cd fusion
ligand_lst=`ls | grep *.pdb`
for ligand in $ligand_lst; do
	prepare_ligand.py $ligand --output $ligand.pdbqt
done
