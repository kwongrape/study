import os
import pymol
import requests
import subprocess

from collections import Counter
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from pymol import cmd
from modeller import *
from modeller.automodel import *
from pathlib import Path                    


def run_msa(fasta, HETATM=False, HETATM_template):

    '''
    example 
    fasta = Path('D:\HM\job\seqeunce.fasta')
    '''

    muscle_exe = Path('D:/HM/muscle5.1.win64.exe')
    msa_in = Path(fasta)
    msa_out = Path(fasta.parent / f'{fasta.stem}_msa.fasta')

    run = subprocess.run([muscle_exe, '-align', msa_in, '-output', msa_out], capture_output=True, text=True)
    print(run.stdout + run.stderr)
    
    with open(msa_out, 'r') as f:
        lines = f.readlines()
        
        msa_seq = []
        for idx, line in enumerate(lines):
            
            if line.find(':') != -1:
                front = line.split('00')[-1]
                msa_seq.append(front)
            
            else:
                msa_seq.append(line)
               
    msa_results = []
    msa = []    
    for line in msa_seq:
        if line.startswith('>'):
            if msa:  
                msa_results.append(msa)
            msa = [line.strip()]  
        else:
            msa.append(line.strip())

    if msa:
        msa_results.append(msa)

    msa_pir = []
    with open(fasta, 'r') as f:
        lines = f.readlines()
        
        titles = []
        for idx, line in enumerate(lines):
            
            if line.startswith('>'):
                titles.append(lines[idx:idx + 2])

    msa_pir = []
    sorting_title = sorted(titles)            
    sorting_msa = sorted(msa_results)
    for idx in range(len(msa_results)):
        msa_pir.append(sorting_title[idx]+sorting_msa[idx][1:])

    with open(fasta.parent / 'processed_msa.ali', 'w') as f:
        
        if HETATM:
            pirs = []
            for msa in msa_pir:
                
                pir = ''.join(msa) + '\n'
                if pir.find(HETATM_template) != -1 or pir.find('target') != -1:
                    pir
                    pir = pir.replace('*', '/..*')
                    
                else:
                    pir = pir.replace('*', '/--*')

                pirs.append(pir)
                
        else:
            pirs = []
            for msa in msa_pir:
                pirs.append(''.join(msa) + '\n')
        
        f.write(''.join(pirs))
                   
def to_fasta(path, PDB, target):
    '''
    to_fasta(Path('D:/HM/job'), ['4jxf_A'], 4jxf_A) 
    '''
    os.chdir(path)

    for ID in PDB:
        pymol.cmd.fetch(ID)
        pymol.cmd.remove('not polymer')
        if target == ID:
            
            if multi_chain:
                pymol.cmd.get_fastastr(ID).split(ID)
            else:
                cif_seq = ''.join(pymol.cmd.get_fastastr(ID).split('\n')[1:])
        pymol.cmd.save(ID+'.fasta', 'polymer')
        pymol.cmd.save(ID+'.pdb', 'polymer')
        pymol.cmd.delete('all')
    
    env = Environ()
    aln = Alignment(env)
    for ID in PDB:
        m = Model(env, file=ID)
        aln.append_model(m, align_codes=ID)
    
    # ~ Target sequence
    if cif_seq:
        aln.append_sequence(cif_seq)
        aln[len(PDB)].code='target'
        aln.write(file='seqeunce.fasta')
    
    multi_chain = True 
    if multi_chain:
        cif_seq.split('>')

uniprot_ID = 'Q9NWZ3'
def uniprot_seq(uniprot_ID):
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_ID}.txt'
    
    try:
        response = requests.get(url)

        if response.status_code == 200:
            text = response.text
        
        else:
            print(f'failed to load {uniprot_ID}.txt')
        
    except requests.exceptions.RequestException as e:
         print(f'Check ur {uniprot_ID}')
    
    seq = []
    domain = []
    switch = False
    for line in text.split('\n'):
        
        if switch:

            if line.startswith('//'):
                switch = False
            else:
                seq.append(line.replace(' ', ''))  
            
        if line.startswith('SQ'):
            length = line.split()[2]
            switch = True
        
        if line.startswith('DR   PDB; '):
            info = line.split()[-1].split('=')
            # k_chain = info[0]
            domain.append(info[1])
                        
    seq_str = ''.join(seq)
    idx = Counter(domain).most_common(1)[0][0].replace('.', '').split('-')
    recommand_seq = seq_str[int(idx[0])-1:int(idx[1])-1]
    
    return recommand_seq

if __name__ == '__main__':
    pymol._launch_no_gui()
    PDB = ['6qas']
    target = '6qas'
    #! single template
    to_fasta(Path('D:/HM/job'), ['4jxf_A'], '4jxf_A') 
    
    #! multi template
    to_fasta(Path('D:/HM/job'), ['4jxf_A', '3cok_A', '4yur_A'], '4jxf_A') 
  
    #! run msa
    run_msa(Path('D:/HM/job/seqeunce.fasta'))
