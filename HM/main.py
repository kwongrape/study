from utils import *

def TBM_single(ali, template):
    '''
    ali = Path('D:\HM\job\processed_msa.ali')
    template = '4jxf_A'
    '''
    env = Environ()
    env.io.atom_files_directory = [str(ali.parent)]
    
    run = AutoModel(env, alnfile=str(ali),
                    knowns=template, sequence='target',
                    assess_methods=(assess.DOPE))
    run.starting_model = 1
    run.ending_model = 5
    run.make()

def TBM_multimer(ali, template):
    '''
    ali = Path('D:\HM\job\processed_msa.ali')
    template = ('4jxf_A', '3cok_A', '4yur_A')
    '''
    env = Environ()
    env.io.atom_files_directory = [str(ali.parent)]
    
    run = AutoModel(env, alnfile=str(ali),
                    knowns=template, sequence='target')
    run.starting_model = 1
    run.ending_model = 5
    run.make()

ali = Path('D:/HM/job/test.ali')
def TBM_multichain(ali, template, sequence):
    '''
    need to handwork
    If homodimer need to restraint
    class MyModel(AutoModel):
    def special_restraints(self, aln):
        # Constrain the A and B chains to be identical (but only restrain
        # the C-alpha atoms, to reduce the number of interatomic distances
        # that need to be calculated):
        s1 = Selection(self.chains['A']).only_atom_types('CA')
        s2 = Selection(self.chains['B']).only_atom_types('CA')
        self.restraints.symmetry.append(Symmetry(s1, s2, 1.0))
    def user_after_single_model(self):
        # Report on symmetry violations greater than 1A after building
        # each model:
        self.restraints.symmetry.report(1.0)
    '''
    env = Environ()
    env.io.atom_files_directory = [str(ali.parent)]
    
    run = AutoModel(env, alnfile = str(ali),
                    knowns='2abx', sequence='1hc9')
    run.starting_model = 1
    run.ending_model = 5
    run.make()
    
    
def TBM_hetero(ali, template):
    '''
    ali = Path('D:\HM\job\processed_msa.ali')
    template = ('4jxf_A', '3cok_A', '4yur_A') 
    
    if restraints consider below code
    class MyModel(AutoModel):
    def special_restraints(self, aln):
        rsr = self.restraints
        for ids in (('NH1:161:A', 'O1A:336:B'),
                    ('NH2:161:A', 'O1B:336:B'),
                    ('NE2:186:A', 'O2:336:B')):
            atoms = [self.atoms[i] for i in ids]
            rsr.add(forms.UpperBound(group=physical.upper_distance,
                                     feature=features.Distance(*atoms),
                                     mean=3.5, stdev=0.1))   
                                     
    This method need to handwork
    
    1. template pdb file have to contain target HETATM
    2. In ali format file need to change end of residues counts and fix "/.*"
    ps "." mean = number of HETATM
    if 2 HETATM
        /.. necessary to 2 HETATM
        
    / = chain seperator
    . = HETATM indicator
    '''
    env = Environ()
    env.io.hetatm = True
    env.io.atom_files_directory = [str(ali.parent)]
    
    run = AutoModel(env, alnfile=str(ali),
                    knowns='4jxf_A', sequence='target')
    run.starting_model = 1
    run.ending_model = 5
    run.make()
      
def loop_opt(template, level, iters, loop_start, loop_end, chain):
    '''
    template = Path('D:\HM\job\processed_msa.ali')
    '''
    env = Environ()
    env.io.atom_files_directory = [str(template.parent)]

    class MyLoop(LoopModel):
        
        def select_loop_atoms(self):
        
            return Selection(self.residue_range(f'{loop_start}:{chain}', f'{loop_end}:{chain}'))

    m = MyLoop(env,
            inimodel=template.name,
            sequence=template.stem)       

    m.loop.starting_model= 1          
    m.loop.ending_model  = iters          

    if level == 'fast':
        m.loop.md_level = refine.very_fast 
    elif level == 'slow':
        m.loop.md_level = refine.very_slow 
    else:
        print('define refine level')

    m.make()
