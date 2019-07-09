from modeller import *
from modeller.automodel import *

# Override the 'special_restraints' and 'user_after_single_model' methods:

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
    
    def special_patches(self, aln):
        self.rename_segments(segment_ids=('A'), renumber_residues=[1])

env = environ()

env.io.hetatm = True

a = MyModel(env, alnfile='aln_rpb1_6gmh.pir', knowns=('6gmh_A_t'),sequence='mod_rpb1_6gmh',assess_methods=(assess.DOPE,assess.normalized_dope))

a.starting_model = 1
a.ending_model = 50
a.make()                           # do comparative modeling


