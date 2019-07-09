from modeller import *
from modeller.automodel import *

# Override the 'special_restraints' and 'user_after_single_model' methods:

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
    
    def special_patches(self, aln):
        self.rename_segments(segment_ids=('B'), renumber_residues=[24])

env = environ()

env.io.hetatm = True

a = MyModel(env, alnfile='aln_rpb2_4ayb.pir', knowns=('4ayb_B_t'),sequence='mod_rpb2_4ayb',assess_methods=(assess.DOPE,assess.normalized_dope))

a.starting_model = 1
a.ending_model = 50
a.make()                           # do comparative modeling


