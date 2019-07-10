###############################
# Modeling of Rpb1/Rpb2
# dimer using pE-MAP data
#
# iecheverria - Salilab - UCSF
# ignacia@salilab.org
###############################
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.pemap
import IMP.pmi.mmcif

import numpy as np
import sys

top_dir = './'
if len(sys.argv)<2:
    print('USAGE: python mod_pemap_rnap.py mic_file mic_weight')
    exit()

mic_file = sys.argv[1]
w_pemap = float(sys.argv[2])
print('Using file ...', mic_file)

###################### SYSTEM SETUP #####################
mdl = IMP.Model()

topo = top_dir+'topology_rnap.dat'
reader_sys = IMP.pmi.topology.TopologyReader(topo,
                                             pdb_dir = top_dir+'data/',
                                             fasta_dir = top_dir+'data/')

bs_sys = IMP.pmi.macros.BuildSystem(mdl,
                                      resolutions=[1])

##############################
# Generate mmcif file
##############################

if '--mmcif' in sys.argv:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput(open('RNAP_pEMAP_initial.cif', 'w'))
    bs_sys.system.add_protocol_output(po)
    po.system.title = ('Genetic interaction mapping informs integrative structure determination of molecular assemblies')
    # Add publication
    #po.system.citations.append(ihm.Citation.from_pubmed_id(0000))

##############################
# Build state
##############################
    
bs_sys.add_state(reader_sys)

root_hier,  dof = bs_sys.execute_macro(max_rb_trans=2.0,
                                       max_rb_rot=0.2)
mols = bs_sys.get_molecules()[0]
##############################
# Connectivity
##############################
output_objects = [] # keep a list of functions that need to be reported
sample_objects = []
rmf_restraints = []

crs = []
for molname in mols:
    for mol in mols[molname]:
        copy_n = IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.set_label(mol.get_name()+'.'+str(copy_n))
        cr.add_to_model()
        output_objects.append(cr)
        crs.append(cr)

##############################
# Excluded Volume
##############################
evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols.values(),
                                                               resolution=1)
evr1.add_to_model()
evr1.set_weight(1.0)
output_objects.append(evr1)

##############################
# pE-map restraint
##############################
pemap = IMP.pmi.restraints.pemap.pEMapRestraint(root_hier,
                                                top_dir+mic_file,
                                                sigma_init=5.0)


pemap.add_to_model()
pemap.set_weight(w_pemap)
output_objects.append(pemap)
print('output pemap: ', pemap.get_output())
dof.get_nuisances_from_restraint(pemap)

##############################
# COM distance restraint
##############################

dCOM = IMP.pmi.restraints.pemap.COMDistanceRestraint(root_hier,
                                            protein0 = 'rpb1',
                                            protein1 = 'rpb2',
                                            distance = 150.0,
                                            strength = 80.)
                                      
    
dCOM.add_to_model()
dCOM.set_weight(1.0)
output_objects.append(dCOM)
print('output dCOM: ', dCOM.get_output())

##############################    
# Write initial coordinates
##############################
out = IMP.pmi.output.Output()
out.init_rmf("all_ini.rmf3", [root_hier])
out.write_rmf("all_ini.rmf3")

##############################
# Shuffle
##############################    
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    bounding_box=((-80,-80,-80),(80,80,80)))

dof.optimize_flexible_beads(500)

##############################
# MC Sampling
##############################    

rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=root_hier,                          
                                    crosslink_restraints=rmf_restraints,          
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=5.0,
                                    global_output_directory="output/",
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=60000,
                                    number_of_best_scoring_models=0)

rex.execute_macro()
po.flush()
exit()
