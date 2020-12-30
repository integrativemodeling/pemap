################################
# Modeling of H3/H4
# dimer using pE-MAP data
#
# iecheverria - Salilab - UCSF
# ignacia@salilab.org
################################
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

####################### INPUT FILES ######################
print(len(sys.argv))
if len(sys.argv) < 2:
    print(len(sys.argv))
    print('USAGE: python mod_pemap_histones.py ../data/hbhistmic_180904_impute_pcsort975_5_4th_uniq_c03.txt')
    exit()
    
mic_file = sys.argv[1]
print('Using file ...', mic_file)
###################### SYSTEM SETUP ######################

mdl = IMP.Model()
topo = 'top_his_comp_models.dat'
reader_sys = IMP.pmi.topology.TopologyReader(topo,
                                             pdb_dir = '../data/',
                                             fasta_dir = '../data/')

bs_sys = IMP.pmi.macros.BuildSystem(mdl,
                                      resolutions=[1])

if '--dry-run' is sys.argv:
    bs_sys.dry_run = True
else:
    bs_sys.dry_run = False
    

##############################
# Generate mmcif file
##############################

if '--mmcif' in sys.argv:
    print('hola')
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput(open('h3_h4_pEMAP_all.cif', 'w'))
    po.system.title = ('Genetic interaction mapping informs integrative determination of biomolecular assembly structures')
    bs_sys.system.add_protocol_output(po)
    
    # Add publication
    #po.system.citations.append(ihm.Citation.from_pubmed_id(0000))

##############################
# Build state
##############################
    
bs_sys.add_state(reader_sys)

root_hier,  dof = bs_sys.execute_macro(max_rb_trans=0.3,
                                       max_rb_rot=0.1)
mols = bs_sys.get_molecules()[0]
states = IMP.atom.get_by_type(root_hier,IMP.atom.STATE_TYPE)

# Write coordinates
out = IMP.pmi.output.Output()
out.init_rmf("all_ini.rmf3", [root_hier])
out.write_rmf("all_ini.rmf3")
out.close_rmf("all_ini.rmf3")


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

pemap = IMP.pmi.restraints.pemap.PEMAPRestraint(root_hier,
                                                mic_file,
                                                sigma_init=5.0)

    
pemap.add_to_model()
pemap.set_weight(1.0)
output_objects.append(pemap)
print('Output PEMAP restraint: ', pemap.get_output())
dof.get_nuisances_from_restraint(pemap)

##############################
# COM distance restraint
##############################
dCOM = IMP.pmi.restraints.pemap.COMDistanceRestraint(root_hier,
                                            'h3',
                                            'h4',
                                            distance = 20.0,
                                            strength = 30.)


dCOM.add_to_model()
dCOM.set_weight(1.0)
output_objects.append(dCOM)
print('Output dCOM restraint: ', dCOM.get_output())

##############################
# Shuffle
##############################
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    bounding_box=((-80,-80,-80),(80,80,80)))

dof.optimize_flexible_beads(200)
############################# SAMPLING ##############################
# Run replica exchange Monte Carlo sampling

num_frames = 50000
if '--mmcif' in sys.argv or '--test' in sys.argv:
    num_frames=5

print('Number of frames:', num_frames)

rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=root_hier,                          
                                    crosslink_restraints=rmf_restraints,          
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=3.0,
                                    global_output_directory="output/",
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=num_frames,
                                    number_of_best_scoring_models=0)

rex.execute_macro()
############################# mmCIF #################################

class pEMap_restraints(object):

    def __init__(self, asym):
        self.asym = asym

        self.residues = {}

        self.prots = {p.split('.')[0]:po.asym_units[p] for p in po.asym_units}
        
    def upper_bound(self, MIC):
        import math

        k = -0.014744
        n = -0.41 
        
        if MIC <= 0.6:
            return (math.log(MIC) - n) / k
        else: return (math.log(0.6) - n) / k
        

    def add_restraints(self, fname):
        
        """ Parse the MIC scores text file and return a set of restraints"""
        l = ihm.location.InputFileLocation(fname)
        with open(fname) as fh:
            return self.get_pEMAP_restraints(fh)
            
    def get_pEMAP_restraints(self, fh):
        all_rest = []
        for line in fh:
            if line == '\n': continue
            resids = line.split()
            if len(resids) >= 5:
                mic_distance = self.upper_bound(float(resids[4]))
                dist = ihm.restraint.UpperBoundDistanceRestraint(mic_distance)

                r0 = self.prots[resids[0]](int(resids[2]),int(resids[2]))
                r1 = self.prots[resids[1]](int(resids[3]),int(resids[3]))
                
                # Check if already in the DB
                if r0 in self.residues.keys():
                    f0 = self.residues[r0]
                else:
                    f0 = ihm.restraint.ResidueFeature([r0])
                    self.residues[r0] = f0

                if r1 in self.residues.keys():
                    f1 = self.residues[r1]
                else:
                    f1 = ihm.restraint.ResidueFeature([r1])
                    self.residues[r1] = f1
                
            
                rest = ihm.restraint.DerivedDistanceRestraint(dataset=D,
                                                              feature1=f0, feature2=f1, distance=dist,
                                                              probability=1.0)
                all_rest.append(rest)

        return all_rest

A = 0
if A == 100:
#if '--mmcif' in sys.argv:
    #import ihm.cross_linkers
    import ihm.dumper
    import ihm.format
    import ihm.location
    import ihm.representation
    import ihm.startmodel
    import ihm.dataset
    import ihm.protocol
    import ihm.analysis
    import ihm.model
    import ihm.restraint
    import ihm.geometry
    
    fname = '../data/hbhistmic_180904_impute_pcsort975_5_4th_uniq_c03.txt'
        
    D_dump = ihm.dumper._DatasetDumper()
    l = ihm.location.InputFileLocation(fname)
    D = ihm.dataset.Dataset(l)
    po.system.orphan_datasets.append(D)
    D_dump.finalize(po.system)

    PR = pEMap_restraints(po.asym_units)
    all_rest = PR.add_restraints(fname)

    rg = ihm.restraint.RestraintGroup(all_rest)
    po.system.restraint_groups.append(rg)

    
    # Correct number of output models to account for multiple runs
    protocol = po.system.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = 500000

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered the 200000 frames down to one cluster of
    # 100 models:
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=500000,
                            num_models_end=100))
    
    # Create an ensemble for the cluster (warning: _add_simple_ensemble
    # is subject to change in future IMP releases) and deposit a single
    # representative model (let's say it's frame 42 from the output RMF file)
    e = po._add_simple_ensemble(analysis.steps[-1],
                                name="Cluster 0", num_models=20000,
                                drmsd=1.04, num_models_deposited=1,
                                localization_densities={}, ensemble_file=None)
    
    # Add the model from RMF
    rh = RMF.open_rmf_file_read_only('../analysis/clustering/cluster.0/cluster_center_model.rmf3')
    IMP.rmf.link_hierarchies(rh, [root_hier])
    IMP.rmf.load_frame(rh, RMF.FrameID(0))
    del rh
    model = po.add_model(e.model_group)

    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        name = asym.split('.')[0]
        fname = f'../analysis/clustering/cluster.0/LPD_{name}.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)
        

    # Replace local links with DOIs
    #repo = ihm.location.Repository(doi="10.5281/zenodo.2598760", root="../..",
    #              top_directory="salilab-imp_deposition_tutorial-1ad5919",
    #              url="https://zenodo.org/record/2598760/files/salilab/"
    #                  "imp_deposition_tutorial-v0.2.zip")
    #po.system.update_locations_in_repositories([repo])

if '--mmcif' in sys.argv:
    po.flush()

