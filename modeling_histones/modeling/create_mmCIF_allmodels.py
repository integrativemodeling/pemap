################################
# Modeling of H3/H4
# dimer using pE-MAP data
#
# iecheverria - Salilab - UCSF
# ignacia@salilab.org
################################
import sys
sys.path.insert(0,'/Users/iecheverria/opt/anaconda3/envs/imp/lib/python3.7/site-packages')
sys.path.insert(1,'/Users/iecheverria/SOFTW/imp/build/lib')
import ihm
print(ihm.__file__)
import ihm.reference

import IMP
exit()
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
import pprint
pprint.pprint(sys.path)

#print(ihm.__file__)
exit()

pprint.pprint(sys.path)
exit()
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
import ihm.reference

sys.path.append('../../utils/')
import make_archive

class pEMap_restraints(object):
    ''' Add MIC restraints to mmCIF files as upper distance bound '''
    

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
                                                              mic_value=mic_value)
                all_rest.append(rest)

        return all_rest


########################

models = {1: {'MIC_file': '../data/hbhistmic_180904_impute_pcsort975_5_4th_uniq_c03.txt',
              'analysis_dir': '../analysis_100/',
              'n_models':2500000,
              'n_cluster':20000,
              'model_precision':1.04,
              'dcd_file':'',
              'description': 'Full data set'}}
          #2: {'MIC_file': '../data/hbhistmic_180904_impute_pcsort975_5_4th_uniq_03_80.csv',
          #    'analysis_dir': '../analysis_80_1/' ,
          #    'n_models':2500000,
          #    'n_cluster':10000,
          #    'model_precision':1.2,
          #    'dcd_file':'',
          #    'description': '80% of the dataset, set 1'}}


sys.argv = ['', '../data/hbhistmic_180904_impute_pcsort975_5_4th_uniq_c03.txt', '--mmcif']

exec(open('mod_pemap_histones.py').read())
print('print', po.asym_units)


#po.system.title = ('Genetic interaction mapping informs integrative determination of biomolecular assembly structures')
for m, info in models.items():
    #bs_sys.system.add_protocol_output(po)

    # Add database with MIC values
    fname = info['MIC_file']
    D_dump = ihm.dumper._DatasetDumper()
    l = ihm.location.InputFileLocation(fname)
    D = ihm.dataset.Dataset(l, details=f'pE-MAP derived spatial restraints, MIC values for: {info["description"]}')
    D.data_type = 'Other'
    po.system.orphan_datasets.append(D)
    D_dump.finalize(po.system)

    # Add restraints
    PR = pEMap_restraints(po.asym_units)
    all_rest = PR.add_restraints(fname)

    rg = ihm.restraint.RestraintGroup(all_rest)
    po.system.restraint_groups.append(rg)

    # Correct number of output models to account for multiple runs
    protocol = po.system.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = info['n_models']

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered the 200000 frames down to one cluster of
    # N models:
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=200000,
                            num_models_end=9999))
    

    # Create an ensemble for the cluster
    e = po._add_simple_ensemble(analysis.steps[-1],
                                name="Cluster 0", num_models=9999,
                                drmsd=8.3, num_models_deposited=1,
                                localization_densities={}, ensemble_file=None)
    
    # Add the model from RMF
    rh = RMF.open_rmf_file_read_only(f'{info["analysis_dir"]}/clustering/cluster.0/cluster_center_model.rmf3')
    IMP.rmf.link_hierarchies(rh, [root_hier])
    IMP.rmf.load_frame(rh, RMF.FrameID(0))
    del rh
    model = po.add_model(e.model_group)

    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        name = asym.split('.')[0]
        fname = f'{info["analysis_dir"]}/clustering/cluster.0/LPD_{name}.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)
    
    '''
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=info['n_models'],
                            num_models_end=info['n_cluster']))
    
   
    
    e = po.add_clustering_ensemble(analysis.steps[-1],
                                   name=f"Cluster 0, {info['description']}, ", num_models=info['n_cluster'],
                                   precision=info['model_precision'], num_models_deposited=1,
                                   localization_densities={}, ensemble_file=None)
    
    # Add the model from RMF
    rh = RMF.open_rmf_file_read_only(f'{info["analysis_dir"]}/clustering/cluster.0/cluster_center_model.rmf3')
    IMP.rmf.link_hierarchies(rh, [root_hier])
    IMP.rmf.load_frame(rh, RMF.FrameID(0))
    del rh
    model = po.add_model(e.model_group)

    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        name = asym.split('.')[0]
        fname = f'{info["analysis_dir"]}/clustering/cluster.0/LPD_{name}.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)
    '''
    # Add uniprot of proteins
    # name : (uniprot id, mutations, [[db_begin, db_end, entity_begin, entity_end]]
    Uniprot={'H3.0':  ('P61830',[],[]),
             'H4.0':  ('P02309',[],[])}

    for prot, (entry, sd, limits) in Uniprot.items():
        print(prot, entry, sd, limits)
        ref = ihm.reference.UniProtSequence.from_accession(entry)
        for seg in limits:
            ref.alignments.append(ihm.reference.Alignment(
                db_begin=seg[0], db_end=seg[1], entity_begin=seg[2], entity_end=seg[3], seq_dif=sd))
            
        po.asym_units[prot].entity.references.append(ref)
    

    # Replace local links with DOIs
    #repos = []
    #for subdir, zipname in make_archive.ARCHIVES.items():
    #    print('subdir', subdir)
    #    repos.append(ihm.location.Repository(
    #        doi="10.5281/zenodo.3836213", root="../%s" % subdir,
    #        url="https://zenodo.org/record/3836213/files/%s.zip" % zipname,
    #        top_directory=os.path.basename(subdir)))
    # 
    #po.system.update_locations_in_repositories(repos)

    
po.flush()
    
    
