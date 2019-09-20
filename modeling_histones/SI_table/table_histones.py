###################################
# Script to summarize all modeling
# information into one table
#
# iecheverria - Salilab - UCSF
# ignacia@salilab.org
###################################

import pandas as pd
import glob
import os
import sys
import numpy as np

sys.path.append('../../pyext/src/analysis')
from create_summary_table import *
import utils

###########################
# Read files
###########################
modeling_script = '../modeling/mod_pemap_histones.py'
mmcif_file = '../modeling/h3_h4_pEMAP_initial.cif'

analysis_dir = '../analysis/'
clustering_dir = os.path.join(analysis_dir,'clustering')
rmf3 = os.path.join(clustering_dir,'cluster.0','cluster_center_model.rmf3')

I = get_input_information(mmcif_file)
input_information = I.get_dictionaries()

R = get_representation(clustering_dir)
representation = R.get_dictionaries()

S = read_modeling_information(modeling_script,
                              analysis_dir,
                              clustering_dir)

sampling = S.get_dictionaries_sampling()

samples = S.get_dictionaries_models()

clustering = S.get_dictionaries_clustering()

#S.update_mmcif_file(mmcif_file)

V = read_validation_information(clustering_dir)
validation = V.get_dictionaries()

V = read_benchmark_information(clustering_dir)
benchmark = V.get_dictionaries()

SS = get_software_information(mmcif_file)
software = SS.get_dictionaries()

D = get_data_availability(clustering_dir)
data_availability = D.get_dictionaries()

################################################
# Edit dictionaries
# Entries is dictionaries can be edited to add
# other custom information
################################################
input_information['Experimental data'] = ['170 pE-MAP derived distance restraints']

representation['Spatial restraints encoded into scoring function'].append('pE-MAP MIC pair-restraints; applied to the R1 representation')

sampling['CPU time'] = ['4 hours on 20 processors']

validation['Percent pE-MAP restraints satisfied per structure'] = ['87 \%']
validation['Percent of excluded volume restraints satisfied per structure'] = ['99 \%']

benchmark['Structural accuracy (95 \% CI)'] = ['3.8 (2.6-5.1) \AA']
benchmark['PDB used for benchmark'] = ['1ID3']


software['Modeling scripts'] = ['https://github.com/salilab/pemap']
software['Homology detection and structure prediction'] = ['HHPred, version 2.0.16']
software['Visualization and plotting'] = ['UCSF Chimera, version 1.10', 'Matplotlib, version 3.0.3 ']

################################################
# Convert ordered dictionaries 
# into lists
################################################
input_information_list = dict_to_list(input_information)
representation_list = dict_to_list(representation)
sampling_list = dict_to_list(sampling)
samples_list = dict_to_list(samples)
clustering_list = dict_to_list(clustering)
validation_list = dict_to_list(validation)
benchmark_list = dict_to_list(benchmark)
software_list = dict_to_list(software)
data_availability_list = dict_to_list(data_availability)


print(sampling_list)


################################################
# Compile all information
# 
################################################
variable_dict = {'complex': 'Histones H3-H4 dimer',
                 'number':1,
                 'input_information': input_information_list, 
                 'representation': representation_list,
                 'sampling': sampling_list,
                 'samples': samples_list,
                 'clustering':clustering_list,
                 'validation':validation_list,
                 'benchmark':benchmark_list,
                 'software':software_list,
                 'data':data_availability_list}

################################################
# Generate tex, pdf file
################################################
template = utils.get_template('../../pyext/src/analysis/SI_template.tex')
utils.compile_pdf_from_template(template, variable_dict, './table_SI_histones.pdf')

exit()
