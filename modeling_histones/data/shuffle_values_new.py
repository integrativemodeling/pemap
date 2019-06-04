################################################
# Shuffle values for all dataset (no cutoff)
# and then get based on the 0.3 
#################################################

import pandas as pd
import sys
import math
import random
import numpy as np

def shuffle_column(df, column):
    df_new = df.copy()
    #print('----', df_new['mic'])
    df_new[column] = list(np.random.permutation(df_new[column].values))
	#df_new[column].reindex(df_new.loc[column].values) 
	#random.sample(list(df_new.loc[column].values), len(df_new.loc[column].values))
    return df_new

file = sys.argv[1]

DF = pd.read_csv(file, sep=' ', header=None, names=['p1','p2','r1','r2','mic','d'])
print(DF.head())

for i in range(3):
	DF_shuffled = shuffle_column(DF,'mic')
	DF_shuffled_sel = DF_shuffled[(DF_shuffled['mic']>=0.3)]
	print(DF_shuffled.head())
	print(len(DF), len(DF_shuffled), len(DF_shuffled_sel))

	file_out = file.split('.')[0]+'_shuffled_all_%s.csv'%(i+1)
	DF_shuffled_sel.to_csv(file_out, sep=' ', index=False, header=True)

