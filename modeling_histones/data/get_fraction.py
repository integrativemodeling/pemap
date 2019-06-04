import pandas as pd
import sys
import math
file = sys.argv[1]
fraction = float(sys.argv[2])

DF = pd.read_csv(file)
n_rows = len(DF)
n_sel = int(math.ceil(fraction*n_rows))
print('Reading file: ', file, 'selecting ', fraction)
print(DF.head(), n_rows)

DF_sel = DF.sample(n=n_sel)
print(len(DF_sel))

file_out = file.split('.')[0]+'_'+str(int(fraction*100))+'_3.csv'
DF_sel.to_csv(file_out, index=False, header=True)

