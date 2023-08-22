#!/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

csv_file = 'results.csv'

dtypes = {'Implementation' : 'category',
         'R' : 'category',
         'AC': 'category',
         'P' : 'category',
         'AC_NEXT': 'category',
         'Preparation time': 'float64',
         'Threads': 'category',
         'N-RAC': 'int64',
         'M-RAC': 'int64',
         'NNZ-RAC': 'int64',
         'Time-RAC': 'float64',
         'N-RACP': 'int64',
         'M-RACP': 'int64',
         'NNZ-RACP': 'int64',
         'Time-RACP': 'float64'
         }

df = pd.read_csv(csv_file, dtype=dtypes)

grouped = df.groupby(['R', 'Threads', 'Implementation'])

# Name of R should be in the following shape :
# data/<collection>/<size>/<smoothing>/dump_lev_d_p0_l<level>_r.mtx
collections_and_sizes = list(set(['/'.join(s.split('/')[1:3]) for s in df['R'].unique()]))
print(collections_and_sizes)

# sizes = list(set([s.split('/')[2] for s in df['R'].unique()]))

# data = sns.load_dataset()
data = grouped['Time-RAC'].mean().reset_index()

# sns.relplot(data=data, x='Threads', y='Time-RAC', hue='Implementation')
