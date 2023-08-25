#!/bin/python3

import pandas as pd

csv_file = 'results/results.csv'

dtypes = {'Implementation' : 'category',
         'Collection' : 'category',
         'Size': 'category',
         'Smoothing' : 'category',
         'Level': 'category',
         'Preparation time': 'float64',
         'Threads': 'category',
         'N-1': 'int64',
         'M-1': 'int64',
         'NNZ-1': 'int64',
         'Time-1': 'float64',
         'N-2': 'int64',
         'M-2': 'int64',
         'NNZ-2': 'int64',
         'Time-2': 'float64'
         }

df = pd.read_csv(csv_file, dtype=dtypes)

# Calculate a combined performance score for each row (implementation) based on Time-1 and Time-2
df['Combined-Time'] = df['Time-1'] + df['Time-2']

for i in df.sort_values(by=['Combined-Time']).groupby(['Collection', 'Size', 'Smoothing']):
    print(i)