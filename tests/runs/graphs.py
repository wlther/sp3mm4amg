#!/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
sns.set_palette("colorblind")

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

grouped = df.groupby(['Collection', 'Size', 'Smoothing', 'Threads', 'Implementation', 'Level'])

data = grouped[['Time-1', 'Time-2']].mean().reset_index()

data['Time-Combined'] = data['Time-1'] + data['Time-2']

for _,row in df[['Collection', 'Size', 'Smoothing']].drop_duplicates().iterrows():
    collection = row['Collection']
    size = row['Size']
    smoothing = row['Smoothing']
    
    subset_data = data[(data["Collection"] == collection) & (data["Size"] == size) & (data["Smoothing"] == smoothing)]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 12))  # Create three subplots
    
    sns.lineplot(
        x="Threads",
        y="Time-1",
        hue="Implementation",
        style="Level",
        data=subset_data,
        ax=axes[0],  # the left subplot
    )
    axes[0].set_title(f"Time-1\nCollection: {collection}, Size: {size}, Smoothing: {smoothing}")
    axes[0].set_xlabel("Threads")
    axes[0].set_ylabel("Time-1")
    
    sns.lineplot(
        x="Threads",
        y="Time-2",
        hue="Implementation",
        style="Level",
        data=subset_data,
        ax=axes[1],  # middle subplot
    )
    axes[1].set_title(f"Time-2\nCollection: {collection}, Size: {size}, Smoothing: {smoothing}")
    axes[1].set_xlabel("Threads")
    axes[1].set_ylabel("Time-2")
    
    sns.lineplot(
        x="Threads",
        y="Time-Combined",
        hue="Implementation",
        style="Level",
        data=subset_data,
        ax=axes[2],  # right subplot
    )
    axes[2].set_title(f"Combined Time\nCollection: {collection}, Size: {size}, Smoothing: {smoothing}")
    axes[2].set_xlabel("Threads")
    axes[2].set_ylabel("Combined Time")
    
    plt.tight_layout()
    
    pdf_filename = f"results/comparison_plot_{collection}_{size}_{smoothing}.pdf"
    plt.savefig(pdf_filename)
    
    plt.close()

print(data.sort_values(by=['Time-Combined']))