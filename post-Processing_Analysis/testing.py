
import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import seaborn as sb
import scipy as scipy
from scipy import stats
import glob
import seaborn as sns

# Define constants
RUs_canonical = ["A", "G", "C", "T", "U"]
RUs_non_canonical = ["TAP-$C^5$", "TAP-N", "BA-$C^5$", "BA-N", "CA", "MM"]



# Concatenate all the databases 
all_csv_files = glob.glob('*.{}'.format('csv'))
df_concat = pd.concat([pd.read_csv(file) for file in all_csv_files], ignore_index=True)


# Get vacuum and water entries
vacuum_df = df_concat[df_concat['Env:']== 'vacuum'] 
water_df = df_concat[df_concat['Env:'] == 'water']

# From each get the entries for canonical and non-canonical bases
vacuum_df_can = vacuum_df[vacuum_df['RU_name'].isin(RUs_canonical)]
vacuum_df_noncan = vacuum_df[vacuum_df['RU_name'].isin(RUs_non_canonical)]
water_df_can = water_df[water_df['RU_name'].isin(RUs_canonical)]
water_df_noncan = water_df[water_df['RU_name'].isin(RUs_non_canonical)]


# Add field that says canonical or non
vacuum_df_can['Type_RU'] = ['canonical']*len(vacuum_df_can['RU_name'])
vacuum_df_noncan ['Type_RU'] = ['non-canonical']*len(vacuum_df_noncan['RU_name'])
water_df_can ['Type_RU'] = ['canonical']*len(water_df_can['RU_name'])
water_df_noncan ['Type_RU'] = ['non-canonical']*len(water_df_noncan['RU_name'])

vacuum_final = [vacuum_df_can, vacuum_df_noncan]
vacuum = pd.concat(vacuum_final)
water_final = [water_df_can, water_df_noncan]
water = pd.concat(water_final)


fig, axn = plt.subplots(2, 3, sharex=False, sharey=False, figsize=(24,16))
for i, (ax, count) in enumerate([(axn.flat[0], 1), 
                                            (axn.flat[1], 2), 
                                            (axn.flat[2], 3),
                                            (axn.flat[3], 4),
                                            (axn.flat[4], 5),
                                            (axn.flat[5], 6)]):
    
    # Plot the histogram with the distribution
    sns.set_context("paper", font_scale=2.5)
    y = sns.histplot(data=vacuum, x='ΔG°$_'+str(count)+'$', hue='Type_RU', multiple='dodge', ax=ax, kde = True,
                    palette=['green', 'gray'], legend=False, line_kws={'lw': 3, 'ls': '--'})
    y.bar_label(y.containers[0])
    y.bar_label(y.containers[1])
    ax.legend(['α', 'β'], fontsize=26, bbox_to_anchor=(1.5, 0.00), loc="lower right")
plt.savefig('ahora.pdf')
