# Import libraries
import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import matplotlib
import math
from scipy.stats import multivariate_normal
from scipy.stats import circvar
from scipy.stats import circmean
from scipy.stats import circstd
import scipy.stats as st
import matplotlib.patches as mpatches
import shutil
import glob

environment = ["vacuum","water"]
TC_names = ["2deoxy_ribofuranose", "ribofuranose", "2deoxy_ribopyranose", "ribopyranose", "threose"]
alpha = ["2endo_alpha", "3endo_alpha", "D1C4_alpha", "D4C1_alpha"]
beta = ["2endo_beta", "3endo_beta", "D1C4_beta", "D4C1_beta"]

#####################################
#                                   #
# OBTAIN FREQUENCY TABLE FOR ANGLES #
#                                   #
#####################################

def frequencies_for_hist(angles_dataframe, lower_value, upper_value, step):
    # Create the bins or intervals
    bins = pd.interval_range(lower_value, upper_value, freq=step)    
    # Assign each value of angle to bin and count bin occurences, angles is the corresponding column of dataframe
    counts = pd.cut(angles_dataframe, bins).value_counts(sort=False)
    # Create frequency table
    freq_table = pd.Series(counts.values, index=bins)
    freq_table = pd.DataFrame({'Bins': freq_table.index, 'Count': freq_table.values})
    frequencies = np.array(freq_table['Count'])
    return (freq_table, frequencies)
    





# Concatenate all the dataframes
cwd = os.getcwd()
for i, folder in enumerate(TC_names):
    source = str(cwd)+"/"+str(folder)+"/summary_db_RU_pos.csv"
    dest = str(cwd)+"/summary_db_RU_pos"+str(i)+".csv"
    shutil.copy(source, dest)

all_csv_files = glob.glob('*.{}'.format('csv'))
df = pd.concat([pd.read_csv(file) for file in all_csv_files], ignore_index=True)
df.to_csv('summary_db_RU_pos.csv')

# Copy file that has IL
source = str(cwd)+"/ribofuranose/IL_list_names.txt"
dest = str(cwd)+"/IL_list_names.txt"
shutil.copy(source, dest)
IL = np.loadtxt("IL_list_names.txt", dtype='str')

# Load the csv database
db_all = pd.read_csv('summary_db_RU_pos.csv')
    
# Obtain entries for vacuum and water environment from db
db_vacuum = db_all[db_all['Environment'] == "vacuum"]
db_water = db_all[db_all['Environment'] == "water"]
    
# Obtain from db_env the entries for canonical and for not_canonical bases
# List of RUs
canonical_RU = ['adenine', 'guanine', 'cytosine', 'thymine', 'uracil']
non_canonical_RU = ['TARC_cbond', 'TARC_nbond', 'BA_cbond', 'BA_nbond', 'melamine', 'CA']

#vacuum
db_vacuum_canonical = db_vacuum[db_vacuum['RU'].isin(canonical_RU)]
chi_vacuum_canonical = db_vacuum_canonical['Chi']
db_vacuum_non_canonical = db_vacuum[db_vacuum['RU'].isin(non_canonical_RU)]
chi_vacuum_non_canonical = db_vacuum_non_canonical['Chi']

# Water
db_water_canonical = db_water[db_water['RU'].isin(canonical_RU)]
chi_water_canonical = db_water_canonical['Chi']
db_water_non_canonical = db_water[db_water['RU'].isin(non_canonical_RU)]
chi_water_non_canonical = db_water_non_canonical['Chi']
       
# Counts for alpha and beta
# Vacuum
all_vacuum_alpha = db_vacuum[db_vacuum['Ring_conf'].isin(alpha)]
all_vacuum_beta = db_vacuum[db_vacuum['Ring_conf'].isin(beta)]   
can_vacuum_alpha = db_vacuum_canonical[db_vacuum_canonical['Ring_conf'].isin(alpha)]
can_vacuum_beta = db_vacuum_canonical[db_vacuum_canonical['Ring_conf'].isin(beta)]
non_can_vacuum_alpha = db_vacuum_non_canonical[db_vacuum_non_canonical['Ring_conf'].isin(alpha)]
non_can_vacuum_beta = db_vacuum_non_canonical[db_vacuum_non_canonical['Ring_conf'].isin(beta)]
    
# Water
all_water_alpha = db_water[db_water['Ring_conf'].isin(alpha)]
all_water_beta = db_water[db_water['Ring_conf'].isin(beta)]   
can_water_alpha = db_water_canonical[db_water_canonical['Ring_conf'].isin(alpha)]
can_water_beta = db_water_canonical[db_water_canonical['Ring_conf'].isin(beta)]
non_can_water_alpha = db_water_non_canonical[db_water_non_canonical['Ring_conf'].isin(alpha)]
non_can_water_beta = db_water_non_canonical[db_water_non_canonical['Ring_conf'].isin(beta)]
     
################### Obtain values of phi2 angle #########################
chi_all_vacuum_alpha = all_vacuum_alpha['Chi']
chi_all_vacuum_beta = all_vacuum_beta['Chi']
chi_can_vacuum_alpha = can_vacuum_alpha ['Chi']
chi_can_vacuum_beta = can_vacuum_beta ['Chi']
chi_non_can_vacuum_alpha = non_can_vacuum_alpha ['Chi']
chi_non_can_vacuum_beta = non_can_vacuum_beta ['Chi']

# Water
chi_all_water_alpha = all_water_alpha['Chi']
chi_all_water_beta = all_water_beta['Chi']
chi_can_water_alpha = can_water_alpha ['Chi']
chi_can_water_beta = can_water_beta ['Chi']
chi_non_can_water_alpha = non_can_water_alpha ['Chi']
chi_non_can_water_beta = non_can_water_beta ['Chi']
       
############### Obtain the frequency tables and counts ##################
# Vacuum
(freq_df_all_vacuum_alpha, count_chi_all_vacuum_alpha) = frequencies_for_hist(chi_all_vacuum_alpha, 0, 360, 90)
(freq_df_all_vacuum_beta, count_chi_all_vacuum_beta) = frequencies_for_hist(chi_all_vacuum_beta, 0, 360, 90)  
                                                                                                                                                                                                                                                                                                                                           
(freq_df_can_vacuum_alpha, count_chi_can_vacuum_alpha) = frequencies_for_hist(chi_can_vacuum_alpha, 0, 360, 90)
(freq_df_can_vacuum_beta, count_chi_can_vacuum_beta) = frequencies_for_hist(chi_can_vacuum_beta, 0, 360, 90)                                                                             
(freq_df_non_can_vacuum_alpha, count_chi_non_canonical_vacuum_alpha) = frequencies_for_hist(chi_non_can_vacuum_alpha, 0, 360, 90)
(freq_df_non_can_vacuum_beta, count_chi_non_canonical_vacuum_beta) = frequencies_for_hist(chi_non_can_vacuum_beta, 0, 360, 90) 
                                                                                                                                                                                             
# Water
(freq_df_all_water_alpha, count_chi_all_water_alpha) = frequencies_for_hist(chi_all_water_alpha, 0, 360, 90)
(freq_df_all_water_beta, count_chi_all_water_beta) = frequencies_for_hist(chi_all_water_beta, 0, 360, 90)  
                                                                                                                                                                                                                                                                                                                                           
(freq_df_can_water_alpha, count_chi_can_water_alpha) = frequencies_for_hist(chi_can_water_alpha, 0, 360, 90)
(freq_df_can_water_beta, count_chi_can_water_beta) = frequencies_for_hist(chi_can_water_beta, 0, 360, 90)                                                                             
(freq_df_non_can_water_alpha, count_chi_non_canonical_water_alpha) = frequencies_for_hist(chi_non_can_water_alpha, 0, 360, 90)
(freq_df_non_can_water_beta, count_chi_non_canonical_water_beta) = frequencies_for_hist(chi_non_can_water_beta, 0, 360, 90) 
                                                                                            
                                                                                        
# Obtain the circ mean, std and total datapoints 
# Vacuum
chi_can_vacuum_rad = np.radians(chi_vacuum_canonical)
(circ_mean_vacuum_canonical, circ_std_vacuum_canonical) = (np.degrees(circmean(chi_can_vacuum_rad)),
                                                           np.degrees(circstd(chi_can_vacuum_rad)))
                                                                      
chi_non_can_vacuum_rad = np.radians(chi_vacuum_non_canonical)    
(circ_mean_vacuum_non_canonical, circ_std_vacuum_non_canonical) = (np.degrees(circmean(chi_non_can_vacuum_rad)),
                                                                   np.degrees(circstd(chi_non_can_vacuum_rad)))

chi_all_vacuum_rad = np.radians(db_vacuum['Chi'])    
(circ_mean_all_vacuum, circ_std_all_vacuum) = (np.degrees(circmean(chi_all_vacuum_rad)),
                                                                   np.degrees(circstd(chi_all_vacuum_rad)))

# Water
chi_can_water_rad = np.radians(chi_water_canonical)
(circ_mean_water_canonical, circ_std_water_canonical) = (np.degrees(circmean(chi_can_water_rad)),
                                                         np.degrees(circstd(chi_can_water_rad)))
                                                                      
chi_non_can_water_rad = np.radians(chi_water_non_canonical)    
(circ_mean_water_non_canonical, circ_std_water_non_canonical) = (np.degrees(circmean(chi_non_can_water_rad)),
                                                                 np.degrees(circstd(chi_non_can_water_rad)))

chi_all_water_rad = np.radians(db_water['Chi'])    
(circ_mean_all_water, circ_std_all_water) = (np.degrees(circmean(chi_all_water_rad)),
                                                                 np.degrees(circstd(chi_all_water_rad)))
                                                                     
## Define vectors with all info for rose diagrams
#Vacuum
info_can_vacuum = [count_chi_can_vacuum_alpha, count_chi_can_vacuum_beta, circ_mean_vacuum_canonical, circ_std_vacuum_canonical]
info_non_can_vacuum = [count_chi_non_canonical_vacuum_alpha, count_chi_non_canonical_vacuum_beta, circ_mean_vacuum_non_canonical, circ_std_vacuum_non_canonical]
info_all_vacuum = [count_chi_all_vacuum_alpha, count_chi_all_vacuum_beta, circ_mean_all_vacuum, circ_std_all_vacuum]

# Water
info_can_water = [count_chi_can_water_alpha, count_chi_can_water_beta, circ_mean_water_canonical, circ_std_water_canonical]
info_non_can_water = [count_chi_non_canonical_water_alpha, count_chi_non_canonical_water_beta, circ_mean_water_non_canonical, circ_std_water_non_canonical]
info_all_water = [count_chi_all_water_alpha, count_chi_all_water_beta, circ_mean_all_water, circ_std_all_water]       
    
    
    
########## Create the rose diagrams #############
fig, axn = plt.subplots(2, 3, sharey=False, sharex=False, subplot_kw=dict(projection="polar"), figsize=(20,14))
for i, (ax, count_alpha, count_beta, mean, std) in enumerate([(axn.flat[0], info_can_vacuum[0], info_can_vacuum[1], 
                                                               info_can_vacuum[2], info_can_vacuum[3]), 
                                                              (axn.flat[1], info_non_can_vacuum[0], 
                                                               info_non_can_vacuum[1], info_non_can_vacuum[2],
                                                               info_non_can_vacuum[3]), 
                                                              (axn.flat[2], info_all_vacuum[0], 
                                                               info_all_vacuum[1], info_all_vacuum[2],
                                                               info_all_vacuum[3]), 
                                                              (axn.flat[3], info_can_water[0], info_can_water[1], 
                                                               info_can_water[2], info_can_water[3]), 
                                                              (axn.flat[4], info_non_can_water[0], 
                                                               info_non_can_water[1], info_non_can_water[2],
                                                               info_non_can_water[3]), 
                                                              (axn.flat[5], info_all_water[0], info_all_water[1], 
                                                               info_all_water[2], info_all_water[3])]):
    # Obtain the quadrants
    theta = np.radians(np.arange(0, 360, 90)) 
    mean = round(mean, 1)
    std = round(std, 1)
    
    # Obtain the range for freqs
    width = np.radians(90)  
    
    # Create bar plot
    bars = ax.bar(theta, count_alpha, width = width, edgecolor = 'black', align = 'edge', facecolor='gray', 
                  linewidth = 2.0, label = 'Alpha')
    
    bars2 = ax.bar(theta, count_beta, width = width, edgecolor = 'black', align = 'edge', facecolor='g', 
                   bottom = count_alpha, linewidth = 2.0, label = 'Beta')
    
    # Configure the grids
    #ax.yaxis.grid(True,color='k',linestyle='--', linewidth=0.5)
    #ax.xaxis.grid(True,color='black',linestyle='-', linewidth=0.8)    
    ax.xaxis.set_visible(True)
    ax.set_facecolor('xkcd:white')
    # Modify grid 
    ax.grid(color='grey', linestyle='--')
    ax.grid(True)
    
    # Configure the axis ticks
    ax.tick_params(axis='x', colors='black', labelsize='22')
    ax.tick_params(axis='y', colors='blue', labelsize='22')
    #plt.xticks(weight = 'bold')
    labelss = ax.get_yticklabels()
    for label in labelss:
        label.set_fontweight('bold')
    
    # Change position of y labels
    ax.set_rlabel_position(90)
    
    labels=['0°', "90°", '180°', "270°"]
    ax.set_xticks(theta)
    #ax.set_xticklabels(labels, fontweight='bold') 
    
    x_tail = 0.0
    y_tail = 0.0
    x_head = float(mean)/180.*np.pi
    array = count_alpha + count_beta
    y_head = max(array)
    dx = x_head - x_tail
    dy = y_head - y_tail
    arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (dx, dy), mutation_scale=30, facecolor = 'red', edgecolor = 'black')
    ax.add_patch(arrow)
    
    # Add the title to rose diagram
    if i == 0:
        ax.set_title('Canonical RUs\n\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" +  "±" + str(std) + "°"
                     + ")" + '\n', fontweight="bold", fontsize=28)
    elif i == 1:
        ax.set_title('Non-canonical RUs\n\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°"
                     + ")" + '\n', fontweight="bold", fontsize=28)
    elif i == 2:
        ax.set_title('All RUs\n\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" + 
                     '\n', fontweight="bold", fontsize=28)
    elif i == 3:
        ax.set_title('\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" + 
                     '\n', fontweight="bold", fontsize=28)
    elif i == 4:
        ax.set_title('\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" 
                     + '\n', fontweight="bold", fontsize=28)
    elif i == 5:
        ax.set_title('\n(' + r'$\bar\chi$' + '=' + str(mean) + "°" + "±" + str(std) + "°" + ")" +
                     '\n', fontweight="bold", fontsize=28)
        ax.legend(['α', 'β'], fontsize=24, bbox_to_anchor=(1.3, 0.00), loc="lower right")
        leg = ax.get_legend()
        leg.legendHandles[0].set_color('grey')
        leg.legendHandles[1].set_color('green') 
     
    # little space between the bar and the label
    labelPadding = 4
    angles = np.radians(np.arange(45, 405, 90)) 
    
    # Set the bar numbers (First way)
    # Set the bar numbers (alpha anomers)
    for bar, angle, height, label in zip(bars, angles, count_alpha, count_alpha):
        # Labels are rotated. Rotation must be specified in degrees :(
        rotation = np.rad2deg(angle)        
        # Flip some labels upside down
        alignment = ""
        if angle >= np.pi/2 and angle < 3*np.pi/2:
            alignment = "right"
            rotation = rotation + 180
        else: 
            alignment = "left"
        if height != 0 and label != 0:
            ax.text(x=angle, 
                y=bar.get_height(), 
                s=label, 
                ha=alignment, 
                va='center', 
                rotation=rotation, 
                rotation_mode="anchor",
                fontsize=18,
                fontweight='bold')
                    
    # Set the bar numbers (First way)
    for bar, angle, height, label, ad_height in zip(bars2, angles, count_beta, count_beta, count_alpha):        
        # Labels are rotated. Rotation must be specified in degrees :(
        rotation = np.rad2deg(angle)        
        # Flip some labels upside down
        alignment = ""
        if angle >= np.pi/2 and angle < 3*np.pi/2:
            alignment = "right"
            rotation = rotation + 180
        else: 
            alignment = "left"
        if height != 0 and label != 0:
            ax.text(
                x=angle, 
                y=bar.get_height()+ad_height, 
                s=label, 
                ha=alignment, 
                va='center', 
                rotation=rotation, 
                rotation_mode="anchor",
                fontsize=18,
                fontweight='bold')
       
plt.tight_layout()       
plt.savefig('Rose_diagram_'+str(IL[1])+'.pdf')                                                  
      
