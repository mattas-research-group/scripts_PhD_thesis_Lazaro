import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import seaborn as sb
import scipy as scipy
from scipy import stats
from matplotlib.gridspec import SubplotSpec

IL_list = ['phosphate', 'arsenate']
TC_list = ['2deoxy_ribofuranose', 'ribofuranose', 'threose', '2deoxy_ribopyranose', 'ribopyranose']


################# Merge dE, dE(ZPE) and dG values in one cell for tables in paper  
for IL in IL_list:
    
    os.chdir(IL) # Open IL folder
    
    for TC in TC_list:
        
        os.chdir(TC) # Open TC folder
        
        # Read csv file    
        df = pd.read_csv('db_'+str(TC)+'_ring_puckering_params.csv')
        
        # Loop through the different TCs
        if TC == "2deoxy_ribofuranose" or TC == "ribofuranose" or TC == "threose": # Case of furanoses

            cols = ['phi2(Phase_angle)_Cremer&Pople', 'Q(Total_amplitude)_Cremer&Pople']

        else: # Case of pyranoses
        
            cols = ['phi2(phase_angle)_Cremer&Pople', 'theta(phase_angle)_Cremer&Pople', 'Q(Total_amplitude)_Cremer&Pople']
    
        # Round values
        v = df[cols].round(1)
        df[cols] = v
        df['combined_CP_params'] = df[cols].apply(lambda row: '\n'.join(row.values.astype(str)), axis=1) 
        df.to_csv('df_with_CPvalues_to_copy.csv')    
        
        os.chdir("..")
    os.chdir("..")
