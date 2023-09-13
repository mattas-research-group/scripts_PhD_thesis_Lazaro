import sys, os
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from numpy import loadtxt
import seaborn as sb
import scipy as scipy
from scipy import stats
from itertools import repeat
from matplotlib.gridspec import SubplotSpec
import shutil


def calc (df, list_conf):
    # Get unique values of TC
    TC = df['TC'].unique()
    list_compared_systems = []
    list_E = []
    list_corr = []
    list_G = []
    for conf1 in list_conf:
        for conf2 in list_conf:
            for conf3 in list_conf:
                for conf4 in list_conf:
                    first_t = df[df['TC'].isin([TC[0]])&df['RU'].isin(['thymine'])&df['Initial ring conf'].isin([conf1])]
                    second_t = df[df['TC'].isin([TC[1]])&df['RU'].isin(['uracil'])&df['Initial ring conf'].isin([conf2])]
                    third_t = df[df['TC'].isin([TC[0]])&df['RU'].isin(['uracil'])&df['Initial ring conf'].isin([conf3])]
                    fourth_t = df[df['TC'].isin([TC[1]])&df['RU'].isin(['thymine'])&df['Initial ring conf'].isin([conf4])]
            
                    # Obtain energy
                    dE = first_t['zpE_nucleoside(kJ/mol)'].to_numpy() + second_t['zpE_nucleoside(kJ/mol)'].to_numpy() - third_t['zpE_nucleoside(kJ/mol)'].to_numpy() - fourth_t['zpE_nucleoside(kJ/mol)'].to_numpy()
                    dcorrE = first_t['corrE_nucleoside(kJ/mol)'].to_numpy() + second_t['corrE_nucleoside(kJ/mol)'].to_numpy() - third_t['corrE_nucleoside(kJ/mol)'].to_numpy() - fourth_t['corrE_nucleoside(kJ/mol)'].to_numpy()            
                    dG = first_t['G_nucleoside(kJ/mol)'].to_numpy() + second_t['G_nucleoside(kJ/mol)'].to_numpy() - third_t['G_nucleoside(kJ/mol)'].to_numpy() - fourth_t['G_nucleoside(kJ/mol)'].to_numpy()
                  
                    # Obtain compared_systems
                    compared_sys = [TC[0]+'_'+conf1+'_T' + ' + ' + TC[1]+'_'+conf2+'_U' + 
                                    ' - ' + TC[0]+'_'+conf3+'_U' + ' - ' +  TC[1]+'_'+conf4+'_T']
                    
                    
                    # Append to final empty vectors
                    list_compared_systems.append(TC[0]+'_'+conf1+'_T' + ' + ' + TC[1]+'_'+conf2+'_U' + 
                                    ' - ' + TC[0]+'_'+conf3+'_U' + ' - ' +  TC[1]+'_'+conf4+'_T')
                    list_E.append(float(dE))
                    list_corr.append(float(dcorrE))
                    list_G.append(float(dG))                   
    
    return [list_compared_systems, list_E, list_corr, list_G]
    



# Read general dataframe
initial_df = pd.read_csv('condensation_reaction_allData_nucleosides.csv')


# Obtain from the dataframe the entries that contain T and U
RU_list = ["thymine", "uracil"]
T_U_df = initial_df[initial_df['RU'].isin(RU_list)]

# Obtain from T_U_df two dataframes from vacuum and water
vacuum_T_U_df = T_U_df[T_U_df['Environment'].isin(['vacuum'])]
water_T_U_df = T_U_df[T_U_df['Environment'].isin(['water'])]

# Obtain from the different environment dataframes the entries belonging to 2dRibf, Ribf, 2dRib and Rib
fura_list = ["2deoxy_ribofuranose", "ribofuranose"]
fura_vacuum_T_U_df = vacuum_T_U_df[vacuum_T_U_df['TC'].isin(fura_list)]
fura_water_T_U_df = water_T_U_df[water_T_U_df['TC'].isin(fura_list)]

pyra_list = ["2deoxy_ribopyranose", "ribopyranose"]
pyra_vacuum_T_U_df = vacuum_T_U_df[vacuum_T_U_df['TC'].isin(pyra_list)]
pyra_water_T_U_df = water_T_U_df[water_T_U_df['TC'].isin(pyra_list)]

# Obtain entries for alpha and beta
alpha_fura_list = ["2endo_alpha", "3endo_alpha"]
fura_vacuum_T_U_df_alpha = fura_vacuum_T_U_df[fura_vacuum_T_U_df['Initial ring conf'].isin(alpha_fura_list)]
fura_water_T_U_df_alpha = fura_water_T_U_df[fura_water_T_U_df['Initial ring conf'].isin(alpha_fura_list)]
alpha_pyra_list = ["D1C4_alpha", "D4C1_alpha"]
pyra_vacuum_T_U_df_alpha = pyra_vacuum_T_U_df[pyra_vacuum_T_U_df['Initial ring conf'].isin(alpha_pyra_list)]
pyra_water_T_U_df_alpha = pyra_water_T_U_df[pyra_water_T_U_df['Initial ring conf'].isin(alpha_pyra_list)]


beta_fura_list = ["2endo_beta", "3endo_beta"]
fura_vacuum_T_U_df_beta = fura_vacuum_T_U_df[fura_vacuum_T_U_df['Initial ring conf'].isin(beta_fura_list)]
fura_water_T_U_df_beta = fura_water_T_U_df[fura_water_T_U_df['Initial ring conf'].isin(beta_fura_list)]
beta_pyra_list = ["D1C4_beta", "D4C1_beta"]
pyra_vacuum_T_U_df_beta = pyra_vacuum_T_U_df[pyra_vacuum_T_U_df['Initial ring conf'].isin(beta_pyra_list)]
pyra_water_T_U_df_beta = pyra_water_T_U_df[pyra_water_T_U_df['Initial ring conf'].isin(beta_pyra_list)]




list_thermo = ['zpE_nucleoside(kJ/mol)', 'corrE_nucleoside(kJ/mol)', 'G_nucleoside(kJ/mol)']
final_df = pd.DataFrame({"env": [], "compared_systems": [], "ΔE": [], "ΔE(ZPE)": [], "ΔG": []})

#list_df = [fura_vacuum_T_U_df_alpha, fura_water_T_U_df_alpha, pyra_vacuum_T_U_df_alpha, pyra_water_T_U_df_alpha,
#          fura_vacuum_T_U_df_beta, fura_water_T_U_df_beta, pyra_vacuum_T_U_df_beta, pyra_water_T_U_df_beta]
list_df = [fura_vacuum_T_U_df_alpha, fura_water_T_U_df_alpha, pyra_vacuum_T_U_df_alpha, pyra_water_T_U_df_alpha, fura_vacuum_T_U_df_beta, fura_water_T_U_df_beta, pyra_vacuum_T_U_df_beta, pyra_water_T_U_df_beta]

for i, df in enumerate(list_df):
    env = df['Environment'].iloc[0]
    #print(env)
    conf_list = df['Initial ring conf'].unique()
    # Obtain values and append to final dataframe
    [comp_sys, dE, dcorrE, dG] = calc (df[['Environment', 'TC', 'Initial ring conf', 'RU', list_thermo[0], list_thermo[1], list_thermo[2]]], conf_list)
    df2 = pd.DataFrame({"env": df['Environment'].iloc[0],
                        "compared_systems": comp_sys, 
                        "ΔE": dE,
                        "ΔE(ZPE)": dcorrE,
                        "ΔG": dG})
                        
    final_df = pd.concat([final_df, df2])

# Convert final dataframe to csv file
final_df.to_csv('T_U_comparison.csv', index=False)

