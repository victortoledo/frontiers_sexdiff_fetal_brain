#!/opt/python-3.6.1/bin/python3
#
# importing modules/packages
# 
print ("Importing modules/packages:")
import os
os.getcwd()
from netZooPy.panda.panda import Panda
import pandas as pd
import matplotlib.pyplot as plt
print ("Modules/packages loaded!")

#
# defining paths
# 
male_data='/home/users/victor/projects/sex_differences_brain_networks/panda_symbol/_h/male_exp.txt'
female_data='/home/users/victor/projects/sex_differences_brain_networks/panda_symbol/_h/female_exp.txt'
motif_data='/home/users/victor/projects/sex_differences_brain_networks/panda_symbol/_h/motif.txt'
ppi_data='/home/users/victor/projects/sex_differences_brain_networks/panda_symbol/_h/ppi.txt'
panda_male_output='/home/users/victor/projects/sex_differences_brain_networks/panda_symbol/_m/output_panda_male.txt'
panda_female_output='/home/users/victor/projects/sex_differences_brain_networks/panda_symbol/_m/output_panda_female.txt'
print ("Paths defined!")

print ("Load the motif data:")
motif_data=pd.read_csv(motif_data,sep="\t",header=None)
motif_data[0].unique().size
motif_data[1].unique().size
print ("Load the PPI data:")
ppi_data=pd.read_csv(ppi_data,sep="\t",header=None)
ppi_data.shape
#print ("Read the male expression data:")
#male_data=pd.read_csv(male_data,sep="\t",header=None)
#print ("Read the female expression data:") 
#female_data=pd.read_csv(female_data,sep="\t",header=None)

#
# calling panda
#
print ("Running PANDA for male samples:")
panda_obj = Panda(male_data, motif_data, ppi_data, save_tmp=True,save_memory = False,  remove_missing=False, keep_expression_matrix = False)
print ("Male PANDA network done! Saving!")
panda_obj.save_panda_results(panda_male_output)
print ("Running PANDA for female samples:")
panda_obj = Panda(female_data, motif_data, ppi_data, save_tmp=True,save_memory = False, remove_missing=False, keep_expression_matrix = False)
print ("Female PANDA network done! Saving!")
panda_obj.save_panda_results(panda_female_output)

print ("Code done!")
