import pandas as pd
import os
import numpy as np

dataDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Data'

# pathT: ex-vivo histopathology Data (Quantification) / %AO
new_pathT = pd.read_excel(os.path.join(dataDir, 'NewFTDData', 'FTLD_Library_3-31-22_update.xlsx'), dtype={'INDDID': str})

# INDDID list
new_path_INDDID_list = new_pathT['INDDID'].tolist()
# Make it unique
new_path_INDDID_list = np.unique(np.array(new_path_INDDID_list))
new_path_INDDID_list = list(map(str, new_path_INDDID_list)) # cast to str

# MRI Thickness and Volume values for Patient Cortical Atlas / 'id' includes the case with 108783x09 -> so set the dtype as str
# Same format as Original --> No need to change
thicknessAllraw = pd.read_csv(os.path.join(dataDir, 'NewFTDData', 'invivoPathCohort_quantsSubSesSchaefer400_tian12.csv'), dtype={'id': str})

# INDDID list
thicknessAllraw_INDDID_list = thicknessAllraw['id'].tolist()
# Make it unique
thicknessAllraw_INDDID_list = np.unique(np.array(thicknessAllraw_INDDID_list))
thicknessAllraw_INDDID_list = list(map(str, thicknessAllraw_INDDID_list)) # cast to str

# Spreadsheet with a column that you can use to assign tau/tdp/hc for the thickness data
path_typeT = pd.read_excel(os.path.join(dataDir, 'NewFTDData', 'InvivoPathCohort_03172023.xls'), dtype={'INDDID': str})

# INDDID list
path_type_INDDID_list = path_typeT['INDDID'].tolist()
# Make it unique
path_type_INDDID_list = np.unique(np.array(path_type_INDDID_list))
path_type_INDDID_list = list(map(str, path_type_INDDID_list)) # cast to str


path_type = path_typeT.groupby('Group')

type_HC = path_type.get_group('HC')
type_tau = path_type.get_group('tau')
type_tdp = path_type.get_group('tdp')

type_HC_INDDID_list = type_HC['INDDID'].tolist()
type_tau_INDDID_list = type_tau['INDDID'].tolist()
type_tdp_INDDID_list = type_tdp['INDDID'].tolist()

# print(type_HC_INDDID_list)
# print(type_tau_INDDID_list)
# print(type_tdp_INDDID_list)
# print(len(type_HC_INDDID_list))
# print(len(type_tau_INDDID_list))
# print(len(type_tdp_INDDID_list))
# print(len(path_typeT.index))

# print(len(new_path_INDDID_list)) # Pathology
# print(len(thicknessAllraw_INDDID_list)) # MR Thickness
# print(len(path_type_INDDID_list)) # MR Thickness - Type

#print((new_path_INDDID_list)) # Pathology
#print((thicknessAllraw_INDDID_list)) # MR Thickness
#print((path_type_INDDID_list)) # MR Thickness - Type

# Check if MR Thickness - Type INDDID vs Pathology INDDID
test1 = np.sum(np.in1d(path_type_INDDID_list, new_path_INDDID_list)) # Number of INDDID of Type spreadsheet in Pathology INDDID / out of 114 unique INDDID
test2 = np.sum(np.in1d(new_path_INDDID_list, path_type_INDDID_list)) # Number of INDDID of Pathology Dataset in Type spreadsheet INDDID / out of 180 unique INDDID

# Check if MR Thickness - Type INDDID vs MR INDDID
test3 = np.sum(np.in1d(path_type_INDDID_list, thicknessAllraw_INDDID_list)) # Number of INDDID of Type spreadsheet in MR Thickness INDDID / out of 114 unique INDDID
test4 = np.sum(np.in1d(thicknessAllraw_INDDID_list, path_type_INDDID_list)) # Number of INDDID of MR Thickness in Type spreadsheet INDDID / out of 111 unique INDDID

print('Number of INDDID of Type spreadsheet in Pathology INDDID: ' + str(test1) + ' out of 114 unique INNDDID')
print('The missing INDDID is: ' + str(np.take(path_type_INDDID_list, np.argwhere(np.in1d(path_type_INDDID_list, new_path_INDDID_list) == False).flatten())))
print('\n')

print('Number of INDDID of Pathology Dataset in Type spreadsheet INDDID: ' + str(test2) + ' out of 180 unique INNDDID')
print('The missing INDDID is: ' + str(np.take(new_path_INDDID_list, np.argwhere(np.in1d(new_path_INDDID_list, path_type_INDDID_list) == False).flatten())))
print('\n')

print('Number of INDDID of Type spreadsheet in MR Thickness INDDID: ' + str(test3) + ' out of 114 unique INNDDID')
print('The missing INDDID is: ' + str(np.take(path_type_INDDID_list, np.argwhere(np.in1d(path_type_INDDID_list, thicknessAllraw_INDDID_list) == False).flatten())))
print('\n')

print('Number of INDDID of MR Thickness in Type spreadsheet INDDID: ' + str(test4) + ' out of 111 unique INNDDID')
print('The missing INDDID is: ' + str(np.take(thicknessAllraw_INDDID_list, np.argwhere(np.in1d(thicknessAllraw_INDDID_list, path_type_INDDID_list) == False).flatten())))




