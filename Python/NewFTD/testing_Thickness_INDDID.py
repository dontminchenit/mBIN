import pandas as pd
import os
import numpy as np

dataDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Data'

# demoT: MRI Data - Demographics 
# measureT: clinical test, measurement (Traditional  Neurological tests in measuring disease severity)
demoT = pd.read_excel(os.path.join(dataDir, 'LBDData', 'MRI Schaeffer Demographics classification.xlsx'))
measuresT = pd.read_excel(os.path.join(dataDir, 'LBDData', 'LBD neuropsych annualized change no formulas.xlsx')) 

# --------------------- THICKNESS PART ---------------------
# MRI Thickness value for All Subjects - schaefer400x7
thicknessAllraw = pd.read_csv(os.path.join(dataDir, 'NewFTDData', 'invivoPathCohort_quantsSubSesSchaefer400_tian12.csv'), dtype={'id': str})

# Look Up Table for Type of MRI Thickness Subjects
thicknessPathLUT = pd.read_excel(os.path.join(dataDir, 'NewFTDData', 'InvivoPathCohort_03172023.xls'), dtype={'INDDID': str})

# Join the Above two dataframes on INDDID (Keep only the ones that INDDID are overlapping)
thicknessAll = pd.merge(thicknessAllraw, thicknessPathLUT, left_on='id', right_on='INDDID', how='inner') # We only lose INDDID 108783x09 in the thicknessAllraw (849 rows lost)

# Group by path type
thickness_path_type = thicknessAll.groupby('Group')

# MRI Thickness values for Healthy Control
thicknessHC = thickness_path_type.get_group('HC')

# MRI Thickness values for Patient (TAU)
thicknessPatientTau = thickness_path_type.get_group('tau')

# MRI Thickness values for Patient (TDP)
thicknessPatientTdp = thickness_path_type.get_group('tdp')

# Measured INDDID
measure_IDs = measuresT.INDDID

# IDs
HC_IDs = np.unique(thicknessHC.INDDID)
Tau_IDs = np.unique(thicknessPatientTau.INDDID)
Tdp_IDs = np.unique(thicknessPatientTdp.INDDID)

print(len(measure_IDs)) # 117
print(len(HC_IDs)) # 54
print(len(Tau_IDs)) # 26
print(len(Tdp_IDs)) # 30

print(np.in1d(HC_IDs, measure_IDs).shape)

test1 = np.sum(np.in1d(HC_IDs, measure_IDs)) # Number of Thickness HC INDDID in Measure INDDID 
test2 = np.sum(np.in1d(Tau_IDs, measure_IDs)) # Number of Thickness TAU INDDID in Measure INDDID 
test3 = np.sum(np.in1d(Tdp_IDs, measure_IDs)) # Number of Thickness TDP INDDID in Measure INDDID 

print('Number of Thickness HC INDDID in Measure INDDID : ' + str(test1) + ' out of 54 unique INNDDID')
#print('The missing INDDID is: ' + str(np.take(HC_IDs, np.argwhere(np.in1d(HC_IDs, measure_IDs) == False).flatten())))
print('\n')

print('Number of Thickness TAU INDDID in Measure INDDID : ' + str(test2) + ' out of 26 unique INNDDID')
#print('The missing INDDID is: ' + str(np.take(Tau_IDs, np.argwhere(np.in1d(Tau_IDs, measure_IDs) == False).flatten())))
print('\n')

print('Number of Thickness TDP INDDID in Measure INDDID : ' + str(test3) + ' out of 30 unique INNDDID')
#print('The missing INDDID is: ' + str(np.take(Tdp_IDs, np.argwhere(np.in1d(Tdp_IDs, measure_IDs) == False).flatten())))
print('\n')

