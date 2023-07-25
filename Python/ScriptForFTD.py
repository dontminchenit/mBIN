import os
import scipy.io
import sys
import time

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NewFTD/")
from LoadDataFTD import loaddataFTD
from PathCovFTD import pathCovFTD
from ThicknessCovFTD import thicknessCovFTD
from ThicknessAtPath import thicknessAtPath
from distributionAnalysis import thicknessDist

# Change these directories to match the folder in your computer's path location of the github respository
# Only used to load the FTDGeneralData_20221114.mat file --> Saved as NetworkDataGeneral
baseDir='/Users/hyung/Research23_Network_Analysis/mBIN/Matlab'

# Location of the data folder
dataDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Data'

# Directory path where results will get written to. 
outputDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Output_Python_NewFTD'

# loads the preconstructed atlas data
NetworkDataGeneral = scipy.io.loadmat(os.path.join(baseDir, 'NetworkAnalysisGeneral', 'FTDGeneralData_20221114.mat'))


t0 = time.time()

# print( NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0])
# print("----")
# print( NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['R0_L1_Index'][0, 0])
# print("----")
# print( NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['R0_L1_Index'][0, 0].shape)
# print("----")
# print( NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['cdata'][0, 0])
# print("----")
# print( NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['cdata'][0, 0].shape)
# print("----")

# Loading Data into desired format
print('LOADING START')
pathCoM, pathT_GM, pathT_WM, pathLUT, HCResults, PatientTAUResults, PatientTDPResults, HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath = loaddataFTD(baseDir, dataDir, outputDir, NetworkDataGeneral)
print('LOADING DONE')

# # Thickness Data Distribution
# print('thicknessDist START')
# thicknessDist(outputDir, HCResults, PatientTAUResults, PatientTDPResults, NetworkDataGeneral)
# print('thicknessDist END')

# Pathology Analysis
print('PATH START')
LabelNames_Path, TAU_missing_index_GM, TDP_missing_index_GM, TAU_missing_index_WM, TDP_missing_index_WM = pathCovFTD(outputDir, NetworkDataGeneral, pathCoM, pathT_GM, pathT_WM, plotON = True)
print('PATH DONE') 

# # Thickness Analysis
# print('THICKNESS START')
# thicknessCovFTD(NetworkDataGeneral, pathLUT, HCResults, PatientTAUResults, PatientTDPResults, outputDir, plotON = True)
# print('THICKNESS DONE')

# # # LabelNames_Path = ['ANG_L', 'ATC_L', 'HIP_L', 'IFC_L', 'M1_L', 'MFC_L', 'OFC_L', 'PC_L', 'S1_L', 'SMTC_L', 'SPC_L', 'V1_L', 'aCING_L', 'aINS_L', 'aITC_L', 'dlPFC_L', 'iPFC_L', 'mPFC_L', 'pCING_L', 'pSTC_L', 'ANG_R', 'ATC_R', 'HIP_R', 'IFC_R', 'M1_R', 'MFC_R', 'OFC_R', 'PC_R', 'S1_R', 'SMTC_R', 'SPC_R', 'V1_R', 'aCING_R', 'aINS_R', 'aITC_R', 'dlPFC_R', 'iPFC_R', 'mPFC_R', 'pCING_R', 'pSTC_R']
# # Thickness Values Matching Pahtology Regions Analysis
# print('THICKNESS AT PATH START')
# thicknessAtPath(outputDir, NetworkDataGeneral, pathCoM, LabelNames_Path, HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath, HCResults, PatientTAUResults, PatientTDPResults, TAU_missing_index_GM, TDP_missing_index_GM, plotON = True)
# print('THICKNESS AT PATH DONE')

t1 = time.time()

total_n = t1-t0
print("time: " + str(total_n))