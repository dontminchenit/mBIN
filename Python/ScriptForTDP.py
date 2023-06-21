import os
import scipy.io
import sys
import time

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NewFTD/")
from LoadDataFTD import loaddataFTD
from PathCovFTD import pathCovFTD
from ThicknessCovFTD import thicknessCovFTD
from ThicknessAtPath import thicknessAtPath

# change these directories to match the folder in your computer's path location of the github respository / Only used to load the FTDGeneralData_20221114.mat file
baseDir='/Users/hyung/Research23_Network_Analysis/mBIN/Matlab'

# Location of the data folder
dataDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Data'

# Directory path where results will get written to. 
outputDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Output_Python_NewFTD'

# loads the preconstructed atlas data
NetworkDataGeneral = scipy.io.loadmat(os.path.join(baseDir, 'NetworkAnalysisGeneral', 'FTDGeneralData_20221114.mat'))


# Manual downsampling of GII vertices and faces
# NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['giiSurface_Both'][0,0]['vertices'][0, 0] = NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['giiSurface_Both'][0,0]['vertices'][0, 0][::2]
# NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['giiSurface_Both'][0,0]['faces'][0, 0] = NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['giiSurface_Both'][0,0]['faces'][0, 0][::2]


t0 = time.time()

# Start by trying to run these script (they're located in the "LBD" folder)
print('LOADING START')
pathCoM, pathT_GM, pathT_WM, pathNames_3D_Map, sn, pathLUT, HCResults, PatientTAUResults, PatientTDPResults, TAUthicknessAtPath, TDPthicknessAtPath, HCthicknessAtPath, TAUvolumeAtPath, TDPvolumeAtPath, HCvolumeAtPath  = loaddataFTD(baseDir, dataDir, outputDir, NetworkDataGeneral)
print('LOADING DONE')


print('PATH START')
LabelNames_Path, TAU_missing_index_GM, TDP_missing_index_GM, TAU_missing_index_WM, TDP_missing_index_WM = pathCovFTD(outputDir, NetworkDataGeneral, pathCoM, pathT_GM, pathT_WM, pathNames_3D_Map, sn, plotON = True)
print('PATH DONE')

# print('THICKNESS START')
# thicknessCovFTD(NetworkDataGeneral, pathLUT, HCResults, PatientTAUResults, PatientTDPResults, outputDir, plotON = True)
# print('THICKNESS DONE')

# LabelNames_Path = ['ANG_L', 'ATC_L', 'HIP_L', 'IFC_L', 'M1_L', 'MFC_L', 'OFC_L', 'PC_L', 'S1_L', 'SMTC_L', 'SPC_L', 'V1_L', 'aCING_L', 'aINS_L', 'aITC_L', 'dlPFC_L', 'iPFC_L', 'mPFC_L', 'pCING_L', 'pSTC_L', 'ANG_R', 'ATC_R', 'HIP_R', 'IFC_R', 'M1_R', 'MFC_R', 'OFC_R', 'PC_R', 'S1_R', 'SMTC_R', 'SPC_R', 'V1_R', 'aCING_R', 'aINS_R', 'aITC_R', 'dlPFC_R', 'iPFC_R', 'mPFC_R', 'pCING_R', 'pSTC_R']
print('THICKNESS AT PATH START')
thicknessAtPath(outputDir, NetworkDataGeneral, pathCoM, LabelNames_Path, TAUthicknessAtPath, TDPthicknessAtPath, HCthicknessAtPath, TAUvolumeAtPath, TDPvolumeAtPath, HCvolumeAtPath, HCResults, PatientTAUResults, PatientTDPResults, TAU_missing_index_GM, TDP_missing_index_GM, sn, plotON = True) #sn = 20
print('THICKNESS AT PATH DONE')

t1 = time.time()

total_n = t1-t0
print("time: " + str(total_n))