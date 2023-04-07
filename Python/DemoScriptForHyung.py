import os
import scipy.io
import sys

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/LBD/")
import LoadDataLBD
import PathCovLBD_noAmy

# change these directories to match the folder in your computer's path location of the github respository:
baseDir='/Users/hyung/Research23_Network_Analysis/mBIN/Matlab'

# Location of the data folder
dataDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Data'

# Directory path where results will get written to. 
outputDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Output_Python'

# loads the preconstructed atlas data
NetworkDataGeneral = scipy.io.loadmat(os.path.join(baseDir, 'NetworkAnalysisGeneral', 'FTDGeneralData_20221114.mat'))

# Start by trying to run these script (they're located in the "LBD" folder)
# pathCoM, pathNamesRaw, pathDataGM, LBD_yesADIndx, LBD_noADIndx, sn = LoadDataLBD.loaddataLBD(baseDir, dataDir, outputDir, NetworkDataGeneral)

# PathCovLBD_noAmy.pathCovLBD_noAmy(outputDir, NetworkDataGeneral, pathCoM, pathNamesRaw, pathDataGM, LBD_yesADIndx, LBD_noADIndx, sn)




temp = NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]
print(temp)

# print(type(NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['cdata'][0, 0]))
# print((NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['cdata'][0, 0]).shape)