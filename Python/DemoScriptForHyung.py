import os
import scipy.io
import sys

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/LBD/")
import LoadDataLBD
import PathCovLBD_noAmy

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/NetworkAnalysisGeneral/Code/SupportFunctions/")
from plotNetwork3 import plotNetwork3

# change these directories to match the folder in your computer's path location of the github respository:
baseDir='/Users/hyung/Research23_Network_Analysis/mBIN/Matlab'

# Location of the data folder
dataDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Data'

# Directory path where results will get written to. 
outputDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Output_Python'

# loads the preconstructed atlas data
NetworkDataGeneral = scipy.io.loadmat(os.path.join(baseDir, 'NetworkAnalysisGeneral', 'FTDGeneralData_20221114.mat'))

# Start by trying to run these script (they're located in the "LBD" folder)
pathCoM, pathNamesRaw, pathDataGM, LBD_yesADIndx, LBD_noADIndx, sn = LoadDataLBD.loaddataLBD(baseDir, dataDir, outputDir, NetworkDataGeneral)

PathCovLBD_noAmy.pathCovLBD_noAmy(outputDir, NetworkDataGeneral, pathCoM, pathNamesRaw, pathDataGM, LBD_yesADIndx, LBD_noADIndx, sn)



# ---------------------------Testing Plot3D---------------------------
# import numpy as np
# import matplotlib.pyplot as plt
# covMat_yesAD = np.array([[np.nan, 0.40864279, 0.82826751, 0.85717921, 0.87531807],
#                          [0.40864279, np.nan, 0.40310775, 0.42895138, 0.53611699],
#                          [0.82826751, 0.40310775, np.nan, 0.54380405, 0.78690313],
#                          [0.85717921, 0.42895138, 0.54380405, np.nan, 0.8192508 ],
#                          [0.87531807, 0.53611699, 0.78690313, 0.8192508, np.nan]])

# currCoM = np.array([[ -9.09673548,  23.98545742,  25.98642635],
#                     [-43.86865997,  -7.98697188,  41.49483967],
#                     [-38.51215363,  13.03796983,  37.25460958],
#                     [-56.05636501, -22.91237593,   0.32856068],
#                     [-32.29449158,  13.55492001,   2.34797442]])

# LabelNames = ['aCING', 'M1', 'MFC', 'SMTC', 'STRI']
# cRange = [0, 1]

#  # Log %AO of LBD with/without AD
# path_yesAD = np.array([[-5.38159242, -3.81825948, -6.08981196, -7.55075267, np.nan, np.nan],
#                        [np.nan, -2.59552094, -5.56703049, np.nan, np.nan, -3.11550729],
#                        [-2.71794362, -2.92779329, np.nan, -3.76465793, -4.14361923, np.nan],
#                        [-4.66398291, -3.51839434, -5.82119475,         np.nan, -6.33750394, -4.84813386],
#                        [-5.4196171,  -2.63413546, -5.84365942, -6.8656541,  -6.42612516, -5.29154038],
#                        [-4.81822448, -4.45262079, -5.14399051,         np.nan, -4.5354183,  -3.77300047],
#                        [np.nan, -2.96901245, -3.60893702, -5.39989615, -3.69086382, -3.15733398],
#                        [-3.90855403,         np.nan,         np.nan, -6.38339836, -5.01572128, -4.3246638 ],
#                        [-6.05476384,         np.nan,         np.nan, -6.41156582, -7.17782149, -5.3491897 ],
#                        [-3.71195858, -3.59451518, -5.33370817, -4.78568005, -4.33937127,         np.nan],
#                        [-3.97373695, -3.40255831, -6.61403591, -5.10527198, -4.71588674, -3.77178292],
#                        [-3.54693666, -3.42129805, -4.1053954,  -3.99821551, -3.9554539,  -3.51056592],
#                        [-3.43755753, -2.95581949, -6.23870035, -5.30737859, -3.50816349, -4.50713645],
#                        [-3.46577203, -2.98321099, -3.88159335, -5.78933326, -3.68076933, -3.59218841],
#                        [-2.71194339, -3.15858067, -6.04904281, -4.51839304, -3.0252348,  -3.21757667],
#                        [        np.nan, -3.36488474, -5.65441153, -7.33817037, -4.3623129,  -4.17729148],
#                        [-2.43984301, -2.15330339, -4.18278983, -3.02214089, -2.52160035, -2.47725996],
#                        [-4.81094276, -3.99921622,         np.nan, -8.00370446, -4.11561228, -4.85645604],
#                        [-5.11715865, -3.59878284, -5.78764839, -6.61745758, -5.03596637, -5.13127867],
#                        [-3.73242241,         np.nan, -4.715651,   -5.47295651, -4.00003792, -3.67320297]])

# path_noAD = np.array([[-5.0960741,  -5.77990787, -7.66868526,         np.nan,         np.nan,         np.nan],
#                       [-6.7893036,          np.nan, -7.59088613, -8.46299439, -7.88384217, -6.58133338],
#                       [        np.nan, -3.98540048, -7.00409099, -5.88031004, -6.480055,   -4.22810455],
#                       [-5.80145125,         np.nan,         np.nan,         np.nan, -7.99826714,         np.nan],
#                       [-6.21239181, -4.68410517,         np.nan,         np.nan, -6.26443368, -6.20218558],
#                       [        np.nan,         np.nan, -5.21488426, -4.91697211, -4.03055678, -4.3574511 ],
#                       [-3.67446517,         np.nan, -6.11914642, -6.65999007, -4.91341784, -5.30051979],
#                       [-5.55892244,         np.nan, -6.78384679, -7.83750741, -6.73176324, -5.12688825],
#                       [        np.nan,         np.nan, -6.59895404, -7.48336637,         np.nan, -4.88554568],
#                       [-3.72108176, -2.37127087, -4.11769508,         np.nan, -3.24705756,         np.nan],
#                       [-6.38270136, -5.4005716,  -7.68936943, -7.95302256, -7.50880958, -5.79853281],
#                       [-5.65965495,         np.nan, -4.45903833, -6.22522952,         np.nan,         np.nan],
#                       [-4.56559655, -3.99122856, -5.47127449, -5.9103916,  -4.38299186, -3.58311502],
#                       [        np.nan, -5.19576078, -6.22059529,         np.nan,         np.nan, -5.87172787],
#                       [-5.39205841, -3.50223393,         np.nan,         np.nan,         np.nan, -5.03410865],
#                       [-5.90295212, -3.58383499,         np.nan, -7.66985774, -8.05815057, -5.66938105],
#                       [-5.11668464, -3.26838124,         np.nan, -5.58225408, -6.59175911, -5.0635615 ],
#                       [-4.12558576, -3.87598041, -4.83691311,         np.nan, -5.96051974, -4.62333416],
#                       [-4.91935921,         np.nan, -5.44878683,         np.nan, -4.68980111,         np.nan],
#                       [-5.7764028,  -3.87030541, -6.10130354, -7.90903464, -5.82350129, -6.05376557],
#                       [-3.65864863, -3.14595069,         np.nan, -4.90543925, -5.11006896, -4.50387793],
#                       [-5.59907035, -5.11167184, -5.71051476, -6.57538968, -5.35026772, -4.8757985 ],
#                       [-3.95349024,         np.nan, -5.98929418, -6.97662393, -4.86255785, -4.10870721],
#                       [-5.59728123,         np.nan,         np.nan, -8.06423746, -6.07154945,         np.nan]])

# path_yesAD_exp = path_yesAD.copy()
# path_noAD_exp = path_noAD.copy()

# # Get min/max %AO of LBD
# minPath = np.nanmin(np.vstack([path_yesAD_exp, path_noAD_exp]), axis=0)
# maxPath = np.nanmax(np.vstack([path_yesAD_exp, path_noAD_exp]) - minPath + 0.0015, axis=0)

# MakerVecYes = np.nanmean(path_yesAD_exp, axis=0)
# MakerVecYes = 3 * (MakerVecYes - minPath) / maxPath

# MakerVecNo = np.nanmean(path_noAD_exp, axis=0)
# MakerVecNo = 3 * (MakerVecNo - minPath) / maxPath

# MakerVecYes = np.delete(MakerVecYes, 1)
# MakerVecNo = np.delete(MakerVecNo, 1)
 
# pSelectCol = 1

# sn = 6

# colorVec = np.ones(sn-1)


# # Define figure
# fig_atlas = plt.figure()

# pSelectCol = 1
# plotNetwork3(fig_atlas, covMat_yesAD, NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0], currCoM, LabelNames, cRange, MakerVecYes, pSelectCol, colorVec, 3)

# plt.savefig(outputDir + "/TESTATLAS_TAUGTTDP2.png", dpi=400)