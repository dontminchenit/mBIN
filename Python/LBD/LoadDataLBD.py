import os
import pandas as pd
import numpy as np
import sys

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/InvivoAnalysis/")
from LoadNetworkDataByID import LoadNetworkDataByID
sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/LBD/")
from findPathCoM import findPathCoM

def loaddataLBD(baseDir, dataDir, outputDir, NetworkDataGeneral):
    
    # Make directory for saving outputs
    if (os.path.exists(outputDir)): # if the directory already exists skip
        pass
    else:
        os.mkdir(outputDir)

    # Read Tables as Pandas
    # demoT: MRI Data - Demographics 
    # pathT: ex-vivo histopathology Data (Quantification) / %AO
    # measureT: clinical test, measurement (Traditional  Neurological tests in measuring disease severity)
    demoT = pd.read_excel(os.path.join(dataDir, 'LBDData', 'MRI Schaeffer Demographics classification.xlsx'))
    pathT = pd.read_csv(os.path.join(dataDir, 'LBDData', 'LBD_AO_aSYN.csv'))
    measuresT = pd.read_excel(os.path.join(dataDir, 'LBDData', 'LBD neuropsych annualized change no formulas.xlsx')) 

    # Index for the case with AD or No Ad for patients with LBD / Pathology Data
    LBD_yesADIndx = (pathT.Group == 'AD+LBD') & (pathT.CDx != 'AD') # False or True
    LBD_noADIndx = (pathT.Group == 'LBD') # False or True

    #------------------------------------------------PATHOLOGY PART-------------------------------------------------------------
    # Get Log %AO of 6 anatomical regions of the brain
    # aCING, ... --> Name of regions (anatomical region) that does not directly match with Atlas. Therefore we match them later.
    pathDataGM = np.log(0.01 * pathT.iloc[:, 1:7].values + 0.00015)

    # Get names of 6 regions of the brain as list.
    pathNamesRaw = list(pathT.columns[1:7])
    sn = len(pathNamesRaw)

    # Mapping from Anatomical regions to Atlas regions (This could be one to many, ex-aCING  multiple ACC regions)
    pathLUT = pd.read_csv(os.path.join(dataDir,'schaefer_path_20210719_20220328.csv'))

    # Read Lookup table to match anatomical regions in the brain to the Atlas region
    AtlasToPathLUT = pd.read_excel(os.path.join(dataDir,'LBDData','PathToAtlasLUT_10_7_2022.xlsx'))

    # Get match between pathology (anatomical regions) to CoM and Atlas Index (unordered)
    pathCoMunordered, pathToAtlasIndexunordered = findPathCoM(pathLUT, AtlasToPathLUT, NetworkDataGeneral['NetworkDataGeneral'][0,0]['Schaefer400x7']['CoM'][0, 0])

    # Ordering the match between pathology (anatomical regions) to CoM and Atlas Index
    pathCoM = np.empty((sn,3,2)) # One path regions corresponds to multiple atlas region
    pathToAtlasIndex = [[None, None] for _ in range(sn)]

    for s in range(sn):
        idx = AtlasToPathLUT[AtlasToPathLUT.PathSpreadSheetNames == pathNamesRaw[s]].index[0]
        pathCoM[s,:,:] = pathCoMunordered[idx, :, :]
        pathToAtlasIndex[s] = pathToAtlasIndexunordered[idx]

    #-------------------------------------------------THICKNESS PART---------------------------------------------------------
    # MRI Thickness value for Control Cortical Atlas - schaefer400x7
    thicknessCtrlraw = pd.read_csv(os.path.join(dataDir, 'LBDData', 'ctrl_4007.csv'))
    
    # Get unique IDs of Healthy Controls
    ctrlIDs = np.unique(thicknessCtrlraw['id'])

    # list of demoT.INDDID - MRI Demographics [Patients]
    thickIDs = demoT.INDDID

    # MRI Thickness and Volume values for PATIENTS Cortical Atlas 
    thicknessAllraw = pd.read_csv(os.path.join(dataDir, 'LBDData', 'ftdc_dlb_ctx_volume.csv'))

    # MRI Thickness and Volume values - Subcortex (Not Used currently)
    volAllraw = pd.read_csv(os.path.join(dataDir, 'LBDData', 'ftdc_dlb_tian_s1.csv'))

    # idx: 1 - if thickIDs are in measuresT.INDDID / 0 - otherwise
    # idx2: The index where thickIDs are in measuresT.INDDID (there is a offset of 1 since python index starts from 0)
    idx = thickIDs.isin(measuresT.INDDID)
    idx2 = pd.Index(measuresT.INDDID).get_indexer(thickIDs[idx])
    idx = idx.astype('uint8') # Boolean to Int

    # Get the measuresT Dataframe with only the rows with match (demoT.INDDID are in measuresT.INDDID)
    measuresT_matched = measuresT.iloc[idx2]

    # list of pathT.INDDID (Common ID of the subject specific)
    AllIDs = pathT.INDDID

    # numLab
    numLab=400

    # Make directory for proccessedResultsLBD
    if (os.path.exists(os.path.join(dataDir, 'proccessedResultsLBD'))): # if the directory already exists skip
        pass
    else:
        os.mkdir(os.path.join(dataDir, 'proccessedResultsLBD'))

    # .mat file paths to save
    thicknessSaveFile = os.path.join(dataDir, 'proccessedResultsLBD', 'thicknessVales.mat')
    subVolSaveFile = os.path.join(dataDir, 'proccessedResultsLBD', 'subVolSaveFile.mat')
    thicknessSaveFileCtrl = os.path.join(dataDir, 'proccessedResultsLBD', 'thicknessValesCtrl.mat')
    subVolSaveFileCtrl = os.path.join(dataDir, 'proccessedResultsLBD', 'subVolSaveFileCtrl.mat')

    # Get thickness mean and volume total values for [Patient MRI data IDs x labels (numLab = 400 regions in the sch region)]
    AllResults = LoadNetworkDataByID(thickIDs, thicknessAllraw, thicknessSaveFile,'Schaefer400x7v1', 400, ICV = False)

    # Get Volume Total and ICV values for [Patient MRI data IDs x labels (numLab)] - SubCortex
    SubCorticalVol = LoadNetworkDataByID(thickIDs, volAllraw, subVolSaveFile,'Tian_Subcortex_S1_3T', 16, ICV = True)

    # Get thickness mean and volume total values for [Control MRI data IDs x labels (numLab)]
    CtrlResults = LoadNetworkDataByID(ctrlIDs, thicknessCtrlraw, thicknessSaveFileCtrl,'schaefer400x7', 400, ICV = False)

    #-------------------------------------------------------------------------------
    
    # Matching Atlas(in-vivo) to pathology data
    # N = 75 (AllResults.Thickness.Mean.shape -> 75 x 400) Number of Patients of MRI Data
    N = AllResults['Thickness']['Mean'].shape[0]
    # PN = 6 (pathDataGM.shape -> 59 x 6) Log %AO of anatomical regions
    PN = pathDataGM.shape[1]

    AllthicknessAtPath = np.empty((N,PN,2))
    CtrlthicknessAtPath = np.empty((N,PN,2))

    for n in range(N):
        for p in range(PN):
            for r in range(2):
                curIdx = pathToAtlasIndex[p][r]
                AllthicknessAtPath[n,p,r] = np.mean(AllResults['Thickness']['Mean'][n,curIdx])
                CtrlthicknessAtPath[n,p,r] = np.mean(CtrlResults['Thickness']['Mean'][n,curIdx])

    return pathCoM, pathNamesRaw, pathDataGM, LBD_yesADIndx, LBD_noADIndx, sn, pathLUT, CtrlResults, AllResults, demoT, measuresT_matched