import os
import pandas as pd
import numpy as np
import sys

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/InvivoAnalysis/")
from LoadNetworkDataByID import LoadNetworkDataByID
from findPathCoM import findPathCoM

def pathLoadData(dataDir, NetworkDataGeneral):
    """This method will load Pathology Data to our desired format

    Args:
        dataDir (str): Path location of the data
        NetworkDataGeneral (obj): Preconstructed atlas data
    
    Returns:
        pathT_GM (DataFrame): Pathology %AO data for GM
        pathT_WM (DataFrame): Pathology %AO data for WM
        pathLUT (DataFrame): Look up table matching Atlas region names to Atlas Labels (Index) in NetworkDataGeneral Object
        sn (int): Number of areas we are able to map to 3D Atlas
        pathCoM (ndarray): Center of Mass of Pathology Regions (Mean of multiple corresponding Atlas regions' Center of Mass)
        pathToAtlasIndex (ndarray of list): ndarray of List of Atlas regions' index corrresponding to Pathology regions (this could be multiple Atlas regions)
    """

    ####################################################################### Loading Pathology %AO #######################################################################
    # Load new_pathT: ex-vivo histopathology Data (Quantification) / %AO for pathology regions
    new_pathT = pd.read_excel(os.path.join(dataDir, 'NewFTDData', 'FTLD Library 4-25-23 update.xlsx'), dtype={'INDDID': str, 'Tau1_TDP2': str})
    
    #### Some columns are not included #### --> Need to include this
        # SlideID, StainingRun,	NumofTiles, Uni.U_Bi.B,	Uni1_Bi2, Old1_New2 --> Columns are missing since there are multiple values for each INDDID
        # AutopsyIDRegion, AutopsyIDHemiRegion, AutopsyIDHemiRegionAnalysisRegion, SlideIDRegion, SlideIDHemiRegion, UniqueID-SlideIDHemiRegionAnalysisRegion --> These columns are just adding {L, R} or {GM, WM} (therefore columns not added)

    # Format the new_pathT to Desired Format - For each INDDID divided into {GM, WM} and {L, R} (maximum 4 rows per INDDID)
    pathT_WMGM = pd.pivot_table(new_pathT, values='AvgPercentAO', index=['INDDID', 'FullAutopsyID', 'AutopsyIDNumOnly', 'Tau1_TDP2', 'Hemisphere_by_slide', 'AnalysisRegion'], columns=['Region'], aggfunc=np.sum)

    # Unstacking the Index --> Need a way to solve this without saving to csv format
    pathT_WMGM.to_csv(os.path.join(dataDir, 'NewFTDData', 'new_pathT(GMWM).csv'))
    pathT_WMGM = pd.read_csv(os.path.join(dataDir, 'NewFTDData', 'new_pathT(GMWM).csv'))

    # Divide the pathT into GM and WM (Still divided into {L, R})
    pathT_WMGM_type = pathT_WMGM.groupby('AnalysisRegion')

    # This contains 2 seperate rows for {L, R}
    pathT_GM_LR = pathT_WMGM_type.get_group('GM')
    pathT_WM_LR = pathT_WMGM_type.get_group('WM')
    
    print("Total Unique INNDID in whole dataset")
    print(len(pd.unique(pathT_WMGM['INDDID'])))
    print("Unique INDDID in GM")
    print(len(pd.unique(pathT_GM_LR['INDDID'])))
    print("Unique INDDID in WM")
    print(len(pd.unique(pathT_WM_LR['INDDID'])))

    # Combine 2 Rows for {L, R} into a single row
    pathT_GM_LR_type = pathT_GM_LR.groupby('Hemisphere_by_slide')
    pathT_GM_L = pathT_GM_LR_type.get_group('L')
    pathT_GM_R = pathT_GM_LR_type.get_group('R')
    pathT_GM = pd.merge(pathT_GM_L, pathT_GM_R, left_on=['INDDID', 'FullAutopsyID', 'AutopsyIDNumOnly', 'Tau1_TDP2', 'AnalysisRegion'], right_on=['INDDID', 'FullAutopsyID', 'AutopsyIDNumOnly', 'Tau1_TDP2', 'AnalysisRegion'], how='outer', suffixes=('_L', '_R')) 

    pathT_WM_LR_type = pathT_WM_LR.groupby('Hemisphere_by_slide')
    pathT_WM_L = pathT_WM_LR_type.get_group('L')
    pathT_WM_R = pathT_WM_LR_type.get_group('R')
    pathT_WM = pd.merge(pathT_WM_L, pathT_WM_R, left_on=['INDDID', 'FullAutopsyID', 'AutopsyIDNumOnly', 'Tau1_TDP2', 'AnalysisRegion'], right_on=['INDDID', 'FullAutopsyID', 'AutopsyIDNumOnly', 'Tau1_TDP2', 'AnalysisRegion'], how='outer', suffixes=('_L', '_R'))

    # Drop Hemisphere_by_slide {L, R} Columns
    pathT_GM = pathT_GM.drop(columns=['Hemisphere_by_slide_L', 'Hemisphere_by_slide_R'])
    pathT_WM = pathT_WM.drop(columns=['Hemisphere_by_slide_L', 'Hemisphere_by_slide_R']) 

    #### [pathT_GM / pathT_WM] The Region Columns are in alphabetical order (from left to right) ####
    
    ####################################################################### Mapping Pathology Regions to Atlas regions #######################################################################
    # Load the Look up table matching Atlas region names to Atlas Labels(Index)
    pathLUT = pd.read_csv(os.path.join(dataDir,'schaefer_path_20210719_20220328.csv'))

    # Load the Look up table matching Pathology region names to Atlas region names
    AtlasToPathLUT = pd.read_excel(os.path.join(dataDir,'NewFTDData','PathToAtlasLUT_5_10_2023(mePFC_PFC_Ignored).xlsx'))

    # Get match between pathology (pathology regions) to CoM and Atlas Index (unordered) / Atlas index 1~400 regions
    pathCoMunordered, pathToAtlasIndexunordered = findPathCoM(pathLUT, AtlasToPathLUT, NetworkDataGeneral['NetworkDataGeneral'][0,0]['Schaefer400x7']['CoM'][0, 0])

    #### Get List of all regions of pathology we can map to 3D Atlas (out of 22) #### 
    # THIS HAS TO MATCH THE ORDER OF Pathology Regions in the Pathology Dataset (Alphabetical order from Left to Right) 
    # ['ANG', 'ATC', 'HIP', 'IFC', 'M1', 'MFC', 'OFC', 'PC', 'S1', 'SMTC', 'SPC', 'V1', 'aCING', 'aINS', 'aITC', 'dlPFC', 'iPFC', 'mPFC', 'pCING', 'pSTC']
    pathNames_3D_Map = np.sort(AtlasToPathLUT["PathSpreadSheetNames"].values)
    
    # sn - denote the number of areas we are able to map to 3D Atlas
    sn = len(pathNames_3D_Map)

    # Ordering the CoM so that it matches the order of Regions in the Pathology Dataset (Columns)
    pathCoM = np.empty((sn,3,2)) # One path regions corresponds to multiple atlas region
    pathToAtlasIndex = [[None, None] for _ in range(sn)]

    for s in range(sn):
        idx = AtlasToPathLUT[AtlasToPathLUT.PathSpreadSheetNames == pathNames_3D_Map[s]].index[0] 
        pathCoM[s,:,:] = pathCoMunordered[idx, :, :]
        pathToAtlasIndex[s] = pathToAtlasIndexunordered[idx]

    # pathCoM and pathToAtlasIndex are ordered by the order of pathNames_3D_Map (= Ordering of regions same as in PathT Dataset Columns Left to Right)

    # Drop Columns in pathT_GM / pathT_GM Where we cannot map to 3D Atlas, using AtlasToPathLUT (+5, for index offset)
    pathT_GM = pathT_GM.drop(pathT_GM.columns[[i + 5 for i, e in enumerate(pathT_GM.columns.values[5:]) if e.split("_")[0] not in pathNames_3D_Map]], axis = 1)
    pathT_WM = pathT_WM.drop(pathT_WM.columns[[i + 5 for i, e in enumerate(pathT_WM.columns.values[5:]) if e.split("_")[0] not in pathNames_3D_Map]], axis = 1)

    # Save pathT GM/WM to csv
    pathT_GM.to_csv(os.path.join(dataDir, 'NewFTDData', 'new_pathT(GM).csv'), index=False)
    pathT_WM.to_csv(os.path.join(dataDir, 'NewFTDData', 'new_pathT(WM).csv'), index=False)

    return pathT_GM, pathT_WM, pathLUT, sn, pathCoM, pathToAtlasIndex

def thickLoadData(dataDir):
    """This method will load Thickness Data to our desired format. Divided into 3 parts. Healthy Control, TAU, TDP

    Args:
        dataDir (str): Path location of the data

    Returns:
        HCResults (dict): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of Healthy Control Subjects
        PatientTAUResults (dict): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of TAU Patients
        PatientTDPResults (dict): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of TDP Patients
    """
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
    thicknessPatientTAU = thickness_path_type.get_group('tau')
    # MRI Thickness values for Patient (TDP)
    thicknessPatientTDP = thickness_path_type.get_group('tdp')

    # IDs
    HC_IDs = np.unique(thicknessHC.INDDID)
    TAU_IDs = np.unique(thicknessPatientTAU.INDDID)
    TDP_IDs = np.unique(thicknessPatientTDP.INDDID)

    print("Length of HC IDs in Thickness: " + str(len(HC_IDs)))
    print("Length of TAU IDs in Thickness: " + str(len(TAU_IDs)))
    print("Length of TDP IDs in Thickness: " + str(len(TDP_IDs)))
    
    # numLab (Number of Label in Schaefer400x7 Atlas)
    numLab=400

    # Make directory for proccessedResultsLBD
    if (os.path.exists(os.path.join(dataDir, 'proccessedResultsNewFTD'))): # if the directory already exists skip
        pass
    else:
        os.mkdir(os.path.join(dataDir, 'proccessedResultsNewFTD'))

    # .mat file paths to save
    thicknessSaveHCFile = os.path.join(dataDir, 'proccessedResultsNewFTD', 'HC_thicknessVales.mat')
    thicknessTAUSaveFile = os.path.join(dataDir, 'proccessedResultsNewFTD', 'Patient_TAU_thicknessVales.mat')
    thicknessTDPSaveFile = os.path.join(dataDir, 'proccessedResultsNewFTD', 'Patient_TDP_thicknessVales.mat')
    
    # Length of HC IDs: 54
    # Length of TAU IDs: 26
    # Length of TDP IDs: 30

    # Get thickness mean and volume total values for [Control MRI data IDs (54) x lables (numLab)] 26 x 400
    HCResults = LoadNetworkDataByID(HC_IDs, thicknessHC, thicknessSaveHCFile,'Schaefer400x7v1', numLab = 400)

    # Get thickness mean and volume total values for [Patient (TAU) MRI data IDs (26) x lables (numLab = 400 regions in the sch region)] / 26 x 400
    PatientTAUResults = LoadNetworkDataByID(TAU_IDs, thicknessPatientTAU, thicknessTAUSaveFile,'Schaefer400x7v1', numLab = 400)

    # Get thickness mean and volume total values for [Patient (TDP) MRI data IDs (30) x lables (numLab = 400 regions in the sch region)] / 26 x 400
    PatientTDPResults = LoadNetworkDataByID(TDP_IDs, thicknessPatientTDP, thicknessTDPSaveFile,'Schaefer400x7v1', numLab = 400)

    # Sanity Check
    assert(PatientTAUResults['Thickness']['Mean'].shape == (26, 400))
    assert(PatientTDPResults['Thickness']['Mean'].shape == (30, 400))
    assert(HCResults['Thickness']['Mean'].shape == (54, 400))

    return HCResults, PatientTAUResults, PatientTDPResults

def thickAtPathLoadData(HCResults, PatientTAUResults, PatientTDPResults, sn, pathToAtlasIndex):
    """_summary_

    Args:
        HCResults (dict): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of Healthy Control Subjects
        PatientTAUResults (dict): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of TAU Patients
        PatientTDPResults (dict): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of TDP Patients
        sn (int): Number of areas we are able to map to 3D Atlas
        pathToAtlasIndex (list): List of Atlas regions index corrresponding to Pathology regions (this could be multiple Atlas regions per Pathology Region)

    Returns:
        HCthicknessAtPath (ndarray): Thickness values of Healthy Control Subjects Matching Pathology Regions
        TAUthicknessAtPath (ndarray): Thickness values of TAU Patient Subjects Matching Pathology Regions
        TDPthicknessAtPath (ndarray): Thickness values of TDP Patient Subjects Matching Pathology Regions
        HCnormVolumeAtPath (ndarray): Normalized Volume (by ICV) values of Healthy Control Subjects Matching Pathology Regions
        TAUnormVolumeAtPath (ndarray): Normalized Volume (by ICV) values of TAU Patient Subjects Matching Pathology Regions
        TDPnormVolumeAtPath (ndarray): Normalized Volume (by ICV) values of TDP Patient Subjects Matching Pathology Regions
    """
    #### Matching Atlas(in-vivo) to pathology data ####

    # Number of subject for HC, TAU, TDP Data
    n_HC = HCResults['Thickness']['Mean'].shape[0]
    n_TAU = PatientTAUResults['Thickness']['Mean'].shape[0]
    n_TDP = PatientTDPResults['Thickness']['Mean'].shape[0]

    # Shape --> N x sn x 2 / N - Number of subjects, sn - Nunber of Path Regions, 2 - {L, R}
    HCthicknessAtPath = np.empty((n_HC,sn,2))
    TAUthicknessAtPath = np.empty((n_TAU,sn,2))
    TDPthicknessAtPath = np.empty((n_TDP,sn,2))
    
    HCnormVolumeAtPath = np.empty((n_HC,sn,2))
    TAUnormVolumeAtPath = np.empty((n_TAU,sn,2))
    TDPnormVolumeAtPath = np.empty((n_TDP,sn,2))
    
    # Weighted Mean (NEED TO IMPLEMENT)
    HCthicknessAtPath_Weighted = np.empty((n_HC,sn,2))
    TAUthicknessAtPath_Weighted = np.empty((n_TAU,sn,2))
    TDPthicknessAtPath_Weighted = np.empty((n_TDP,sn,2))
    
    for n in range(n_HC):
        for p in range(sn):
            for r in range(2):
                curIdx = pathToAtlasIndex[p][r]
                HCthicknessAtPath[n,p,r] = np.mean(HCResults['Thickness']['Mean'][n,curIdx])
                HCnormVolumeAtPath[n,p,r] = np.mean(HCResults['Volume']['Normalized'][n,curIdx])

    for n in range(n_TAU):
        for p in range(sn):
            for r in range(2):
                curIdx = pathToAtlasIndex[p][r]
                TAUthicknessAtPath[n,p,r] = np.mean(PatientTAUResults['Thickness']['Mean'][n,curIdx])
                TAUnormVolumeAtPath[n,p,r] = np.mean(PatientTAUResults['Volume']['Normalized'][n,curIdx]) # Normalized Volume
    
    for n in range(n_TDP):
        for p in range(sn):
            for r in range(2):
                curIdx = pathToAtlasIndex[p][r]
                TDPthicknessAtPath[n,p,r] = np.mean(PatientTDPResults['Thickness']['Mean'][n,curIdx])
                TDPnormVolumeAtPath[n,p,r] = np.mean(PatientTDPResults['Volume']['Normalized'][n,curIdx])
    
    return HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath

# Main Function Loop
def loaddataFTD(baseDir, dataDir, outputDir, NetworkDataGeneral):
    """This method will be used to load Pathology, Thickness, Thickness at Pathology region data.
    Calls 3 functions sequentially. pathLoadData(), thickLoadData(), thickAtPathLoadData()

    Args:
        baseDir (str): Path location of the github respository
        dataDir (str): Path location of the data
        outputDir (str): Path location to save the analysis results
        NetworkDataGeneral (obj): Preconstructed atlas data

    Returns:
        pathCoM (ndarray): Center of Mass of Pathology Regions (Mean of multiple corresponding Atlas regions' Center of Mass)
        pathT_GM (DataFrame): Pathology %AO data for GM
        pathT_WM (DataFrame): Pathology %AO data for WM
        sn (int): Number of areas we are able to map to 3D Atlas / 20 
        pathLUT (DataFrame): Look up table matching Atlas region names to Atlas Labels (Index) in NetworkDataGeneral Object
        HCResults (dict of ndarrays): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of Healthy Control Subjects
        PatientTAUResults (dict of ndarrays): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of TAU Patients
        PatientTDPResults (dict of ndarrays): Dictionary containing Thickness (Mean, Median), Volume (Total, ICV, Normalized), Age, Sex of TDP Patients
        HCthicknessAtPath (ndarray): Thickness values of Healthy Control Subjects Matching Pathology Regions
        TAUthicknessAtPath (ndarray): Thickness values of TAU Patient Subjects Matching Pathology Regions
        TDPthicknessAtPath (ndarray): Thickness values of TDP Patient Subjects Matching Pathology Regions
        HCnormVolumeAtPath (ndarray): Normalized Volume (by ICV) values of Healthy Control Subjects Matching Pathology Regions
        TAUnormVolumeAtPath (ndarray): Normalized Volume (by ICV) values of TAU Patient Subjects Matching Pathology Regions
        TDPnormVolumeAtPath (ndarray): Normalized Volume (by ICV) values of TDP Patient Subjects Matching Pathology Regions
    """

    # Make directory for saving outputs
    if (os.path.exists(outputDir)): # if the directory already exists skip
        pass
    else:
        os.mkdir(outputDir)

    # PATHOLOGY PART
    pathT_GM, pathT_WM, pathLUT, sn, pathCoM, pathToAtlasIndex = pathLoadData(dataDir, NetworkDataGeneral)

    # THICKNESS PART
    HCResults, PatientTAUResults, PatientTDPResults = thickLoadData(dataDir)
    
    # THICKNESS DATA MATCHING TO PATHOLOGY DATASET
    HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath = thickAtPathLoadData(HCResults, PatientTAUResults, PatientTDPResults, sn, pathToAtlasIndex)

    return pathCoM, pathT_GM, pathT_WM, pathLUT, HCResults, PatientTAUResults, PatientTDPResults, HCthicknessAtPath, TAUthicknessAtPath, TDPthicknessAtPath, HCnormVolumeAtPath, TAUnormVolumeAtPath, TDPnormVolumeAtPath