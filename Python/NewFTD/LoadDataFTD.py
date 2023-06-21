import os
import pandas as pd
import numpy as np
import sys

sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/InvivoAnalysis/")
from LoadNetworkDataByID import LoadNetworkDataByID
sys.path.insert(1, "/Users/hyung/Research23_Network_Analysis/mBIN/Python/LBD/")
from findPathCoM import findPathCoM


def loaddataFTD(baseDir, dataDir, outputDir, NetworkDataGeneral):
    
    # Make directory for saving outputs
    if (os.path.exists(outputDir)): # if the directory already exists skip
        pass
    else:
        os.mkdir(outputDir)

    # Read Tables as Pandas
    # demoT: MRI Data - Demographics 
    # measureT: clinical test, measurement (Traditional  Neurological tests in measuring disease severity)
    demoT = pd.read_excel(os.path.join(dataDir, 'LBDData', 'MRI Schaeffer Demographics classification.xlsx'))
    measuresT = pd.read_excel(os.path.join(dataDir, 'LBDData', 'LBD neuropsych annualized change no formulas.xlsx')) 


    # ------------------------------------------------------------------------------------PATHOLOGY PART---------------------------------------------------------------------------------------
    # pathT: ex-vivo histopathology Data (Quantification) / %AO
    # FTLD_Library_3-31-22_update --> OLD
    # FTLD Library 4-25-23 update --> NEW!!
    new_pathT = pd.read_excel(os.path.join(dataDir, 'NewFTDData', 'FTLD Library 4-25-23 update.xlsx'), dtype={'INDDID': str, 'Tau1_TDP2': str})
    
    #### Also this has further divided into GM and WM / L and R Hemisphere

    # SlideID, StainingRun,	NumofTiles, Uni.U_Bi.B,	Uni1_Bi2, Old1_New2 --> Columns are missing since there are multiple values for each INDDID
    # AutopsyIDRegion, AutopsyIDHemiRegion, AutopsyIDHemiRegionAnalysisRegion, SlideIDRegion, SlideIDHemiRegion, UniqueID-SlideIDHemiRegionAnalysisRegion --> These columns are just adding {L, R} or {GM, WM}

    # Format the new_pathT to Desired Format - For each INDDID divided into {GM, WM} and Also {L, R} (maximum 4 rows per INDDID)
    pathT_WMGM = pd.pivot_table(new_pathT, values='AvgPercentAO', index=['INDDID', 'FullAutopsyID', 'AutopsyIDNumOnly', 'Tau1_TDP2', 'Hemisphere_by_slide', 'AnalysisRegion'], columns=['Region'], aggfunc=np.sum)

    # Unstacking the Index --> Need a way to solve this without saving to csv format
    pathT_WMGM.to_csv(os.path.join(dataDir, 'NewFTDData', 'new_pathT(GMWM).csv'))
    pathT_WMGM = pd.read_csv(os.path.join(dataDir, 'NewFTDData', 'new_pathT(GMWM).csv'))

    print("Total Unique INNDID in whole dataset")
    print(len(pd.unique(pathT_WMGM['INDDID'])))

    # Divide the pathT into GM and WM (Still divided into {L, R})
    pathT_WMGM_type = pathT_WMGM.groupby('AnalysisRegion')

    # This contains 2 seperate rows for {L, R}
    pathT_GM_LR = pathT_WMGM_type.get_group('GM')
    pathT_WM_LR = pathT_WMGM_type.get_group('WM')
    
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

    # Save pathT GM/WM to csv
    pathT_GM.to_csv(os.path.join(dataDir, 'NewFTDData', 'new_pathT(GM).csv'), index=False)
    pathT_WM.to_csv(os.path.join(dataDir, 'NewFTDData', 'new_pathT(WM).csv'), index=False)

    # Mapping from Anatomical regions to Atlas regions (This could be one to many, ex - aCING -> multiple ACC regions) / Connected to Schaffer Atlas
    pathLUT = pd.read_csv(os.path.join(dataDir,'schaefer_path_20210719_20220328.csv'))

    # Read Lookup table to match anatomical regions in the brain to the Atlas region
    # PathToAtlasLUT_5_10_2023(mePFC_PFC_Ignored) --> NEW!!
    AtlasToPathLUT = pd.read_excel(os.path.join(dataDir,'NewFTDData','PathToAtlasLUT_5_10_2023(mePFC_PFC_Ignored).xlsx'))

    # Get match between pathology (anatomical regions) to CoM and Atlas Index (unordered) / Atlas index 1~400 regions
    pathCoMunordered, pathToAtlasIndexunordered = findPathCoM(pathLUT, AtlasToPathLUT, NetworkDataGeneral['NetworkDataGeneral'][0,0]['Schaefer400x7']['CoM'][0, 0])

    #### Get List of all regions of pathology we can map to 3D Atlas (out of 22) #### CURRENTLY MANUAL
    # THIS HAS TO MATCH THE ORDER OF Pathology Regions in the Pathology Dataset (Left to Right)
    pathNames_3D_Map  = ['ANG', 'ATC', 'HIP', 'IFC', 'M1', 'MFC', 'OFC', 'PC', 'S1', 'SMTC', 'SPC', 'V1', 'aCING', 'aINS', 'aITC', 'dlPFC', 'iPFC', 'mPFC', 'pCING', 'pSTC']
    # sn - denote the number of areas we are able to map to 3D Atlas
    sn = 20
    ########

    # Ordering the CoM so that it matches the order of Regions in the Pathology Dataset (Columns)
    pathCoM = np.empty((sn,3,2)) # One path regions corresponds to multiple atlas region
    pathToAtlasIndex = [[None, None] for _ in range(sn)]

    for s in range(sn):
        idx = AtlasToPathLUT[AtlasToPathLUT.PathSpreadSheetNames == pathNames_3D_Map[s]].index[0] 
        pathCoM[s,:,:] = pathCoMunordered[idx, :, :]
        pathToAtlasIndex[s] = pathToAtlasIndexunordered[idx]

    # pathCoM and pathToAtlasIndex are ordered by the order of pathNames_3D_Map (= Ordering of regions in PathT Dataset Columns Left to Right)

    # ------------------------------------------------------------------------------------THICKNESS PART------------------------------------------------------------------------------------
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

    ########### MEASURE VS THICKNESS PART NOT YET IMPLEMENTED (MISSING DATA) ###########
    
    # numLab
    numLab=400

    # Make directory for proccessedResultsLBD
    if (os.path.exists(os.path.join(dataDir, 'proccessedResultsNewFTD'))): # if the directory already exists skip
        pass
    else:
        os.mkdir(os.path.join(dataDir, 'proccessedResultsNewFTD'))

    # .mat file paths to save
    thicknessTAUSaveFile = os.path.join(dataDir, 'proccessedResultsNewFTD', 'Patient_TAU_thicknessVales.mat')
    thicknessTDPSaveFile = os.path.join(dataDir, 'proccessedResultsNewFTD', 'Patient_TDP_thicknessVales.mat')
    thicknessSaveHCFile = os.path.join(dataDir, 'proccessedResultsNewFTD', 'HC_thicknessVales.mat')

    # Length of HC IDs: 54
    # Length of TAU IDs: 26
    # Length of TDP IDs: 30
    # Get thickness mean and volume total values for [Patient (TAU) MRI data IDs (26) x lables (numLab = 400 regions in the sch region)] / 26 x 400
    PatientTAUResults = LoadNetworkDataByID(TAU_IDs, thicknessPatientTAU, thicknessTAUSaveFile,'Schaefer400x7v1', numLab = 400, ICV = False)

    # Get thickness mean and volume total values for [Patient (TDP) MRI data IDs (30) x lables (numLab = 400 regions in the sch region)] / 26 x 400
    PatientTDPResults = LoadNetworkDataByID(TDP_IDs, thicknessPatientTDP, thicknessTDPSaveFile,'Schaefer400x7v1', numLab = 400, ICV = False)  # NOT IMPLEMENTED SOME ISSUE

    # Get thickness mean and volume total values for [Control MRI data IDs (54) x lables (numLab)] 26 x 400
    HCResults = LoadNetworkDataByID(HC_IDs, thicknessHC, thicknessSaveHCFile,'Schaefer400x7v1', numLab = 400, ICV = False)

    assert(PatientTAUResults['Thickness']['Mean'].shape == (26, 400))
    assert(PatientTDPResults['Thickness']['Mean'].shape == (30, 400))
    assert(HCResults['Thickness']['Mean'].shape == (54, 400))
    

    # --------------------------------------MATCHING THICKNESS DATA SET TO PATH DATASET--------------------------------------
    # Matching Atlas(in-vivo) to pathology data

    # pathToAtlasIndex --> 20 x 2 (sn)

    # N = 75 (AllResults.Thickness.Mean.shape -> 75 x 400) Number of Patients of MRI Data
    n_TAU = PatientTAUResults['Thickness']['Mean'].shape[0]
    n_TDP = PatientTDPResults['Thickness']['Mean'].shape[0]
    n_HC = HCResults['Thickness']['Mean'].shape[0]

    # sn = 20 Log %AO of anatomical regions

    # Shape --> N x 20 x 2 / N - Number of subjects, 20 - Nunber of Path Regions, 2 - {L, R}
    # N --> Ordered by the 
    TAUthicknessAtPath = np.empty((n_TAU,sn,2))
    TDPthicknessAtPath = np.empty((n_TDP,sn,2))
    HCthicknessAtPath = np.empty((n_HC,sn,2))

    TAUvolumeAtPath = np.empty((n_TAU,sn,2))
    TDPvolumeAtPath = np.empty((n_TDP,sn,2))
    HCvolumeAtPath = np.empty((n_HC,sn,2))

    # Weighted Mean (NEED TO IMPLEMENT)
    TAUthicknessAtPath_Weighted = np.empty((n_TAU,sn,2))
    TDPthicknessAtPath_Weighted = np.empty((n_TDP,sn,2))
    HCthicknessAtPath_Weighted = np.empty((n_HC,sn,2))

    for n in range(n_TAU):
        for p in range(sn):
            for r in range(2):
                curIdx = pathToAtlasIndex[p][r]
                TAUthicknessAtPath[n,p,r] = np.mean(PatientTAUResults['Thickness']['Mean'][n,curIdx])
                TAUvolumeAtPath[n,p,r] = np.mean(PatientTAUResults['Volume']['Normalized'][n,curIdx]) # Normalized Volume
    
    for n in range(n_TDP):
        for p in range(sn):
            for r in range(2):
                curIdx = pathToAtlasIndex[p][r]
                TDPthicknessAtPath[n,p,r] = np.mean(PatientTDPResults['Thickness']['Mean'][n,curIdx])
                TDPvolumeAtPath[n,p,r] = np.mean(PatientTDPResults['Volume']['Normalized'][n,curIdx])
    
    for n in range(n_HC):
        for p in range(sn):
            for r in range(2):
                curIdx = pathToAtlasIndex[p][r]
                HCthicknessAtPath[n,p,r] = np.mean(HCResults['Thickness']['Mean'][n,curIdx])
                HCvolumeAtPath[n,p,r] = np.mean(HCResults['Volume']['Normalized'][n,curIdx])

    return pathCoM, pathT_GM, pathT_WM, pathNames_3D_Map, sn, pathLUT, HCResults, PatientTAUResults, PatientTDPResults, TAUthicknessAtPath, TDPthicknessAtPath, HCthicknessAtPath, TAUvolumeAtPath, TDPvolumeAtPath, HCvolumeAtPath