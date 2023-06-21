import numpy as np
from scipy.io import savemat

# CtrlResults = LoadNetworkDataByID(ctrlIDs, thicknessCtrlraw, thicknessSaveFileCtrl,'schaefer400x7', 400, ICV = False)
# Have a general dataloader with data types (9 types)
def LoadNetworkDataByID(ID, allThickVals, outMatSaveFile, atlasName, numLab, ICV):
    """
    A general dataloader to load Thickness (Mean, Median) / Volume (Total, ICV) values for [IDs x numLab]

    Attributes
    ----------
    ID : list of int
        list of INDDID for [Control/Patients]
    allThickVals : 
        Thickness and Volumes from MRI Data of the Cortex [Control/Patients]
    outMatSaveFile : 
        .mat file to save the loaded values, Thickness (Mean, Median) / Volume (Total, ICV) values for [IDs x numLab]
    atlasName : 
        Atlas Name
    numLab : 
        number of labels (regions?)
    ICV : # FOR THE LewyBody / The data is bit different from FTD --> Different way of loading the dataset!!
        Boolean value denoting if we want to load ICV (IntraCranial Volume)
    
    Returns
    ----------
    the loaded values, Thickness (Mean, Median) / Volume (Total, ICV) values for [IDs x numLab]
    """

    # Length of IDs List
    N = len(ID)

    # AllResults Dictionary
    AllResults = {}

    AllResults['Age'] = np.full((N), np.nan)
    AllResults['Sex'] = [None for _ in range(N)]

    AllResults['Thickness'] = {}
    AllResults['Thickness']['Mean'] = np.full((N, numLab), np.nan)
    AllResults['Thickness']['Median'] = np.full((N, numLab), np.nan)

    AllResults['Volume'] = {}
    AllResults['Volume']['Total'] = np.full((N, numLab), np.nan)
    AllResults['Volume']['ICV'] = np.full((N), np.nan)
    AllResults['Volume']['Normalized'] = np.full((N, numLab), np.nan)

    if(ICV):
        AllResults['Volume']['ICV'] = np.full((N, numLab), np.nan)

    for i in range(N):
        print(i+1)

        currID = ID[i]
        currDates = allThickVals[allThickVals['id'] == currID]['date']
        uniqDates = np.unique(currDates) # Date is to specify a image in a sequence of MRI images --> Need a way to specify this part. Currently the first date is used.

        print(str(currID))

        if (not ICV): # List of length 400
            currThickness = allThickVals[(allThickVals['id'] == currID) & 
                                        (allThickVals['date'] == uniqDates[0]) & 
                                        (allThickVals['system'] == atlasName) & 
                                        (allThickVals['measure'] == 'thickness') & 
                                        (allThickVals['metric'] == 'mean')]['value']

            currVolume = allThickVals[(allThickVals['id'] == currID) & 
                                  (allThickVals['date'] == uniqDates[0]) & 
                                  (allThickVals['system'] == atlasName) & 
                                  (allThickVals['measure'] == 'volume') & 
                                  (allThickVals['metric'] == 'numeric')]['value']
            
            currICV = allThickVals[(allThickVals['id'] == currID) & 
                                  (allThickVals['date'] == uniqDates[0]) & 
                                  (allThickVals['system'] == 'AntsBrain') & 
                                  (allThickVals['measure'] == 'volume') & 
                                  (allThickVals['metric'] == 'numeric')]['value']
            
            AllResults['Volume']['ICV'][i] = currICV.tolist()[0]

            currAge = allThickVals[(allThickVals['id'] == currID) & 
                                        (allThickVals['date'] == uniqDates[0]) & 
                                        (allThickVals['system'] == atlasName) & 
                                        (allThickVals['measure'] == 'thickness') & 
                                        (allThickVals['metric'] == 'mean')]['AgeAtMR']
            
            # Sanity Check, that Age for IDs is unique.
            assert(currAge.shape == (numLab,))
            assert(len(np.unique(currAge)) == 1)
            
            AllResults['Age'][i] = currAge.tolist()[0]

            currSex = allThickVals[(allThickVals['id'] == currID) & 
                                        (allThickVals['date'] == uniqDates[0]) & 
                                        (allThickVals['system'] == atlasName) & 
                                        (allThickVals['measure'] == 'thickness') & 
                                        (allThickVals['metric'] == 'mean')]['Sex']
            
            # Sanity Check, that Sex for IDs is unique.
            assert(currSex.shape == (numLab,))
            assert(len(np.unique(currSex)) == 1)

            AllResults['Sex'][i] = currSex.tolist()[0]

            
            # Label index - for 400 regions of the brain
            currLabelsThick = allThickVals[(allThickVals['id'] == currID) & 
                                       (allThickVals['date'] == uniqDates[0]) & 
                                       (allThickVals['system'] == atlasName) & 
                                       (allThickVals['measure'] == 'thickness') & 
                                       (allThickVals['metric'] == 'mean')]['label']
            
            currLabelsVol = allThickVals[(allThickVals['id'] == currID) & 
                                     (allThickVals['date'] == uniqDates[0]) & 
                                     (allThickVals['system'] == atlasName) & 
                                     (allThickVals['measure'] == 'volume') & 
                                     (allThickVals['metric'] == 'numeric')]['label']      
            
            for L in range(numLab):
                if (np.sum(currLabelsThick == L+1) == 1): # 2 ways to fail / 1) couldn't find it 2) Multiple matches. 
                    if (len(currThickness) > 0):
                        AllResults['Thickness']['Mean'][i, L] = currThickness[currLabelsThick == L+1]
                else:
                    print('Missing Thickness Label ' + str(L+1) + ' ID ' + str(currID))

                if (np.sum(currLabelsVol == L+1) == 1):
                    if (len(currVolume) > 0):
                        AllResults['Volume']['Total'][i, L] = currVolume[currLabelsVol == L+1]
                else:
                    print('Missing Volume Label ' + str(L+1) + ' ID ' + str(currID))
        
        if (ICV):
            currVolume = allThickVals[(allThickVals['id'] == currID) & 
                                  (allThickVals['date'] == uniqDates[0]) & 
                                  (allThickVals['system'] == atlasName)]['volume']
            
            currICV = allThickVals[(allThickVals['id'] == currID) & 
                                  (allThickVals['date'] == uniqDates[0]) & 
                                  (allThickVals['system'] == atlasName)]['icv']
            
            currLabelsVol = allThickVals[(allThickVals['id'] == currID) & 
                                     (allThickVals['date'] == uniqDates[0]) & 
                                     (allThickVals['system'] == atlasName)]['label']
            
            for L in range(numLab):
                if (np.sum(currLabelsVol == L+1) == 1):
                    if (len(currVolume) > 0):
                        AllResults['Volume']['Total'][i, L] = currVolume[currLabelsVol == L+1]
                        AllResults['Volume']['ICV'][i, L] = currICV[currLabelsVol == L+1]
                else:
                    print('Missing Volume Label ' + str(L+1) + ' ID ' + str(currID))
        
    # Get the Normalized Volume value (Volume / ICV) per each subject
    AllResults['Volume']['Normalized'] = AllResults['Volume']['Total'] / AllResults['Volume']['ICV'][:, None]
    
    
    savemat(outMatSaveFile, AllResults)
    return AllResults