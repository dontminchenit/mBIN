import numpy as np
from scipy.io import savemat

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
    ICV : 
        Boolean value denoting what kind of Network Data we are trying to load
    
    Returns
    ----------
    the loaded values, Thickness (Mean, Median) / Volume (Total, ICV) values for [IDs x numLab]
    """

    # Length of IDs List
    N = len(ID)

    # AllResults Dictionary
    AllResults = {}

    AllResults['Thickness'] = {}
    AllResults['Thickness']['Mean'] = np.full((N, numLab), np.nan)
    AllResults['Thickness']['Median'] = np.full((N, numLab), np.nan)

    AllResults['Volume'] = {}
    AllResults['Volume']['Total'] = np.full((N, numLab), np.nan)
    if(ICV):
        AllResults['Volume']['ICV'] = np.full((N, numLab), np.nan)

    for i in range(N):
        print(i+1)

        currID = ID[i]
        currDates = allThickVals[allThickVals['id'] == currID]['date']
        uniqDates = np.unique(currDates) # Date is to specify a image in a sequence of MRI images --> Need a way to specify this part. Currently the first date is used.

        print(str(currID))

        if (not ICV):
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

    savemat(outMatSaveFile, AllResults)
    return AllResults