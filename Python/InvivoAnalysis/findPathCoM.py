import numpy as np
import pandas as pd

# findPathCoM(pathLUT, AtlasToPathLUT, NetworkDataGeneral['NetworkDataGeneral'][0,0]['Schaefer400x7']['CoM'][0, 0])
def findPathCoM(pathLUT, AtlasToPathLUT, AtlasCoM):
    """Function to find the Center of Mass for Pathology Regions in Atlas

    Args:
        pathLUT (Dataframe): Look up table matching Atlas region names to Atlas Labels(Index)
        AtlasToPathLUT (Dataframe): Look up table matching Pathology region names to Atlas region names
        AtlasCoM (obj): Table containing CoM for each Atlas Labeles(Index)

    Returns:
        pathCoM (ndarray): Center of Mass of Pathology Regions (Mean of multiple corresponding Atlas regions Center of Mass)
        pathToAtlasIndex (ndarray): ndarray of list of Atlas regions index corrresponding to Pathology regions (this could be multiple Atlas regions)
    """

    # Get number of Pathology regions to match to Atlas
    NL = len(AtlasToPathLUT)

    # numpy to save the Center of Mass of the Pathology regions
    pathCoM = np.zeros((NL, 3, 2))

    # List of matching Pathology region to atlas index
    pathToAtlasIndex = [[None, None] for _ in range(NL)]

    hemi = ['L', 'R']
    
    for l in range(NL):
        for r in range(2):
            # Get the matched Atlas region
            AOName = AtlasToPathLUT.Matched[l]

            # Get row indexes of where pathname match Atlas region and {L, R} hemisphere
            # (modified csv file column name Label.Hemi --> Label_Hemi)
            pathLabelIdx = np.where((pathLUT.PathName.str.match(AOName) & pathLUT.Label_Hemi.str.match(hemi[r])))[0]

            if not (pathLabelIdx.size == 0): # If there is a match
                # Get Label_ID400x7 for the matched rows (modified csv file column name Label.ID400x7 --> Label_ID400x7)
                # Need to offset by 1 (index difference between matlab and python)
                atlasIdx = (pathLUT.Label_ID400x7.loc[pathLabelIdx] - 1).to_numpy()
                
                # Get CoM of Pathology Regions (Mean of multiple corresponding Atlas regions CoMs)
                pathCoM[l,:,r] = np.nanmean(AtlasCoM[atlasIdx,:], axis=0) # Approximation. Works if the regions are the similar volume. 

                # Get list of Atlas regions index corrresponding to Pathology regions (this could be multiple Atlas regions)
                pathToAtlasIndex[l][r] = atlasIdx
    
    return pathCoM, pathToAtlasIndex
