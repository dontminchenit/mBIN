import pandas as pd
import os
import numpy as np

# Location of the data folder
dataDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Data'

# pathT: ex-vivo histopathology Data (Quantification) / %AO
new_pathT = pd.read_excel(os.path.join(dataDir, 'NewFTDData', 'FTLD_Library_3-31-22_update.xlsx'), dtype={'INDDID': str})

# Region list
new_pathT_Region_list = new_pathT['Region'].tolist()
# Make it unique
new_pathT_Region_list = np.unique(np.array(new_pathT_Region_list))
new_pathT_Region_list = list(map(str, new_pathT_Region_list)) # cast to str

# Mapping from Anatomical regions to Atlas regions (This could be one to many, ex-aCING  multiple ACC regions)
pathLUT = pd.read_csv(os.path.join(dataDir,'schaefer_path_20210719_20220328.csv'))

# Look Up Table Region list
pathLUT_Region_list = pathLUT['PathName'].tolist()
# Make it unique
pathLUT_Region_list = np.unique(np.array(pathLUT_Region_list))
pathLUT_Region_list = list(map(str, pathLUT_Region_list)) # cast to str

print(new_pathT_Region_list)
#print(len(new_pathT_Region_list))
print(pathLUT_Region_list)
#print(len(pathLUT_Region_list))

# Check if what regions of new_Path is in Current LUT Region
test1 = np.sum(np.in1d(new_pathT_Region_list, pathLUT_Region_list)) # Number of Regions of new_Path in Path LUT Region / out of 22 unique new_Path LUT

print('Number of Regions of new_Path in Path LUT Region: ' + str(test1) + ' out of 22 unique Regions of New_Path')
print('The Overlapping Region is: ' + str(np.take(new_pathT_Region_list, np.argwhere(np.in1d(new_pathT_Region_list, pathLUT_Region_list) == True).flatten())))
print('The missing Region is: ' + str(np.take(new_pathT_Region_list, np.argwhere(np.in1d(new_pathT_Region_list, pathLUT_Region_list) == False).flatten())))

# Read Lookup table to match anatomical regions in the brain to the Atlas region
AtlasToPathLUT = pd.read_excel(os.path.join(dataDir,'LBDData','PathToAtlasLUT_10_7_2022.xlsx'))

AtlasToPathLUT_Path_list = AtlasToPathLUT['PathSpreadSheetNames'].tolist()
AtlasToPathLUT_Atlas_list = AtlasToPathLUT['Matched'].tolist()

# Check if what new_path Regions are in path regions of AtlasToPathLUT
test2 = np.sum(np.in1d(new_pathT_Region_list, AtlasToPathLUT_Path_list)) # Number of Regions of new_Path in path regions of AtlasToPathLUT / out of 22 unique new_Path LUT

print('Number of Regions of new_Path in path regions of AtlasToPathLUT: ' + str(test2) + ' out of 22 unique Regions of New_Path')
print('The Overlapping Region is: ' + str(np.take(new_pathT_Region_list, np.argwhere(np.in1d(new_pathT_Region_list, AtlasToPathLUT_Path_list) == True).flatten())))
print('The missing Region is: ' + str(np.take(new_pathT_Region_list, np.argwhere(np.in1d(new_pathT_Region_list, AtlasToPathLUT_Path_list) == False).flatten())))
