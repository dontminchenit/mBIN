import numpy as np
import os
import sys
import pickle
from functools import reduce
import matplotlib.pyplot as plt
from ipywidgets import interactive, FloatSlider, IntSlider
from scipy import stats
import scipy
from matplotlib.animation import FuncAnimation
import imageio
from scipy.sparse import lil_matrix, triu, coo_matrix


gitRepoDir = "/Users/hyroh/Desktop/FTD_Research/mBIN"
dataDir='/Users/hyroh/Desktop/FTD_Research/Data'
baseDir = os.path.join(gitRepoDir, 'Matlab')
loadData_hf = os.path.join(gitRepoDir, 'FTD_JupyterNotebook/HelperFunctions/Load_Dataset_HelperFunctions')
analysis_hf = os.path.join(gitRepoDir, 'FTD_JupyterNotebook/HelperFunctions/Analysis_HelperFuntions')
path_dataDir= os.path.join(gitRepoDir, 'FTD_JupyterNotebook/Data/Pathology_Data')
path_CalData = os.path.join(gitRepoDir, 'FTD_JupyterNotebook/Pathology_Analysis/Calculated_Data')
thick_dataDir= os.path.join(gitRepoDir, 'FTD_JupyterNotebook/Data/Thickness_Data')
thick_CalData = os.path.join(gitRepoDir, 'FTD_JupyterNotebook/Thickness_Analysis/Calculated_Data')
thickAtPath_dataDir= os.path.join(gitRepoDir, 'FTD_JupyterNotebook/Data/ThicknessAtPath_Data')
thickAtPath_CalData = os.path.join(gitRepoDir, 'FTD_JupyterNotebook/Thickness_At_Pathology_Analysis/Calculated_Data')

sys.path.insert(0, analysis_hf)
import nodalStrength as ns
import normalize as norm


# Function to create adjacency matrix
def create_adjacency_matrix(vertices, faces):
    # Initialize adj Matrix
    num_vertices = len(vertices)
    adjacency_matrix = np.zeros((num_vertices, num_vertices), dtype=int)
    
    for face in faces: # for each face
        for i in range(3): # For each pair of vertice indices (total of 3 pairs)
            v1 = face[i]
            v2 = face[(i + 1) % 3]
            # Update the Adjacency matrix
            adjacency_matrix[v1, v2] = 1
            adjacency_matrix[v2, v1] = 1

    return adjacency_matrix

# Function to create adjacency matrix
def create_adjacency_matrix_sparse(vertices, faces):
    # Initialize adj Matrix
    num_vertices = len(vertices)
    adjacency_matrix = lil_matrix((num_vertices, num_vertices), dtype=int)
    
    for face in faces: # for each face
        for i in range(3): # For each pair of vertice indices (total of 3 pairs)
            v1 = face[i]
            v2 = face[(i + 1) % 3]
            # Update the Adjacency matrix
            adjacency_matrix[v1, v2] = 1
            adjacency_matrix[v2, v1] = 1

    return adjacency_matrix

def sliceAndOrderAdjMax(adjMat, vert_index_toThickCoM, thickCoM_ordered_list):
    """
    args:
        adjMat (ndarray): 2D Adjacent Matrix of Vertices (327684 x 327684)
        
        vert_index_toThickCoM (list): contains list of vertices index that match the 400 thickness regions in order. 
        Therefore, vert_index_toThickCoM[0] is the vertice index that matches to Thickness CoM index 0.
        
        thickCoM_ordered_list (list): List containing the sliced vertices indices in order (of Thickness CoM Indices)
    """
    # slice the adjcaency matrix (unordered)
    adjMat = adjMat[vert_index_toThickCoM, :][:, vert_index_toThickCoM]
    
    # Order the sliced adjcaency matrix
    adjMat = adjMat[thickCoM_ordered_list][:, thickCoM_ordered_list]
    
    return adjMat

def strongConnection_adjMat(pathToAtlasIndex_list, adj_mat):
    
    strong_connection_at_path_list = []

    for i in range(len(pathToAtlasIndex_list)):
        # For each Pathology Node
        connection_list = []
        for j in pathToAtlasIndex_list[i]: # Get indices of equivalent Thickness region (out of 400)
            # Get indices of thickness regions that is connected in adjacency matrix
            connections = np.nonzero(adj_mat[j])[0].tolist()

            # Save each connections into a list
            connection_list.append(connections)

        # Get connection that unique 
        strong_connection_list = np.unique(np.array(connection_list))

        # Get the indicies of the pathology regions (that is strongly connected)
        strong_connection_at_path = []
        for index in range(len(pathToAtlasIndex_list)):
            for el in strong_connection_list:
                if el in pathToAtlasIndex_list[index]:
                    strong_connection_at_path.append(index)

        # Get only the unique strongly connected pathology regions             
        strong_connection_at_path = np.unique(np.array(strong_connection_at_path))

        strong_connection_at_path_list.append(strong_connection_at_path)

    return strong_connection_at_path_list


# loads the preconstructed Atlas data
NetworkDataGeneral = scipy.io.loadmat(os.path.join(baseDir, 'NetworkAnalysisGeneral', 'FTDGeneralData_20221114.mat'))

with open(os.path.join(path_dataDir, 'pathToAtlasIndex.pkl'), 'rb') as f:
    pathToAtlasIndex = pickle.load(f)
f.close()

# cov_volAtPath_w_dict_Drop
with open(os.path.join(thickAtPath_CalData, 'cov_volAtPath_w_dict_Drop.pkl'), 'rb') as f:
    cov_volAtPath_w_dict_Drop = pickle.load(f)
f.close()

# path_TAU_Drop
with open(os.path.join(path_dataDir, 'path_TAU_Drop.pkl'), 'rb') as f:
    path_TAU_Drop = pickle.load(f)
f.close()

# path_TDP_Drop
with open(os.path.join(path_dataDir, 'path_TDP_Drop.pkl'), 'rb') as f:
    path_TDP_Drop = pickle.load(f)
f.close()

# pathCoM
with open(os.path.join(path_dataDir, 'pathCoM.pkl'), 'rb') as f:
    pathCoM = pickle.load(f)
f.close()
pathCoM = np.vstack((pathCoM[:, :, 0], pathCoM[:, :, 1]))

# CoM_TAU_Drop
with open(os.path.join(path_dataDir, 'CoM_TAU_Drop.pkl'), 'rb') as f:
    CoM_TAU_Drop = pickle.load(f)
f.close()

# CoM_TDP_Drop
with open(os.path.join(path_dataDir, 'CoM_TDP_Drop.pkl'), 'rb') as f:
    CoM_TDP_Drop = pickle.load(f)
f.close()

# LabelNames
with open(os.path.join(path_CalData, 'LabelNames.pkl'), 'rb') as f:
    LabelNames = pickle.load(f)
f.close()

# pathNames_TAU_Drop
with open(os.path.join(path_CalData, 'pathNames_TAU_Drop.pkl'), 'rb') as f:
    pathNames_TAU_Drop = pickle.load(f)
f.close()

# pathNames_TDP_Drop
with open(os.path.join(path_CalData, 'pathNames_TDP_Drop.pkl'), 'rb') as f:
    pathNames_TDP_Drop = pickle.load(f)
f.close()

# TAU_missing_index
with open(os.path.join(path_CalData, 'TAU_missing_index.pkl'), 'rb') as f:
    TAU_missing_index = pickle.load(f)
f.close()

# TDP_missing_index
with open(os.path.join(path_CalData, 'TDP_missing_index.pkl'), 'rb') as f:
    TDP_missing_index = pickle.load(f)
f.close()

# Min/Max Range of Normalizing
t_min = -1
t_max = 1

# normalizing TAU EXCLUDING NaN!
path_TAU_Drop_Norm = norm.normalize2d(path_TAU_Drop, t_min, t_max)

# normalizing TDP EXCLUDING NaN!
path_TDP_Drop_Norm = norm.normalize2d(path_TDP_Drop, t_min, t_max)

pathOrig = np.concatenate((np.nanmean(path_TAU_Drop_Norm, axis=0), 
                           np.nanmean(path_TDP_Drop_Norm, axis=0)))

ymin = np.min(pathOrig)
ymax = np.max(pathOrig)

pathToAtlasIndex_list = []

for i in range(len(pathToAtlasIndex)):
    curr_l = pathToAtlasIndex[i][0]
    pathToAtlasIndex_list.append(curr_l)

for i in range(len(pathToAtlasIndex)):
    curr_r = pathToAtlasIndex[i][1]
    pathToAtlasIndex_list.append(curr_r)

pathToAtlasIndex_list_TAU = pathToAtlasIndex_list.copy()

for i in TAU_missing_index[::-1]:
    del pathToAtlasIndex_list_TAU[i]

pathToAtlasIndex_list_TDP = pathToAtlasIndex_list.copy()

for i in TDP_missing_index[::-1]:
    del pathToAtlasIndex_list_TDP[i]



mesh_vertices = NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['giiSurface_Both'][0,0]['vertices'][0, 0]
mesh_faces = NetworkDataGeneral['NetworkDataGeneral'][0, 0]['Schaefer400x7']['GII'][0, 0]['giiSurface_Both'][0,0]['faces'][0, 0] - 1

thick_CoM = NetworkDataGeneral['NetworkDataGeneral'][0,0]['Schaefer400x7']['CoM'][0, 0]

# Load the vert_index_toThickCoM
with open(os.path.join(thickAtPath_CalData, 'vert_index_toThickCoM.pkl'), 'rb') as f:
    vert_index_toThickCoM = pickle.load(f)
f.close()

# Create a mapping from thickness region to indices
region_to_indices = {region: None for region in np.unique(vert_index_toThickCoM)}
for i, region in enumerate(vert_index_toThickCoM):
    region_to_indices[region] = i
    
# Get keys in ascending order
sorted_keys = sorted(region_to_indices.keys())

# Retrieve values based on sorted keys
thickCoM_ordered_list = [region_to_indices[key] for key in sorted_keys]





# Create adjacency matrix - Sparse
vert_adj_matrix_sparse = create_adjacency_matrix_sparse(mesh_vertices, mesh_faces)
# vert_adj_matrix = vert_adj_matrix_sparse.toarray()

# # 1 Connection
# vert_adj_matrix_toThick = sliceAndOrderAdjMax(vert_adj_matrix, vert_index_toThickCoM, thickCoM_ordered_list)
# print(vert_adj_matrix_toThick)


# 2 Connection
vert_adj_matrix_sparse_2 = (triu(vert_adj_matrix_sparse) @ vert_adj_matrix_sparse).astype(bool).astype(int)
vert_adj_matrix_2 = vert_adj_matrix_sparse_2.toarray()

vert_adj_matrix_toThick_2 = sliceAndOrderAdjMax(vert_adj_matrix_2, vert_index_toThickCoM, thickCoM_ordered_list)
print(vert_adj_matrix_toThick_2)

print(strongConnection_adjMat(pathToAtlasIndex_list_TAU, vert_adj_matrix_toThick_2))