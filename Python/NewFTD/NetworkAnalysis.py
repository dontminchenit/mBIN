import pickle
import networkx as nx
from networkx.generators.random_graphs import erdos_renyi_graph
import time

outputDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Output_Python_NewFTD'

minv_TAU_GT_TDP = 0
minv_TDP_GT_TAU = 0
Random_minv = 0

# Two Random Graph for reference
n = 8
p = 0.5
gr1 = erdos_renyi_graph(n, p)
gr2 = erdos_renyi_graph(n, p)

for v in nx.optimize_graph_edit_distance(gr1, gr2):
    print(v)
    Random_minv = v

print("Random Reference N=5:")
print(Random_minv)


# load graph object from file

# Pathology TAU > TDP
Path_TAU_GT_TDP = pickle.load(open(outputDir + '/Graph_Objects/Pathology/FTD_TAU_GT_TDP_GM.pickle', 'rb'))

# Thickness AT Path TAU > TDP
Thick_TAU_GT_TDP = pickle.load(open(outputDir + '/Graph_Objects/Thickness_At_Path/cmpCovTAU_gt_TDP_Path_pthresh_05.pickle', 'rb'))

t0 = time.time()
for v in nx.optimize_graph_edit_distance(Path_TAU_GT_TDP, Thick_TAU_GT_TDP):
    print(v)
    minv_TAU_GT_TDP = v
print("TAU > TDP:")
print(minv_TAU_GT_TDP)
t1 = time.time()

total_n = t1-t0
print("time: " + str(total_n))

# Pathology TAU < TDP
Path_TDP_GT_TAU = pickle.load(open(outputDir + '/Graph_Objects/Pathology/FTD_TDP_GT_TAU_GM.pickle', 'rb'))

# Thickness AT Path TAU < TDP
Thick_TDP_GT_TAU = pickle.load(open(outputDir + '/Graph_Objects/Thickness_At_Path/cmpCovTDP_gt_TAU_Path_pthresh_05.pickle', 'rb'))

for v in nx.optimize_graph_edit_distance(Path_TDP_GT_TAU, Thick_TDP_GT_TAU):
    print(v)
    minv_TDP_GT_TAU = v

print("TDP > TAU:")    
print(minv_TDP_GT_TAU)

total_n = t1-t0
print("time: " + str(total_n))

# Two Random Graph for reference
n = 40
p = 0.5
gRandom1 = erdos_renyi_graph(n, p)
gRandom2 = erdos_renyi_graph(n, p)

for v in nx.optimize_graph_edit_distance(gRandom1, gRandom2):
    print(v)
    Random_minv = v

print("Random Reference:")
print(Random_minv)