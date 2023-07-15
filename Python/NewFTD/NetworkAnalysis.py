import pickle
import networkx as nx
from networkx.generators.random_graphs import erdos_renyi_graph
import time
from netrd.distance import Hamming

def hammingHelperfunction(G1, G2, G3, G4):
    dist_obj = Hamming()

    dis_1_3 = dist_obj.dist(G1, G3)
    dis_1_2 = dist_obj.dist(G1, G2)
    dis_1_4 = dist_obj.dist(G1, G4)

    print(f"Distance between Pathology TAU > TDP vs Thickness At Path TAU > TDP: {dis_1_3}")
    print(f"Distance between Pathology TAU > TDP vs Pathology TDP > TAU: {dis_1_2}")
    print(f"Distance between Pathology TAU > TDP vs Thickness At Path TDP > TAU FD: {dis_1_4}")

outputDir='/Users/hyung/Research23_Network_Analysis/NetworkAnalysisData-selected/Output_Python_NewFTD'

######### Pathology vs Thickness At Path  #########
# Pathology TAU > TDP  # 1
path_TAU_gt_TDP = pickle.load(open(outputDir + '/Graph_Objects/FTD_TAU_GT_TDP_GM.pickle', 'rb'))
# Pathology TDP > TAU  # 2
path_TDP_gt_TAU = pickle.load(open(outputDir + '/Graph_Objects/FTD_TDP_GT_TAU_GM.pickle', 'rb'))
# Thickness At Path TAU > TDP  # 3
thick_at_path_TAU_gt_TDP = pickle.load(open(outputDir + '/Graph_Objects/cmpCovTAU_gt_TDP_Path_pthresh_05.pickle', 'rb'))
# Thickness At Path TDP > TAU # 4
thick_at_path_TDP_gt_TAU = pickle.load(open(outputDir + '/Graph_Objects/cmpCovTDP_gt_TAU_Path_pthresh_05.pickle', 'rb'))

hammingHelperfunction(path_TAU_gt_TDP, path_TDP_gt_TAU, thick_at_path_TAU_gt_TDP, thick_at_path_TDP_gt_TAU)

######### Pathology FD = 40 vs Thickness At Path FD = 40 #########
# Pathology TAU > TDP FD = 40 edges # 1
path_TAU_gt_TDP_FD_40 = pickle.load(open(outputDir + '/Graph_Objects/covTAU_gt_TDP_FD_40_GM.pickle', 'rb'))
# Pathology TDP > TAU FD = 40 edges # 2
path_TDP_gt_TAU_FD_40 = pickle.load(open(outputDir + '/Graph_Objects/covTDP_gt_TAU_FD_40_GM.pickle', 'rb'))
# Thickness At Path TAU > TDP FD = 40 edges  # 3
thick_at_path_TAU_gt_TDP_FD_40 = pickle.load(open(outputDir + '/Graph_Objects/cmpCovTAU_gt_TDP_FD_40_Path_pthresh_05.pickle', 'rb'))
# Thickness At Path TDP > TAU FD = 40 edges # 4
thick_at_path_TDP_gt_TAU_FD_40 = pickle.load(open(outputDir + '/Graph_Objects/cmpCovTDP_gt_TAU_FD_40_Path_pthresh_05.pickle', 'rb'))

hammingHelperfunction(path_TAU_gt_TDP_FD_40, path_TDP_gt_TAU_FD_40, thick_at_path_TAU_gt_TDP_FD_40, thick_at_path_TDP_gt_TAU_FD_40)

######### Pathology vs Thickness At Path W-Score  #########
# Pathology TAU > TDP  # 1
path_TAU_gt_TDP = pickle.load(open(outputDir + '/Graph_Objects/FTD_TAU_GT_TDP_GM.pickle', 'rb'))
# Pathology TDP > TAU  # 2
path_TDP_gt_TAU = pickle.load(open(outputDir + '/Graph_Objects/FTD_TDP_GT_TAU_GM.pickle', 'rb'))
# Thickness At Path TAU > TDP  / W Score # 3
thick_at_path_TAU_gt_TDP_WScore = pickle.load(open(outputDir + '/Graph_Objects/cmpCovTAU_gt_TDP_Path_W_Path_pthresh_05_WScore.pickle', 'rb'))
# Thickness At Path TDP > TAU  / W Score # 4
thick_at_path_TDP_gt_TAU_WScore = pickle.load(open(outputDir + '/Graph_Objects/cmpCovTDP_gt_TAU_Path_W_Path_pthresh_05_WScore.pickle', 'rb'))

hammingHelperfunction(path_TAU_gt_TDP, path_TDP_gt_TAU, thick_at_path_TAU_gt_TDP_WScore, thick_at_path_TDP_gt_TAU_WScore)

######### Pathology FD = 40 vs Thickness At Path W-Score FD = 40 #########
# Pathology TAU > TDP FD = 40 edges # 1
path_TAU_gt_TDP_FD_40 = pickle.load(open(outputDir + '/Graph_Objects/covTAU_gt_TDP_FD_40_GM.pickle', 'rb'))
# Pathology TDP > TAU FD = 40 edges # 2
path_TDP_gt_TAU_FD_40 = pickle.load(open(outputDir + '/Graph_Objects/covTDP_gt_TAU_FD_40_GM.pickle', 'rb'))
# Thickness At Path TAU > TDP FD = 40 edges / W Score # 3
thick_at_path_TAU_gt_TDP_FD_40_WScore = pickle.load(open(outputDir + '/Graph_Objects/cmpCovTAU_gt_TDP_FD_Path_W_40_Path_pthresh_05_WScore.pickle', 'rb'))
# Thickness At Path TDP > TAU FD = 40 edges / W Score # 4
thick_at_path_TDP_gt_TAU_FD_40_WScore = pickle.load(open(outputDir + '/Graph_Objects/cmpCovTDP_gt_TAU_FD_Path_W_40_Path_pthresh_05_WScore.pickle', 'rb'))

hammingHelperfunction(path_TAU_gt_TDP_FD_40, path_TDP_gt_TAU_FD_40, thick_at_path_TAU_gt_TDP_FD_40_WScore, thick_at_path_TDP_gt_TAU_FD_40_WScore)

