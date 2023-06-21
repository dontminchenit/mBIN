import numpy as np
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

def thickDistCmp(HCData, TAUData, TDPData, numLab):

    datalist = [HCData, TAUData, TDPData]
    dataNamelist = ['HC', 'TAU', 'TDP']
    g1_index = [0, 0, 1]
    g2_index = [1, 2, 2]

    # Get the Thickness value of HC, TAU, and TDP / thickHC, thickTAU, thickTDP
    for i in range(3): # HC vs TAU, HC vs TDP, TAU vs TDP
        g1 = datalist[g1_index[i]]
        g2 = datalist[g2_index[i]]

        testName = dataNamelist[g1_index[i]] + ' VS ' + dataNamelist[g2_index[i]]

        ttest_pVal_all = np.full(numLab, np.nan) # shape (1, numLab)
        ttest_tstat_all = np.full(numLab, np.nan) # shape (1, numLab)

        # do t-test on all labels
        for j in range(numLab):
            t_stat, pval = stats.ttest_ind(g1, g2)

            ttest_pVal_all[j] = pval
            ttest_tstat_all[j] = t_stat
        
        # correct for multiple comparison
        notNaN = ~np.isnan(TTestsAll_pVal) # Get get Boolean value -> 1 for pvalues that are not NaN.
        reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(TTestsAll_pVal[notNaN], alpha = 0.001, method = 'fdr_bh') # Alpha = 0.001 / Method: Benjamini/Hochberg (non-negative) --> Equivalent to 'fdr' in Matlab
        # reject - true for hypothesis that can be rejected for given alpha
        # pvals_corrected - p-values corrected for multiple tests
        # alphacSidak - corrected alpha for Sidak method
        # alphacBonf - corrected alpha for Bonferroni method







i = 4  # thickness
j = 6  # mean

H2 = plt.figure(2)
plt.clf()
panelAll = panel()
panelAll.pack(3, 6)

panelAll.de.margin = 1
panelAll.marginbottom = 1
panelAll.marginleft = -10
panelAll.margin = [0, 0, 0, 0]

measure = AllResults[graphNames[i]][measureNames[j]]
if j == 1:  # if betweenness, then we need to normalize
    measure = measure / ((numLab - 1) * (numLab - 2))

# we have the measure for each subgroup
GroupMeasures = [None] * 3
GroupMeasures[0] = measure[HCIND, :]
GroupMeasures[1] = measure[svIND, :]
GroupMeasures[2] = measure[naIND, :]

saveDirSub = os.path.join(saveDirBase, graphNames[i], measureNames[j])
CRange = [0, np.percentile(measure, 95)]

# We make these comparisons between groups
testG1All = [2, 3, 2, 3]
testG2All = [1, 1, 3, 2]
# plotg1 = [1, 3, 5, -1]
# plotg2 = [2, 4, 6, 7]
plotg1 = [0, 0, 3, 0]
plotg2 = [1, 2, 0, 0]
for g in range(3):
    G1 = testG1All[g]
    G2 = testG2All[g]
    testName = GroupNames[G1] + 'vs' + GroupNames[G2]
    TTestsAll_pVal = np.full(numLab, np.nan)
    TTestsAll_tstat = np.full(numLab, np.nan)
    # do t-test on all labels
    for t in range(numLab):
        h1, p, ci, stats = ttest2(
            GroupMeasures[G1][:, t], GroupMeasures[G2][:, t], 'left'
        )
        TTestsAll_pVal[t] = p
        TTestsAll_tstat[t] = stats.tstat

    # correct for multiple comparison
    notNaN = ~np.isnan(TTestsAll_pVal)
    padj_fdr, alpha_fdr = multicmp(TTestsAll_pVal[notNaN], 'fdr', 0.001)

    allPadj = [padj_fdr]
    allPadjName = [
        'FDR_Benjamini-Hochberg',
        'FWE_Hochberg-Bonferroni',
        'FWE_Holm-Bonferroni',
    ]

    for p in range(1):  # range(3):
        data = TTestsAll_pVal
        data[~notNaN] = 1
        data[notNaN] = allPadj[p]
        CRange = [0, 1]
        saveName = testName
        saveDir = os.path.join(saveDirSub, testName, allPadjName[p])
        # plotBrain2(labelmesh, data, LabelLUT, CRange, saveDir, saveName, 0.001, panelAll, plotg2[g])
        # plotBrain3(Lausanne250SurfaceMeshFromGII, data, CRange, saveDir,
