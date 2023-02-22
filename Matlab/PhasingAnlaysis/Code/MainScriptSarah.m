winBase = fullfile('D:','Min','Dropbox (Personal)','Research','Projects','PhasingAnalysis');
macBase = fullfile('/Users','Min','Dropbox (Personal)','Research','Projects','PhasingAnalysis');
currBase = MacBase;



inInfo = fullfile(currBase,'Data','dti_subs_5subRemoved.csv');
NetDataDir= fullfile(currBase,'Data','pathcsvs20210125','csvsSent');
outTableSaveFile = fullfile(currBase,'Data','Sarah_AllGraphs_7_15_2020.xlsx');
outMatSaveFile = fullfile(currBase,'Results','allMeasures_AllGraphs_Sarah_7_29_2020.mat');

[AllResults, outTable] = LoadNetworkDataSarah(inInfo,NetDataDir,outTableSaveFile,outMatSaveFile);

%PlotMeanAndComparisons
%PlotMeanAndComparisonsHubs