function [AllResults, outTable] = LoadNetworkData(inInfo,NetDataDir,outTableSaveFile,outMatSaveFile)

infoTable = readtable(inInfo);
outTable = infoTable;
ID = infoTable.INDDID;
N = length(ID);
AllResults = struct;
numLab=448;
%numLab=118;

% for GraphEdgeType = 1:2
%     
%     %we have 3 ways to construct a graph
%     switch GraphEdgeType
%         case 1
%             fileInSuffix = '_sc.csv';
%             graphName = 'SCGraph';
%         case 2
%             fileInSuffix = '_MeanTractMedianFA.csv';
%             graphName = 'FAGraph';
%         case 3
%             fileInSuffix = '_MeanTractLength.csv';
%             graphName = 'TractLenGraph';
%     end
%     
%     %initialize all the measures
%     AllResults.(graphName).('WeightMtx') = cell(N,1);
%     AllResults.(graphName).('WeightMtxNormal') = cell(N,1);
%     AllResults.(graphName).('LengthMtxNormal') = cell(N,1);
%     AllResults.(graphName).('BetweenCen') = nan(N,numLab);
%     AllResults.(graphName).('ClusterCoeff') = nan(N,numLab);
%     AllResults.(graphName).('LocalEff') = nan(N,numLab);
%     AllResults.(graphName).('EigenCen') = nan(N,numLab);
%     AllResults.(graphName).('Degree') = nan(N,numLab);
%     AllResults.(graphName).('MeanStr') = nan(N,1);
%     AllResults.(graphName).('CharPathLen') = nan(N,1);
%     AllResults.(graphName).('GlobalEff') = nan(N,1);
%     AllResults.(graphName).('GlobalEffNorm') = nan(N,1);
%     AllResults.(graphName).('NodalEcc') = nan(N,numLab);
%     AllResults.(graphName).('NetRadius') = nan(N,1);
%     AllResults.(graphName).('NetDiameter') = nan(N,1);
%     
%     %initialize table for global measures
%     outTable.([graphName 'MeanStr']) = nan(N,1);
%     outTable.([graphName 'CharPathLen']) = nan(N,1);
%     outTable.([graphName 'GlobalEff']) = nan(N,1);
%     outTable.([graphName 'GlobalEffNorm']) = nan(N,1);
%     
%     %Start Calculations
%     for i = 1:N
%         i
%         %current ID and Data
%         currID = ID(i);
%         curDate = infoTable.MRIDatec(i);
%         
%         %load the adj matrix
%         currMtxFile = dir(fullfile(NetDataDir,[num2str(currID) '_' num2str(curDate) fileInSuffix]));
%         currAdjMtx = table2array(readtable(fullfile(NetDataDir,currMtxFile(1).name)));
%         %set weight and length matrices
%         W=currAdjMtx;
%         W_normalized = weight_conversion(W, 'normalize');
%         L = weight_conversion(W, 'lengths');
%         
%         %calc the local measures
%         AllResults.(graphName).('WeightMtx'){i} = W;
%         AllResults.(graphName).('WeightMtxNormal'){i} = W_normalized;
%         AllResults.(graphName).('LengthMtxNormal'){i} = L;
%         AllResults.(graphName).('BetweenCen')(i,:) = betweenness_wei(L);
%         AllResults.(graphName).('ClusterCoeff')(i,:) = clustering_coef_wu(W_normalized);
%         AllResults.(graphName).('LocalEff')(i,:) = efficiency_wei(W_normalized, 2);
%         AllResults.(graphName).('EigenCen')(i,:) = eigenvector_centrality_und(W);
%         degree = sum(W,2);
%         AllResults.(graphName).('Degree')(i,:) = degree;
%         
%         %calc global measures
%         GlobEffNorm = efficiency_wei(W_normalized);
%         [D, B]=distance_wei(L);
%         [lambda,efficiency,ecc,radius,diameter] = charpath(D);
%         
%         AllResults.(graphName).('MeanStr')(i) = mean(degree);
%         AllResults.(graphName).('CharPathLen')(i) = lambda;
%         AllResults.(graphName).('GlobalEff')(i) = efficiency;
%         AllResults.(graphName).('GlobalEffNorm')(i) = GlobEffNorm;
%         AllResults.(graphName).('NodalEcc')(i,:) = ecc;
%         AllResults.(graphName).('NetRadius')(i) = radius;
%         AllResults.(graphName).('NetDiameter')(i) = diameter;
%         
%         %write global values to table
%         outTable.([graphName 'MeanStr'])(i) = mean(degree);
%         outTable.([graphName 'CharPathLen'])(i) = lambda;
%         outTable.([graphName 'GlobalEff'])(i) = efficiency;
%         outTable.([graphName 'GlobalEffNorm'])(i) = GlobEffNorm;
%     end
% end

 AllResults.Thickness.Mean = nan(N,numLab);
 AllResults.Thickness.Median = nan(N,numLab);
 AllResults.Volume.Total = nan(N,numLab);
 

%Thickness and Volume Values
for i = 1:N
    i
    currID = ID(i);
    curDate = infoTable.MRIDatec(i);
    
    currMeasureFile = dir(fullfile(NetDataDir,[num2str(currID) '_' num2str(curDate) '_lausanne.csv']));
    
    currThickVals = readtable(fullfile(NetDataDir,currMeasureFile(1).name));
    currThickVals = sortrows(currThickVals,'label');
    AllResults.Thickness.Median(i,:) = currThickVals{strcmpi(currThickVals.measure,'thickness') & strcmpi(currThickVals.metric,'median'),'value'};
    AllResults.Thickness.Mean(i,:) = currThickVals{strcmpi(currThickVals.measure,'thickness') & strcmpi(currThickVals.metric,'mean'),'value'};
    AllResults.Volume.Total(i,:) = currThickVals{strcmpi(currThickVals.system,'lausanne250') & strcmpi(currThickVals.measure,'volume') & strcmpi(currThickVals.metric,'numeric'),'value'};
end

save(outMatSaveFile,'AllResults');
writetable(outTable,outTableSaveFile);