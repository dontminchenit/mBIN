function plotGraphs(plotMtx,gmT_Curr,baseSaveDir,saveName,GMWMOffset,MarkerSize,Tau1TDP2)

% G = graph(adjMtx,'upper');
% 
% h=plot(G,'interpreter','none');
% %labelnode(h,1,gmT_Curr.Properties.VariableNames')
% for i = 1:23
% labelnode(h,i,gmT_Curr.Properties.VariableNames{i})
% end

% weight=nonzeros(adjMtx);
% h.EdgeCData=weight; 


laussaneBase = 'D:/Min/Dropbox (Personal)/Research/Projects/PPA/Code/NetworkViewer';
%laussaneBase = '/Users/min/Dropbox (Personal)/Research/Projects/PPA/Code/NetworkViewer';


LaussLabelLUT = readtable('./copath_labels_lausanne250_20191104_V2.lut.csv');
loadLaussMesh = load(fullfile(laussaneBase,'lausanne250_mesh.mat'));
MeshLauss = loadLaussMesh.mesh;

meshlabelsLauss = MeshLauss.attributes.attribute_array;
%find center of mass, to show the location of nodes
LaussCoM=zeros(height(LaussLabelLUT),3);
for i = 1:height(LaussLabelLUT)
    currLabel = LaussLabelLUT{i,1};
    idx = find(meshlabelsLauss==currLabel);
    LaussCoM(i,:)=mean(MeshLauss.points(idx,:),1);
end



PathLabelName=gmT_Curr.Properties.VariableNames;

[Path2LaussIdx, PathCoM, OffsetV1, OffsetV2] = PathLabels2Laussane(PathLabelName,LaussLabelLUT,LaussCoM);


    SN=length(gmT_Curr.Properties.VariableNames);
    PathLabelName = cell(1,SN);
    for j = 1:SN
    nameSplit = strsplit(gmT_Curr.Properties.VariableNames{j},'_');
    PathLabelName{j} = nameSplit{:,3};
    end

    %We use this to only show the Hippocampus Region for the revision
    hipIndex = strcmp(PathLabelName,'HIP')
    plotMtx(~hipIndex,~hipIndex)=0;

figure(3)
%MarkerSize=ones(1,size(plotMtx,1));
plotNetwork(plotMtx, MeshLauss, PathCoM,baseSaveDir,saveName,PathLabelName,[0 1],1,MarkerSize, OffsetV1, OffsetV2,GMWMOffset,Tau1TDP2)

%plotNetwork2(plotMtx, MeshLauss, PathCoM,saveDir,[saveName 'load'],PathLabelName,[-10 0],(k+4),nanmedian(gmT_Curr{:,:},1))
