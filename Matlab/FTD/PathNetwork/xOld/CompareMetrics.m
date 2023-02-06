
saveDir = './PaperExpCore_3_23_2021/NetMetrics/';
mkdir(saveDir);

TauGMMtxM = calcNetworkMetrics(TauGMMtx);
TDPGMMtxM = calcNetworkMetrics(TDPGMMtx);
H=figure(1)
clf
PlotMetricComparisons(TauGMMtxM,TDPGMMtxM)
title('GM network')
print(H,fullfile(saveDir,'GM.png'),'-dpng','-r400')

TauWMMtxM = calcNetworkMetrics(TauWMMtx);
TDPWMMtxM = calcNetworkMetrics(TDPWMMtx);

H=figure(2)
clf
PlotMetricComparisons(TauWMMtxM,TDPWMMtxM)
title('WM network')
print(H,fullfile(saveDir,'WM.png'),'-dpng','-r400')

TauAllMtxM = calcNetworkMetrics(TauAllMtx);
TDPAllMtxM = calcNetworkMetrics(TDPAllMtx);

H=figure(3)
clf
PlotMetricComparisons(TauAllMtxM,TDPAllMtxM)
title('GM+WM network')
print(H,fullfile(saveDir,'GM+WM.png'),'-dpng','-r400')

% 
% 
% %%%%
% winBase = fullfile('D:','Min','Dropbox (Personal)','Research','Projects','PPA');
% macBase = fullfile('/Users','Min','Dropbox (Personal)','Research','Projects','PPA');
% currBase = winBase;
% load(fullfile(currBase,'Data','ProcessedResults','./Lausanne250SurfaceMeshFromGII.mat'))
% 
% 
% laussaneBase = '../../PPA/Code/NetworkViewer';
% LaussLabelLUT = readtable('./copath_labels_lausanne250_20191104_V2.lut.csv');
% loadLaussMesh = load(fullfile(laussaneBase,'lausanne250_mesh.mat'));
% MeshLauss = loadLaussMesh.mesh;
% 
% meshlabelsLauss = MeshLauss.attributes.attribute_array;
% %find center of mass, to show the location of nodes
% LaussCoM=zeros(height(LaussLabelLUT),3);
% for i = 1:height(LaussLabelLUT)
%     currLabel = LaussLabelLUT{i,1};
%     idx = find(meshlabelsLauss==currLabel);
%     LaussCoM(i,:)=mean(MeshLauss.points(idx,:),1);
% end
% 
% 
% 
% PathLabelName=gmwmT_AllTau.Properties.VariableNames;
% 
% [Path2LaussIdx, PathCoM] = PathLabels2Laussane(PathLabelName,LaussLabelLUT,LaussCoM);
% 
% dataTau = size(Path2LaussIdx,2);
% dataTDP = size(Path2LaussIdx,2);
% 
% for i = 1:size(Path2LaussIdx,1)
%     
%     curr = Path2LaussIdx(i,:);
%     idx = curr(~isnan(curr));
%     dataTau(idx) = TauAllMtxM.LocalEff(i);
%     dataTDP(idx) = TDPAllMtxM.LocalEff(i);
% end
% 
% 
% H2=figure(5)
% clf
% panelAll = panel();
% panelAll.pack(2,6);
% 
% panelAll.de.margin = 1;
% panelAll.marginbottom = 0;
% panelAll.margin = [1 1 1 1];
% 
% 
% CRange = [-1 1]

%plotBrain3(Lausanne250SurfaceMeshFromGII,dataTau,CRange,saveDir,saveName,0,panelAll,1)

%plotBrain3(Lausanne250SurfaceMeshFromGII,dataTDP,CRange,saveDir,saveName,0,panelAll,2)