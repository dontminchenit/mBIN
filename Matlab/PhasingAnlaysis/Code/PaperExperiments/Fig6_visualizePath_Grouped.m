%must run the following
%Fig8_relateconnectivityToPath.m
%Fig8_matchSubToPathInfo

H2=figure(2)
clf

plotg1=[1 2 3 4];
panelAll = panel();
panelAll.pack(4,6);

panelAll.de.margin = 1;
panelAll.marginbottom = 0;
panelAll.margin = [1 1 1 1];


Tau1TDP0=strcmp(PathInfo_T.Group,'tau');

PSP1= strcmp(PathInfo_T.subgroup(Tau1TDP0),'psp');
PID1= strcmp(PathInfo_T.subgroup(Tau1TDP0),'pid');
CBD1= strcmp(PathInfo_T.subgroup(Tau1TDP0),'cbd');
%Tau1TDP0=strcmp(PathInfo_T.subgroup,'psp');

%Separate data into Tau/TDP and WM/GM
TauGMPathLabel = GMPathLabel(Tau1TDP0,:);
TauWMPathLabel = GMPathLabel(Tau1TDP0,:);

TDPGMPathLabel = WMPathLabel(~Tau1TDP0,:);
TDPWMPathLabel = WMPathLabel(~Tau1TDP0,:);

TauDiseaseDur=SarahDiseaseDur(Tau1TDP0);
TDPDiseaseDur=SarahDiseaseDur(~Tau1TDP0);

%TauGMPath = log(GMPath(Tau1TDP0,:)+offset);
%TDPGMPath = log(GMPath(~Tau1TDP0,:)+offset);

%TauWMPath = log(WMPath(Tau1TDP0,:)+offset);
%TDPWMPath = log(WMPath(~Tau1TDP0,:)+offset);

%rescale PathData and make log
[TauGMPath, TauWMPath] = rescalePathData(GMPath(Tau1TDP0,:),WMPath(Tau1TDP0,:));
[TDPGMPath, TDPWMPath] = rescalePathData(GMPath(~Tau1TDP0,:),WMPath(~Tau1TDP0,:));

%rescale PathData subgropus
% [TauGMPathPSP, TauWMPathPSP] = rescalePathData(GMPath(PSP1,:),WMPath(PSP1,:));
% [TauGMPathPID, TauWMPathPID] = rescalePathData(GMPath(PID1,:),WMPath(PID1,:));
% [TauGMPathCBD, TauWMPathCBD] = rescalePathData(GMPath(CBD1,:),WMPath(CBD1,:));
% 
% TauGMPath(PSP1,:)=TauGMPathPSP;
% TauGMPath(PID1,:)=TauGMPathPID;
% TauGMPath(CBD1,:)=TauGMPathCBD;
% 
% TauWMPath(PSP1,:)=TauWMPathPSP;
% TauWMPath(PID1,:)=TauWMPathPID;
% TauWMPath(CBD1,:)=TauWMPathCBD;



%PathDataIn_T=rescalePathData(PathDataIn_T);


%find the center location for each path region.
TauGMCenterVal=findPathCenter(TauGMPath, TauGMPathLabel);
TauWMCenterVal=findPathCenter(TauWMPath, TauWMPathLabel);

TDPGMCenterVal=findPathCenter(TDPGMPath, TDPGMPathLabel);
TDPWMCenterVal=findPathCenter(TDPWMPath, TDPWMPathLabel);

%plot Tau
dataGM=mean(TauGMCenterVal,1,'omitnan')';
dataWM=mean(TauWMCenterVal,1,'omitnan')';

dataComb=[dataGM dataWM];
dataComb=rescale(dataComb);
CRange=[0 1];

%GM
%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,dataComb(:,1),CRange,saveDir,saveName,0,panelAll,1,1,25)
%WM
%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,dataComb(:,2),CRange,saveDir,saveName,0,panelAll,2,1,25)


%plot TDP
dataGM=mean(TDPGMCenterVal,1,'omitnan')';
dataWM=mean(TDPWMCenterVal,1,'omitnan')';

dataComb=[dataGM dataWM];
dataComb=rescale(dataComb);
CRange=[0 1];

%GM
%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,dataComb(:,1),CRange,saveDir,saveName,0,panelAll,3,1,25)
%WM
%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,dataComb(:,2),CRange,saveDir,saveName,0,panelAll,4,1,25)


%print(H2,fullfile(saveDirBase,'pathPlotGrouped.tif'),'-dpng','-r400');




%find the mean FA for Group 1
HCmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{1});

%find the mean FA for Group 2
TAUmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{2});

%find the mean FA for Group 3
TDPmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{3});

prctVal=0;

thresh = prctile(CoMDistMatrix(:),prctVal);

threshMatrix = double(CoMDistMatrix > thresh);
threshMatrix(threshMatrix~=1)=nan;

%figure(5)
%imagesc(threshMatrix)

[TAUWeakEdge, TauWeakEdgeReg, TAUIndvZscore, TAUMeanZscore, TAUIndvZscoreReg] = calcWeakEdges(TAUmtxs,HCmtxs,threshMatrix);
[TDPWeakEdge, TDPWeakEdgeReg, TDPIndvZscore, TDPMeanZscore, TDPIndvZscoreReg] = calcWeakEdges(TDPmtxs,HCmtxs,threshMatrix);


TauWeakEdgeReg_CenterGM=findPathCenterCount(TauWeakEdgeReg', TauGMPathLabel);
TauWeakEdgeReg_CenterWM=findPathCenterCount(TauWeakEdgeReg', TauWMPathLabel);

TDPWeakEdgeReg_CenterGM=findPathCenterCount(TDPWeakEdgeReg', TDPGMPathLabel);
TDPWeakEdgeReg_CenterWM=findPathCenterCount(TDPWeakEdgeReg', TDPWMPathLabel);

TauZscoreReg_CenterGM=findPathCenter(TAUIndvZscoreReg', TauGMPathLabel);
TauZscoreReg_CenterWM=findPathCenter(TAUIndvZscoreReg', TauWMPathLabel);

TDPZscoreReg_CenterGM=findPathCenter(TDPIndvZscoreReg', TDPGMPathLabel);
TDPZscoreReg_CenterWM=findPathCenter(TDPIndvZscoreReg', TDPWMPathLabel);


useZscore=1;


for q = 1:4
    
    switch q
        case 1
            CurrPath=TauGMPath;
            DiseaseDur = TauDiseaseDur;
            if(useZscore)
                CurrMetric = TauZscoreReg_CenterGM;
              
                

                edgename = 'Z-score';
            else
                CurrMetric = TauWeakEdgeReg_CenterGM;
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'Tau GM';
            
        case 2
            CurrPath=TauWMPath;
            DiseaseDur = TauDiseaseDur;
            if(useZscore)
                CurrMetric = TauZscoreReg_CenterWM;
                edgename = 'Z-score';
            else
                CurrMetric = TauWeakEdgeReg_CenterWM;
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'Tau WM';
            
            
        case 3
            CurrPath=TDPGMPath;
            DiseaseDur = TDPDiseaseDur;
            if(useZscore)
                CurrMetric = TDPZscoreReg_CenterGM;
                edgename = 'Z-score';
            else
                CurrMetric = TDPWeakEdgeReg_CenterGM;
                DiseaseDur = TDPDiseaseDur;
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'TDP43 GM';
            
        case 4
            CurrPath=TDPWMPath;
            DiseaseDur = TDPDiseaseDur;
            if(useZscore)
                CurrMetric = TDPZscoreReg_CenterWM;
                edgename = 'Z-score';
            else
                CurrMetric = TDPWeakEdgeReg_CenterWM;
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'TDP43 WM';
        case 5
            CurrPath=TauGMPath(PSP1,:);
            
            if(useZscore)
                CurrMetric = TauZscoreReg_CenterGM(PSP1,:);
                edgename = 'Z-score';
            else
                CurrMetric = TauWeakEdgeReg_CenterGM(PSP1,:);
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'Tau GM PSP';
            
        case 6
            CurrPath=TauWMPath(PSP1,:);
            if(useZscore)
                CurrMetric = TauZscoreReg_CenterWM(PSP1,:);
                edgename = 'Z-score';
            else
                CurrMetric = TauWeakEdgeReg_CenterWM(PSP1,:);
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'Tau WM PSP';
        case 7
            CurrPath=TauGMPath(PID1,:);
            
            if(useZscore)
                CurrMetric = TauZscoreReg_CenterGM(PID1,:);
                edgename = 'Z-score';
            else
                CurrMetric = TauWeakEdgeReg_CenterGM(PID1,:);
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'Tau GM PID';
            
        case 8
            CurrPath=TauWMPath(PID1,:);
            if(useZscore)
                CurrMetric = TauZscoreReg_CenterWM(PID1,:);
                edgename = 'Z-score';
            else
                CurrMetric = TauWeakEdgeReg_CenterWM(PID1,:);
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'Tau WM PID';
        case 9
            CurrPath=TauGMPath(CBD1,:);
            
            if(useZscore)
                CurrMetric = TauZscoreReg_CenterGM(CBD1,:);
                edgename = 'Z-score';
            else
                CurrMetric = TauWeakEdgeReg_CenterGM(CBD1,:);
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'Tau GM CBD';
            
        case 10
            CurrPath=TauWMPath(CBD1,:);
            if(useZscore)
                CurrMetric = TauZscoreReg_CenterWM(CBD1,:);
                edgename = 'Z-score';
            else
                CurrMetric = TauWeakEdgeReg_CenterWM(CBD1,:);
                edgename = 'Avg # Weak Edges';
            end
            pathname = 'Tau WM CBD';
            
    end
    CorrPath = nan(1,LabelN);
    CorrPval = nan(1,LabelN);
    RegressB = nan(1,LabelN);
    RegressPaval  = nan(1,LabelN);
    SampleCount = nan(1,LabelN);
    AllEdge = [];
    AllPath = [];
    AllSubID = [];
    AllRegID = [];
    for i = 1:LabelN
        
        edgeCount = CurrMetric(:,i);
        
        %pathVal = TauGMPath(:,i);
        pathVal = CurrPath(:,i);
        
        subID=[1:length(edgeCount)]';
        regID =i*ones(length(edgeCount),1);
        
         AllEdge = [AllEdge ; edgeCount];
         AllPath = [AllPath ; pathVal];
         AllSubID = [AllSubID ; subID];
         AllRegID = [AllRegID ; regID];
   
        
        SampleCount(i) = sum(~isnan(edgeCount) & ~isnan(pathVal));
        
        if(SampleCount(i)>=5)       
        [rho,pval]=corr(pathVal,edgeCount, 'rows','complete');
        
        CorrPath(i)=rho;
        CorrPval(i)=pval;
        
        %[b,bint,r,rint,stats] = regress(pathVal,edgeCount);
        
        mdl = fitlm(edgeCount,pathVal);
        mdl = fitlm(pathVal,edgeCount);
        CorrPath(i)=mdl.Coefficients.Estimate(2);
        CorrPval(i)=mdl.Coefficients.pValue(2);
        
        
        end
    end
    
    data=CorrPath;
    data(CorrPval>=.05)=nan;
    %data(SampleCount<5)=nan;
    CRange=[-1 1];
    
    H2=figure(2)
    %plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,q,2,25)
    
%     
%     H2=figure(3)
%     data=SampleCount;
%     plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,q,2,25)
%  

H10=figure(10)
clf
%AllPath(AllEdge<-6)=[];
%AllEdge(AllEdge<-6)=[];


%AllPath(AllEdge>6)=[];
%AllEdge(AllEdge>6)=[];


%scatter(AllEdge,AllPath,[],AllSubID);
N=size(CurrMetric,1);
AllColors=[]
durThresh=5
durThreshBot=0
CurrPath=[];
CurrEdge=[];
for currSub=1:N
%scatter(AllEdge(AllSubID==currSub),AllPath(AllSubID==currSub),[],AllSubID(AllSubID==currSub));
%colormap(distinguishable_colors(N))
length(AllEdge(AllSubID==currSub))
length(AllPath(AllSubID==currSub))
if(~isnan(DiseaseDur(currSub)) && DiseaseDur(currSub)> durThreshBot && DiseaseDur(currSub)<= durThresh && sum(~isnan(AllEdge(AllSubID==currSub)))>=2 && sum(~isnan(AllPath(AllSubID==currSub)))>=2)
%q=scatter(AllEdge(AllSubID==currSub),AllPath(AllSubID==currSub),[],DiseaseDur(currSub)*ones(length(AllSubID(AllSubID==currSub)),1));
CurrPath=[CurrPath;AllEdge(AllSubID==currSub)];
CurrEdge=[CurrEdge;AllPath(AllSubID==currSub)];
caxis([0 5])%max(DiseaseDur)])
colormap(jet)
%colorbar
hold on
%lsline
currColor = vals2colormap(DiseaseDur(currSub),[], [0 max(DiseaseDur)]);
AllColors=[AllColors;currColor];
end
end
%hlines = lsline;
%for k = 1:numel(hlines)
%    set(hlines(k), 'Color', AllColors(k,:))
%end

scatter(CurrEdge,CurrPath)
xlim([-9 -1])
ylim([-4 1.5])
lsline
[hc, pc]=corr(CurrEdge,CurrPath,'rows','complete')
mdl = fitlm(CurrEdge,CurrPath);
hr=mdl.Coefficients.Estimate(2)
pr=mdl.Coefficients.pValue(2)
ylabel('log(normalized %AO)')
xlabel('z-score')
set(gca,'FontSize',20)
title([pathname ' vs ' edgename '(b=' num2str(hr) ', p=' num2str(pr) ' || r=' num2str(hc) ', p=' num2str(pc) ')']);

[pathname ' vs ' edgename '(b=' num2str(hr) ', p=' num2str(pr) ' || r=' num2str(hc) ', p=' num2str(pc) ')']

%print(H10,fullfile(saveDirBase,[pathname '_vs_' edgename '_threshLT15.tif']),'-dpng','-r400');
print(H10,fullfile(saveDirBase,[pathname '_vs_' edgename '_durThresh' num2str(durThreshBot) 'to' num2str(durThresh) '_colored_wlines.tif']),'-dpng','-r400');
end

if(useZscore)
%    print(H2,fullfile(saveDirBase,['edgeReg_Grouped_zscore_sample5_pval5_prctGreater' num2str(prctVal) '.tif']),'-dpng','-r400');
    %print(H2,fullfile(saveDirBase,'edgeCorr_Grouped_zscore_sample1.tif'),'-dpng','-r400');
else
%    print(H2,fullfile(saveDirBase,'edgeCorr_Grouped_sample5.tif'),'-dpng','-r400');
end


