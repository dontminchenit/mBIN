H2=figure(2)
clf

plotg1=[1 2 3 4];
panelAll = panel();
panelAll.pack(4,6);

panelAll.de.margin = 1;
panelAll.marginbottom = 0;
panelAll.margin = [1 1 1 1];

offset=0.00015;

Tau1TDP0=strcmp(PathInfo_T.Group,'tau');

% TauGMPath = log(GMPath(Tau1TDP0,:)+offset);
% TDPGMPath = log(GMPath(~Tau1TDP0,:)+offset);
% 
% TauWMPath = log(WMPath(Tau1TDP0,:)+offset);
% TDPWMPath = log(WMPath(~Tau1TDP0,:)+offset);

%rescale by all groups
[TauGMPath, TauWMPath] = rescalePathData(GMPath(Tau1TDP0,:),WMPath(Tau1TDP0,:));
[TDPGMPath, TDPWMPath] = rescalePathData(GMPath(~Tau1TDP0,:),WMPath(~Tau1TDP0,:));

%rescale by subgropus



%figure(5)
%imagesc(threshMatrix)

%plot Tau
dataGM=mean(TauGMPath,1,'omitnan')';
dataWM=mean(TauWMPath,1,'omitnan')';

dataComb=[dataGM dataWM];
dataComb=rescale(dataComb);
CRange=[0 1];

%GM
%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,dataComb(:,1),CRange,saveDir,saveName,0,panelAll,1,1,25)
%WM
%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,dataComb(:,2),CRange,saveDir,saveName,0,panelAll,2,1,25)


%plot TDP
dataGM=mean(TDPGMPath,1,'omitnan')';
dataWM=mean(TDPWMPath,1,'omitnan')';

dataComb=[dataGM dataWM];
dataComb=rescale(dataComb);
CRange=[0 1];

%GM
%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,dataComb(:,1),CRange,saveDir,saveName,0,panelAll,3,1,25)
%WM
%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,dataComb(:,2),CRange,saveDir,saveName,0,panelAll,4,1,25)


% print(H2,fullfile(saveDirBase,'pathPlot.tif'),'-dpng','-r400');




%find the mean FA for Group 1
HCmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{1});

%find the mean FA for Group 2
TAUmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{2});

%find the mean FA for Group 3
TDPmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{3});


thresh = 0;%;prctile(CoMDistMatrix(:),50);
threshMatrix = double(CoMDistMatrix > thresh);
threshMatrix(threshMatrix~=1)=nan;

[TAUWeakEdge, TauWeakEdgeReg, TAUIndvZscore, TAUMeanZscore, TAUIndvZscoreReg] = calcWeakEdges(TAUmtxs,HCmtxs,threshMatrix);
[TDPWeakEdge, TDPWeakEdgeReg, TDPIndvZscore, TDPMeanZscore, TDPIndvZscoreReg] = calcWeakEdges(TDPmtxs,HCmtxs,threshMatrix);


useZscore=0;

for q = 1:4
    
    switch q
        case 1
            CurrPath=TauGMPath;
            CurrMetric = TauWeakEdgeReg;
            if(useZscore)CurrMetric = TAUIndvZscoreReg;end
        case 2
            CurrPath=TauWMPath;
            CurrMetric = TauWeakEdgeReg;
            if(useZscore)CurrMetric = TAUIndvZscoreReg;end
        case 3
            CurrPath=TDPGMPath;
            CurrMetric = TDPWeakEdgeReg;
            if(useZscore)CurrMetric = TDPIndvZscoreReg;end
        case 4
            CurrPath=TDPWMPath;
            CurrMetric = TDPWeakEdgeReg;
            if(useZscore)CurrMetric = TDPIndvZscoreReg;end
    end
CorrPath = nan(1,LabelN);
CorrPval = nan(1,LabelN);
   SampleCount = nan(1,LabelN);
   
       AllEdge = [];
    AllPath = [];
for i = 1:LabelN
    
    edgeCount = CurrMetric(i,:)';
    
    %pathVal = TauGMPath(:,i);
    pathVal = CurrPath(:,i);
    
         AllEdge = [AllEdge ; edgeCount];
         AllPath = [AllPath ; pathVal];
    
       SampleCount(i) = sum(~isnan(edgeCount) & ~isnan(pathVal));
    
     if(SampleCount(i)>=5) 
        [rho,pval]=corr(pathVal,edgeCount, 'rows','complete');

        CorrPath(i)=rho;
        CorrPval(i)=pval;
    
    
        mdl = fitlm(edgeCount,pathVal);
    %    CorrPath(i)=mdl.Coefficients.Estimate(2);
    %    CorrPval(i)=mdl.Coefficients.pValue(2);
    end
end

data=CorrPath;
data(CorrPval>=.05)=nan;
%data(SampleCount<5)=nan;
CRange=[-1 1];

%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,q,2,25)


figure(10)
%AllPath(AllEdge<-5)=[];
%AllEdge(AllEdge<-5)=[];

%AllPath(AllEdge>5)=[];
%AllEdge(AllEdge>5)=[];


scatter(AllEdge,AllPath);
corr(AllEdge,AllPath,'rows','complete')

end

print(H2,fullfile(saveDirBase,'edgeCountRegress_sample05_pval.tif'),'-dpng','-r400');
%print(H2,fullfile(saveDirBase,'edgeCorr_sample05_pval5.tif'),'-dpng','-r400');


