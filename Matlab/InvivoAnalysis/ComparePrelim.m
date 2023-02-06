

%convert to z-scores


dorsAttnLabels = contains(LabelLUT.Label_Name2,'DorsAttn')';
ventAttnLabels = contains(LabelLUT.Label_Name2,'VentAttn')'; 


thickCtrl = AllResults.Thickness.Mean(CohortsIndexs==1,:);
thickTau = AllResults.Thickness.Mean(CohortsIndexs==2,:);
thickTDP = AllResults.Thickness.Mean(CohortsIndexs==3,:);

thickTau_z = nan(size(thickTau));
thickTDP_z = nan(size(thickTDP));

cmptstatData=nan(numLab,1);
pathCmppvalData=nan(numLab,1);
Tau_pvalData=nan(numLab,1);
tautstatData=nan(numLab,1);
TDP_pvalData=nan(numLab,1);
tdptstatData=nan(numLab,1);
for i = 1:numLab

    currCtrlMean = mean(thickCtrl(:,i),'omitnan');
    currCtrlStd = std(thickCtrl(:,i),'omitnan');
    
    thickTau_z(:,i)=(thickTau(:,i)-currCtrlMean)/currCtrlStd;
    thickTDP_z(:,i)=(thickTDP(:,i)-currCtrlMean)/currCtrlStd;

    [H,P,CI,STATS]=ttest2(thickCtrl(:,i),thickTau(:,i),'Tail','right');
    tautstatData(i) = STATS.tstat;
    Tau_pvalData(i)=P;

    [H,P,CI,STATS]=ttest2(thickCtrl(:,i),thickTDP(:,i),'Tail','right');
    tdptstatData(i) = STATS.tstat;
    TDP_pvalData(i)=P;

    
    %[H,P,CI,STATS]=ttest2(thickTau_z(:,i),thickTDP_z(:,i));
    [H,P,CI,STATS]=ttest2(thickTau(:,i),thickTDP(:,i));
    cmptstatData(i) = STATS.tstat;
    pathCmppvalData(i)=P;
end

saveDirSub='D:\Min\Dropbox (Personal)\Research\Projects\FTD\InvivoPath\PrelimResults'

H2=figure(2)
        clf
        panelAll = panel();
        panelAll.pack(4,6);

        panelAll.de.margin = 1;
        panelAll.marginbottom = 1;
        panelAll.marginleft = -10;
        panelAll.margin = [0 0 0 0];


testName='prelim';
saveName = testName;
saveDir = fullfile(saveDirSub,testName);
currTData = cmptstatData;
currPData = pathCmppvalData;

% currTData = tautstatData;
% currPData = Tau_pvalData;
% 
% currTData = tdptstatData;
% currPData = TDP_pvalData;

data=currTData .* dorsAttnLabels';
CRange=[-4 4];
plotBrain3Points(CoM,schaefer400x7SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,1,2,10)
 
 data=currPData .* dorsAttnLabels';
 data(dorsAttnLabels==0)=nan;
 CRange=[0 1];
 plotBrain3Points(CoM,schaefer400x7SurfaceMeshFromGII,data,CRange,saveDir,saveName,0.01,panelAll,2,1,10)
 
 data=currTData .* ventAttnLabels';
 CRange=[-4 4];
 plotBrain3Points(CoM,schaefer400x7SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,3,2,10)
% 
data=currPData .* ventAttnLabels';
data(ventAttnLabels==0)=nan;
CRange=[0 1];
plotBrain3Points(CoM,schaefer400x7SurfaceMeshFromGII,data,CRange,saveDir,saveName,0.01,panelAll,4,1,10)

print(H2,fullfile(saveDir,'cmpNetwork.tif'),'-dpng','-r400');

%         H2=figure(2)
%         panelAll = panel();
%         panelAll.pack(1,6);
% 
%         panelAll.de.margin = 1;
%         panelAll.marginbottom = 1;
%         panelAll.marginleft = -10;
%         panelAll.margin = [0 0 0 0];
% 
% 
% data = TTestsAll_pVal
% data(~notNaN) = 1;
% data(notNaN) = allPadj{p};
% CRange=[0 1];
% saveName = testName;
% saveDir = fullfile(saveDirSub,testName,allPadjName{p});
%plotBrain2(labelmesh,data,LabelLUT,CRange,saveDir,saveName,0.001,panelAll,plotg2(g))
%plotBrain3(Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0.001,panelAll,plotg2(g))
%plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0.001,panelAll,1,2,10)






