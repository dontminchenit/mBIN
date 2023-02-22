%LoadDataSarah

%
H=figure(1)
clf
panelAll = panel();
panelAll.pack(2,6);

panelAll.de.margin = 1;
panelAll.marginbottom = 0;
panelAll.margin = [1 1 1 1];

%plot TauStaging

%CRange=[-5 5]
%plotBrain3(Lausanne250SurfaceMeshFromGII,T_tauStagingLaus,CRange,saveDir,saveName,0,panelAll,1)
data=T_tauStagingLaus;
data(data==0) = 6;
CRange = [0 6]
plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,data,0,panelAll,1,7,10)

%plotNetwork3(zeroMtx, Lausanne250SurfaceMeshFromGII, CoM,saveDirBase,'phases',[],[0 1],1, keepIndx+.001,panelAll,s,ColorVec)


%plot TDPStaging
%plotBrain3(Lausanne250SurfaceMeshFromGII,T_tdpStagingLaus,CRange,saveDir,saveName,0,panelAll,2)

data=T_tdpStagingLaus;
data(data==0) = 6;
plotBrain3Points(CoM,Lausanne250SurfaceMeshFromGII,data,0,panelAll,2,7,10)


print(H,fullfile(saveDirBase,'fig1.tif'),'-dpng','-r400'); 



%boxplot of each stage with each measure and tau vs tdp
%Fig7bSarah

%Fig7cSarah