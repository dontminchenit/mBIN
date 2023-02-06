H2=figure(10)
clf
panelAll = panel();
panelAll.pack(2,2);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

currCoM=[pathCoM(:,:,1); pathCoM(:,:,2)];
cRange = [0 1];
%MakerVec=ones(1,N);
PathMakerVecTau=mean(matchedPathDataCat(CohortsIndexs==2,:),1,'omitnan');
PathMakerVecTau(isnan(PathMakerVecTau))=0.1;
PathMakerVecTau=3*PathMakerVecTau/(max(PathMakerVecTau));
PathMakerVecTDP=mean(matchedPathDataCat(CohortsIndexs==3,:),1,'omitnan');
PathMakerVecTDP(isnan(PathMakerVecTDP))=0.1;
PathMakerVecTDP=3*PathMakerVecTDP/(max(PathMakerVecTDP));

ThickMakerVecTau=mean(Tau_W,1,'omitnan');
ThickMakerVecTau(isnan(ThickMakerVecTau))=0.1;
ThickMakerVecTau=.01+3*(ThickMakerVecTau-min(ThickMakerVecTau))/(max(ThickMakerVecTau)-min(ThickMakerVecTau));
ThickMakerVecTDP=mean(TDP_W,1,'omitnan');
ThickMakerVecTDP(isnan(ThickMakerVecTDP))=0.1;
ThickMakerVecTDP=.01+3*(ThickMakerVecTDP-min(ThickMakerVecTDP))/(max(ThickMakerVecTDP)-min(ThickMakerVecTDP));


colorVec=ones(1,N);
colorVec2=3*ones(1,N);
LabelNames=regNames;
%pSelect=1;

pathcmpCovTDP_gt_Tau=zeros(size(pathcmpCovTDP_gt_Tau))
pathcmpCovTau_gt_TDP=zeros(size(pathcmpCovTau_gt_TDP))

plotNetwork3_double(pathcmpCovTDP_gt_Tau, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,PathMakerVecTDP,ThickMakerVecTau,panelAll,1,colorVec,colorVec2,0)
plotNetwork3_double(pathcmpCovTau_gt_TDP, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,PathMakerVecTau,ThickMakerVecTDP,panelAll,2,colorVec,colorVec2,0)
print(H2,fullfile(baseSaveDir,[ saveName '_PathDiff.png']),'-dpng','-r400');