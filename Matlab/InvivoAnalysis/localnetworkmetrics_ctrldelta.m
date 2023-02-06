
thickcovMatTau_delta=cmpCovTau_gt_Ctrl;
thickcovMatTDP_delta=cmpCovTDP_gt_Ctrl;
thickcovMatTau_delta(isnan(thickcovMatTau_delta))=0;
thickcovMatTDP_delta(isnan(thickcovMatTDP_delta))=0;

H2=figure(10)
clf
panelAll = panel();
panelAll.pack(2,2);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

emptymtx=zeros(size(thickcovMatTDP_delta))


%MakerVecTau=mean(matchedPathDataCat(CohortsIndexs==2,:),1,'omitnan');
%MakerVecTDP=mean(matchedPathDataCat(CohortsIndexs==3,:),1,'omitnan');

MakerVecTau=mean(Tau_W,1,'omitnan');
MakerVecTDP=mean(TDP_W,1,'omitnan');

MakerVecTau=MakerVecTau-min(MakerVecTau);
MakerVecTDP=MakerVecTDP-min(MakerVecTDP);

MakerVecTau(isnan(MakerVecTau)| MakerVecTau == 0)=0.1;
MakerVecTau=3*MakerVecTau/(max(MakerVecTau));
MakerVecTDP(isnan(MakerVecTDP)| MakerVecTDP == 0)=0.1;
MakerVecTDP=3*MakerVecTDP/(max(MakerVecTDP));


%plotNetwork3(thickcovMatTDP_delta, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTDP,panelAll,1,colorVec,0)
%plotNetwork3(thickcovMatTau_delta, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTau,panelAll,2,colorVec,0)
%print(H2,fullfile(baseSaveDir,[ saveName 'ThicknessNetwork_delta.png']),'-dpng','-r400');


H2=figure(11)
clf
panelAll = panel();
panelAll.pack(2,2);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];


MakerVecTau=strengths_und(thickcovMatTau_delta);
MakerVecTDP=strengths_und(thickcovMatTDP_delta);

%MakerVecTau=betweenness_wei(thickcovMatTau_delta);
%MakerVecTDP=betweenness_wei(thickcovMatTDP_delta);

%MakerVecTau=abs(clustering_coef_wu(thickcovMatTau_delta));
%MakerVecTDP=abs(clustering_coef_wu(thickcovMatTDP_delta));

%MakerVecTau=eigenvector_centrality_und(thickcovMatTau_delta);
%MakerVecTDP=eigenvector_centrality_und(thickcovMatTDP_delta);

%MakerVecTau=participation_coef(thickcovMatTau_delta);
%MakerVecTDP=participation_coef(thickcovMatTDP_delta);


MakerVecTau(isnan(MakerVecTau)| MakerVecTau == 0)=0.1;
maxval=max(max(MakerVecTau),max(MakerVecTDP));
MakerVecTau=3*MakerVecTau/(maxval);
MakerVecTDP(isnan(MakerVecTDP)| MakerVecTDP == 0)=0.1;
MakerVecTDP=3*MakerVecTDP/(maxval);

plotNetwork3(emptymtx, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTDP,panelAll,1,colorVec,0)
plotNetwork3(emptymtx, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTau,panelAll,2,colorVec,0)

%print(H2,fullfile(baseSaveDir,[ saveName 'strength.png']),'-dpng','-r400');
print(H2,fullfile(baseSaveDir,[ saveName 'betweenness.png']),'-dpng','-r400');
%print(H2,fullfile(baseSaveDir,[ saveName 'clustering.png']),'-dpng','-r400');
%print(H2,fullfile(baseSaveDir,[ saveName 'eigenvector.png']),'-dpng','-r400');
%print(H2,fullfile(baseSaveDir,[ saveName 'participation.png']),'-dpng','-r400');




