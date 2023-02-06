%dorsAttnLabels = contains(LabelLUT.Label_Name2,'DorsAttn')';
%ventAttnLabels = contains(LabelLUT.Label_Name2,'VentAttn')'; 
%limbicAttnLabels = contains(LabelLUT.Label_Name2,'Limbic')'; 


%currNetwork=dorsAttnLabels;
%saveName='DorsalAttn';

%currNetwork=ventAttnLabels;
%saveName='VentralAttn';

%currNetwork=limbicAttnLabels;
saveName='All';



regNames = AtlasToPathLUT.PathSpreadSheetNames;
SN= length(regNames);
N=SN;
regNames_L=cell(1,SN);
regNames_R=cell(1,SN);

    for j = 1:SN
        regNames_L{j} = [regNames{j} '-L'];
        regNames_R{j} = [regNames{j} '-R'];
    end
%regNames=[regNames_L regNames_R];


%make path Cov Matrix
matchedPathDataCat= mean(matchedPathDataGM,3,'omitnan');

pathTau = log(matchedPathDataCat(CohortsIndexs==2,:)+0.00015);
pathTDP = log(matchedPathDataCat(CohortsIndexs==3,:)+0.00015);

H=figure(2)
imagesc(pathTau)
title('Tau')
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
colorbar
%set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_tauPathLog.png']),'-dpng','-r400');


H=figure(3)
imagesc(pathTDP)
title('TDP')
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
colorbar
%set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TPDPathLog.png']),'-dpng','-r400');



N = SN;

pathcovMatTau = nan(N,N);
pathcovMatTDP = nan(N,N);

pathcmpCovTau_gt_TDP = nan(N,N);
pathcmpCovTDP_gt_Tau = nan(N,N);


for i = 1:N
    for j = 1:N
      
        if(i~=j)
        NTau=sum(~isnan(pathTau(:,i)) & ~isnan(pathTau(:,j)));
        NTDP=sum(~isnan(pathTDP(:,i)) & ~isnan(pathTDP(:,j)));

        if(NTau > 3)
            pathcovMatTau(i,j)=corr(pathTau(:,i),pathTau(:,j),'Rows','complete');
        end
        if(NTDP > 3)
            pathcovMatTDP(i,j)=corr(pathTDP(:,i),pathTDP(:,j),'Rows','complete');
        end
        
        
        if(NTau > 3 && NTDP > 3)
            pathcmpCovTau_gt_TDP(i,j) = TestCorrSig(pathcovMatTDP(i,j),pathcovMatTau(i,j),NTDP,NTau)<.05;
            pathcmpCovTDP_gt_Tau(i,j) = TestCorrSig(pathcovMatTau(i,j),pathcovMatTDP(i,j),NTau,NTDP)<.05;
        end
        
        end
        
    end
end

plotON=0
if(plotON)
H=figure(2)
imagesc(pathcovMatTau)
title('Tau')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_tauPathCov.png']),'-dpng','-r400');


H=figure(3)
imagesc(pathcovMatTDP)
title('TDP')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TPDPathCov.png']),'-dpng','-r400');


H=figure(6)
imagesc(pathcmpCovTau_gt_TDP)
title('Tau > TDP')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TAUGTTDP.png']),'-dpng','-r400');
pSelect=3;

H=figure(7)
imagesc(pathcmpCovTDP_gt_Tau)
title('TDP > Tau')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_1TDPGTTAU.png']),'-dpng','-r400');
pSelect=4;

H2=figure(10)
clf
panelAll = panel();
panelAll.pack(2,2);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

currCoM=pathCoM(:,:,1);% pathCoM(:,:,2)];
cRange = [0 1];
%MakerVec=ones(1,N);
MakerVecTau=mean(matchedPathDataCat(CohortsIndexs==2,:),1,'omitnan');
MakerVecTau(isnan(MakerVecTau))=0.1;
MakerVecTau=3*MakerVecTau/(max(MakerVecTau));
MakerVecTDP=mean(matchedPathDataCat(CohortsIndexs==3,:),1,'omitnan');
MakerVecTDP(isnan(MakerVecTDP))=0.1;
MakerVecTDP=3*MakerVecTDP/(max(MakerVecTDP));
colorVec=ones(1,N);
LabelNames=regNames;
%pSelect=1;

%plotNetwork3(pathcmpCovTDP_gt_Tau, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTDP,panelAll,1,colorVec,0)
%plotNetwork3(pathcmpCovTau_gt_TDP, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTau,panelAll,2,colorVec,0)
plotNetwork3(pathcovMatTDP, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTDP,panelAll,1,colorVec,0)
plotNetwork3(pathcovMatTau, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTau,panelAll,2,colorVec,0)
print(H2,fullfile(baseSaveDir,[ saveName '_PathDiff_all.png']),'-dpng','-r400');

end

%make thickness Cov Matrix
matchedThicknessDataCat= mean(matchedThicknessData,3,'omitnan');




N = SN;
thickCtrl = matchedThicknessDataCat(CohortsIndexs==1,:);
thickTau = matchedThicknessDataCat(CohortsIndexs==2,:);
thickTDP = matchedThicknessDataCat(CohortsIndexs==3,:);

Ctrl_Mean = mean(thickCtrl,1,'omitnan');
Ctrl_SD = std(thickCtrl,1,'omitnan');

Tau_W=nan(size(thickTau));
for i = 1:size(thickTau,1)
    Tau_W(i,:)= (thickTau(i,:) - Ctrl_Mean)./Ctrl_SD;
end

TDP_W=nan(size(thickTDP));
for i = 1:size(thickTDP,1)
    TDP_W(i,:)= (thickTDP(i,:) - Ctrl_Mean)./Ctrl_SD;
end
    
mean_Tau_w = mean(Tau_W,1,'omitnan');
mean_TDP_w = mean(TDP_W,1,'omitnan');
    
thickcovMatCtrl = nan(N,N);
thickcovMatTau = nan(N,N);
thickcovMatTDP = nan(N,N);

cmpCovTau_gt_Ctrl = nan(N,N);
cmpCovTDP_gt_Ctrl = nan(N,N);
cmpCovTau_gt_TDP = nan(N,N);
cmpCovTDP_gt_Tau = nan(N,N);


for i = 1:N
    for j = 1:N
      
        if(i~=j)
        %NCtrl=sum(~isnan(thickCtrl(:,i)) & ~isnan(thickCtrl(:,j)));
        NTau=sum(~isnan(Tau_W(:,i)) & ~isnan(Tau_W(:,j)));
        NTDP=sum(~isnan(TDP_W(:,i)) & ~isnan(TDP_W(:,j)));

        %if(NCtrl > 3)
        %    thickcovMatCtrl(i,j)=corr(thickCtrl(:,i),thickCtrl(:,j),'Rows','complete');
        %end
        if(NTau > 3)
            thickcovMatTau(i,j)=corr(Tau_W(:,i),Tau_W(:,j),'Rows','complete');
        end
        if(NTDP > 3)
            thickcovMatTDP(i,j)=corr(TDP_W(:,i),TDP_W(:,j),'Rows','complete');
        end
        
        if(NTau > 3 && NTDP > 3)
           % cmpCovTau_gt_Ctrl(i,j) = TestCorrSig(thickcovMatCtrl(i,j),thickcovMatTau(i,j),NCtrl,NTau)<.05;
           % cmpCovTDP_gt_Ctrl(i,j) = TestCorrSig(thickcovMatCtrl(i,j),thickcovMatTDP(i,j),NCtrl,NTDP)<.05;
            cmpCovTau_gt_TDP(i,j) = TestCorrSig(thickcovMatTDP(i,j),thickcovMatTau(i,j),NTDP,NTau)<.05;
            cmpCovTDP_gt_Tau(i,j) = TestCorrSig(thickcovMatTau(i,j),thickcovMatTDP(i,j),NTau,NTDP)<.05;
        end
        end
        
    end
end

plotON=1
if(plotON)

H=figure(2)
imagesc(Tau_W)
title('Tau')
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
colorbar
%set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_tauW.png']),'-dpng','-r400');


H=figure(3)
imagesc(TDP_W)
title('TDP')
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
colorbar
%set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TPDW.png']),'-dpng','-r400');

    
    
    
H=figure(2)
imagesc(thickcovMatTau)
title('Tau')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_tauWCov.png']),'-dpng','-r400');


H=figure(3)
imagesc(thickcovMatTDP)
title('TDP')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TPDWCov.png']),'-dpng','-r400');


H=figure(6)
imagesc(cmpCovTau_gt_TDP)
title('Tau > TDP')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TAUGTTDP_W.png']),'-dpng','-r400');
pSelect=3;

H=figure(7)
imagesc(cmpCovTDP_gt_Tau)
title('TDP > Tau')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TDPGTTAU_W.png']),'-dpng','-r400');
pSelect=4;

H2=figure(10)
clf
panelAll = panel();
panelAll.pack(2,2);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

currCoM=pathCoM(:,:,1);% pathCoM(:,:,2)];
cRange = [0 1];
MakerVec=ones(1,N);
colorVec=ones(1,N);
LabelNames=regNames;
%pSelect=1;

plotNetwork3(cmpCovTDP_gt_Tau, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTDP,panelAll,1,colorVec,0)
plotNetwork3(cmpCovTau_gt_TDP, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVecTau,panelAll,2,colorVec,0)
print(H2,fullfile(baseSaveDir,[ saveName '_PathDiff_W.png']),'-dpng','-r400');

end



Tau_Degree=sum(thickcovMatTau,1,'omitnan');
TDP_Degree=sum(thickcovMatTDP,1,'omitnan');

plotON=0;
if(plotON)
H2=figure(10)
clf
panelAll = panel();
panelAll.pack(2,4);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

H=figure(1)
imagesc(thickcovMatCtrl)
title('controls')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_ctrl.png']),'-dpng','-r400');



H=figure(2)
imagesc(thickcovMatTau)
title('Tau')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_tau.png']),'-dpng','-r400');


H=figure(3)
imagesc(thickcovMatTDP)
title('TDP')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TDP.png']),'-dpng','-r400');


H=figure(4)
imagesc(cmpCovTau_gt_Ctrl)
title('Tau > Ctrl')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TAUGTCTRL.png']),'-dpng','-r400');

cRange = [0 1];
MakerVec=ones(1,N);
colorVec=ones(1,N);
LabelNames=regNames;
pSelect=1;
plotNetwork3(cmpCovTau_gt_Ctrl, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVec,panelAll,pSelect,colorVec,0)


H=figure(5)
imagesc(cmpCovTDP_gt_Ctrl)
title('TDP > Ctrl')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TDPGTCTRL.png']),'-dpng','-r400');
pSelect=2;
plotNetwork3(cmpCovTDP_gt_Ctrl, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVec,panelAll,pSelect,colorVec,0)
%,UnderlayData)

H=figure(6)
imagesc(cmpCovTau_gt_TDP)
title('Tau > TDP')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TAUGTTDP.png']),'-dpng','-r400');
pSelect=3;
plotNetwork3(cmpCovTau_gt_TDP, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVec,panelAll,pSelect,colorVec,0)



H=figure(7)
imagesc(cmpCovTDP_gt_Tau)
title('TDP > Tau')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_1TDPGTTAU.png']),'-dpng','-r400');
pSelect=4;
plotNetwork3(cmpCovTDP_gt_Tau, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVec,panelAll,pSelect,colorVec,0)

print(H2,fullfile(baseSaveDir,[ saveName '_SigPlots.png']),'-dpng','-r400');
end