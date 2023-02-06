%dorsAttnLabels = contains(LabelLUT.Label_Name2,'DorsAttn')';
%ventAttnLabels = contains(LabelLUT.Label_Name2,'VentAttn')'; 
%limbicAttnLabels = contains(LabelLUT.Label_Name2,'Limbic')'; 


%currNetwork=dorsAttnLabels;
%saveName='DorsalAttn';

%currNetwork=ventAttnLabels;
%saveName='VentralAttn';

%currNetwork=limbicAttnLabels;
saveName='AO';



regNames = AtlasToPathLUT.PathSpreadSheetNames;
SN= length(regNames);
regNames_L=cell(1,SN);
regNames_R=cell(1,SN);

    for j = 1:SN
        regNames_L{j} = [regNames{j} '-L'];
        regNames_R{j} = [regNames{j} '-R'];
    end
regNames=[regNames_L regNames_R];
    
matchedThicknessDataCat= [subMatchedPathDataGM(:,:,1) subMatchedPathDataGM(:,:,2)];

currCoM=[pathCoM(:,:,1); pathCoM(:,:,2)];

N = SN*2;
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
    
covMatCtrl = nan(N,N);
covMatTau = nan(N,N);
covMatTDP = nan(N,N);

cmpCovTau_gt_Ctrl = nan(N,N);
cmpCovTDP_gt_Ctrl = nan(N,N);
cmpCovTau_gt_TDP = nan(N,N);
cmpCovTDP_gt_Tau = nan(N,N);


for i = 1:N
    for j = 1:N
      
        if(i~=j)
        NCtrl=sum(~isnan(thickCtrl(:,i)) & ~isnan(thickCtrl(:,j)));
        NTau=sum(~isnan(thickTau(:,i)) & ~isnan(thickTau(:,j)));
        NTDP=sum(~isnan(thickTDP(:,i)) & ~isnan(thickTDP(:,j)));

        if(NCtrl > 3)
            covMatCtrl(i,j)=corr(thickCtrl(:,i),thickCtrl(:,j),'Rows','complete');
        end
        if(NTau > 3)
            covMatTau(i,j)=corr(thickTau(:,i),thickTau(:,j),'Rows','complete');
        end
        if(NTDP > 3)
            covMatTDP(i,j)=corr(thickTDP(:,i),thickTDP(:,j),'Rows','complete');
        end
        
        if(NCtrl > 3 && NTau > 3 && NTDP > 3)
            cmpCovTau_gt_Ctrl(i,j) = TestCorrSig(covMatCtrl(i,j),covMatTau(i,j),NCtrl,NTau)<.05;
            cmpCovTDP_gt_Ctrl(i,j) = TestCorrSig(covMatCtrl(i,j),covMatTDP(i,j),NCtrl,NTDP)<.05;
            cmpCovTau_gt_TDP(i,j) = TestCorrSig(covMatTDP(i,j),covMatTau(i,j),NTDP,NTau)<.05;
            cmpCovTDP_gt_Tau(i,j) = TestCorrSig(covMatTau(i,j),covMatTDP(i,j),NTau,NTDP)<.05;
        end
        end
        
    end
end

Tau_Degree=sum(covMatTau,1,'omitnan');
TDP_Degree=sum(covMatTDP,1,'omitnan');

plotON=1;
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
imagesc(covMatCtrl)
title('controls')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_ctrl.png']),'-dpng','-r400');



H=figure(2)
imagesc(covMatTau)
title('Tau')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNames)
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_tau.png']),'-dpng','-r400');


H=figure(3)
imagesc(covMatTDP)
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
plotNetwork3(cmpCovTau_gt_Ctrl, schaefer400x7SurfaceMeshFromGII, currCoM,baseSaveDir,saveName,LabelNames,cRange,[],MakerVec,panelAll,pSelect,colorVec)


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
plotNetwork3(cmpCovTDP_gt_Ctrl, schaefer400x7SurfaceMeshFromGII, currCoM,baseSaveDir,saveName,LabelNames,cRange,[],MakerVec,panelAll,pSelect,colorVec)

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
plotNetwork3(cmpCovTau_gt_TDP, schaefer400x7SurfaceMeshFromGII, currCoM,baseSaveDir,saveName,LabelNames,cRange,[],MakerVec,panelAll,pSelect,colorVec)



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
plotNetwork3(cmpCovTDP_gt_Tau, schaefer400x7SurfaceMeshFromGII, currCoM,baseSaveDir,saveName,LabelNames,cRange,[],MakerVec,panelAll,pSelect,colorVec)

print(H2,fullfile(baseSaveDir,[ saveName '_SigPlots.png']),'-dpng','-r400');
end