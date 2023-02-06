dorsAttnLabels = contains(LabelLUT.Label_Name2,'DorsAttn')';
ventAttnLabels = contains(LabelLUT.Label_Name2,'VentAttn')'; 
limbicAttnLabels = contains(LabelLUT.Label_Name2,'Limbic')'; 


%currNetwork=dorsAttnLabels;
%saveName='DorsalAttn';

%currNetwork=ventAttnLabels;
%saveName='VentralAttn';

currNetwork=limbicAttnLabels;
saveName='Limbic';

currNetwork=ones(size(limbicAttnLabels));
saveName='all';


regNamesRaw = LabelLUT.Label_Name2(currNetwork);
SN= length(regNamesRaw);
regNames=cell(1,SN);
    
    for j = 1:SN
        nameSplit = strrep(regNamesRaw{j},'_','-');
        regNames{j} = nameSplit;
    end


currCoM=CoM(currNetwork,:);

N = sum(currNetwork);
%thickCtrl = AllResults.Thickness.Mean(CohortsIndexs==1,currNetwork);
pathTau = pathDataGM(CohortsIndexs==2,currNetwork);
pathTDP = pathDataGM(CohortsIndexs==3,currNetwork);



%covMatCtrl = nan(N,N);
covMatTau = nan(N,N);
covMatTDP = nan(N,N);

%cmpCovTau_gt_Ctrl = nan(N,N);
%cmpCovTDP_gt_Ctrl = nan(N,N);
cmpCovTau_gt_TDP = nan(N,N);
cmpCovTDP_gt_Tau = nan(N,N);


for i = 1:N
    for j = 1:N
      
        %covMatCtrl(i,j)=corr(thickCtrl(:,i),thickCtrl(:,j),'Rows','complete');


        if(i~=j)
        %NCtrl=sum(~isnan(thickCtrl(:,i)) & ~isnan(thickCtrl(:,j)));
        NTau=sum(~isnan(pathTau(:,i)) & ~isnan(pathTau(:,j)))
        NTDP=sum(~isnan(pathTDP(:,i)) & ~isnan(pathTDP(:,j)))

        %cmpCovTau_gt_Ctrl(i,j) = TestCorrSig(covMatCtrl(i,j),covMatTau(i,j),NCtrl,NTau)<.05;
        %cmpCovTDP_gt_Ctrl(i,j) = TestCorrSig(covMatCtrl(i,j),covMatTDP(i,j),NCtrl,NTDP)<.05;
        if(NTau > 3 && NTDP >3)
            
        covMatTau(i,j)=corr(pathTau(:,i),pathTau(:,j),'Rows','complete');
        covMatTDP(i,j)=corr(pathTDP(:,i),pathTDP(:,j),'Rows','complete');
            
        cmpCovTau_gt_TDP(i,j) = TestCorrSig(covMatTDP(i,j),covMatTau(i,j),NTDP,NTau)<.05;
        cmpCovTDP_gt_Tau(i,j) = TestCorrSig(covMatTau(i,j),covMatTDP(i,j),NTau,NTDP)<.05;
        end
        end
        
    end
end


H2=figure(10)
clf
panelAll = panel();
panelAll.pack(2,4);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

% H=figure(1)
% imagesc(covMatCtrl)
% title('controls')
% caxis([0 1])
% colormap('jet')
% set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
% set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
% set(gca,'TickLabelInterpreter','none','FontSize',5)
% print(H,fullfile(baseSaveDir,[ saveName '_ctrl.png']),'-dpng','-r400');


H=figure(2)
imagesc(covMatTau)
title('Tau')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_tau.png']),'-dpng','-r400');


H=figure(3)
imagesc(covMatTDP)
title('TDP')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TDP.png']),'-dpng','-r400');


% H=figure(4)
% imagesc(cmpCovTau_gt_Ctrl)
% title('Tau > Ctrl')
% caxis([0 1])
% colormap('jet')
% set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
% set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
% set(gca,'TickLabelInterpreter','none','FontSize',5)
% print(H,fullfile(baseSaveDir,[ saveName '_TAUGTCTRL.png']),'-dpng','-r400');
% 
 cRange = [0 1];
 MakerVec=ones(1,sum(currNetwork));
 colorVec=ones(1,sum(currNetwork));
LabelNames=regNamesRaw;
 pSelect=1;
% plotNetwork3(cmpCovTau_gt_Ctrl, schaefer400x7SurfaceMeshFromGII, currCoM,baseSaveDir,saveName,LabelNames,cRange,[],MakerVec,panelAll,pSelect,colorVec)
% 
% 
% H=figure(5)
% imagesc(cmpCovTDP_gt_Ctrl)
% title('TDP > Ctrl')
% caxis([0 1])
% colormap('jet')
% set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
% set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
% set(gca,'TickLabelInterpreter','none','FontSize',5)
% print(H,fullfile(baseSaveDir,[ saveName '_TDPGTCTRL.png']),'-dpng','-r400');
% pSelect=2;
% plotNetwork3(cmpCovTDP_gt_Ctrl, schaefer400x7SurfaceMeshFromGII, currCoM,baseSaveDir,saveName,LabelNames,cRange,[],MakerVec,panelAll,pSelect,colorVec)

H=figure(6)
imagesc(cmpCovTau_gt_TDP)
title('Tau > TDP')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TAUGTTDP.png']),'-dpng','-r400');
pSelect=3;
plotNetwork3(cmpCovTau_gt_TDP, schaefer400x7SurfaceMeshFromGII, currCoM,baseSaveDir,saveName,LabelNames,cRange,[],MakerVec,panelAll,pSelect,colorVec)



H=figure(7)
imagesc(cmpCovTDP_gt_Tau)
title('TDP > Tau')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_1TDPGTTAU.png']),'-dpng','-r400');
pSelect=4;
plotNetwork3(cmpCovTDP_gt_Tau, schaefer400x7SurfaceMeshFromGII, currCoM,baseSaveDir,saveName,LabelNames,cRange,[],MakerVec,panelAll,pSelect,colorVec)

print(H2,fullfile(baseSaveDir,[ saveName '_SigPlots.png']),'-dpng','-r400');