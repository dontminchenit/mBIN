
%regNames=cell(1,SN);
    
%    for j = 1:SN
%        nameSplit = strrep(regNamesRaw{j},'_','-');
%        regNames{j} = nameSplit;
%    end

SN=4;


selidx = [1 3 4 5];
currCoM=pathCoM(selidx,:,1);
pathNamesRawSel = pathNamesRaw(selidx);
N = length(pathNamesRawSel);
%thickCtrl = AllResults.Thickness.Mean(CohortsIndexs==1,currNetwork);


data=AllthicknessAtPath(:,:,1);

ctrldata=CtrlthicknessAtPath(:,:,1);

ctrldata_mean=mean(ctrldata,1);
ctrldata_std=std(ctrldata,1);

data_w=(data-ctrldata_mean)./ctrldata_std;

path_yesAD = data(demoT.LBD0_LBDAD1 == 1,selidx);
path_noAD = data(demoT.LBD0_LBDAD1 == 0,selidx);

path_yesAD_w = data_w(demoT.LBD0_LBDAD1 == 1,selidx);
path_noAD_w = data_w(demoT.LBD0_LBDAD1 == 0,selidx);


%covMatCtrl = nan(N,N);
covMat_yesAD = nan(N,N);
covMat_noAD = nan(N,N);

%cmpCovTau_gt_Ctrl = nan(N,N);
%cmpCovTDP_gt_Ctrl = nan(N,N);
cmpCovYes_gt_No = nan(N,N);
cmpCovNo_gt_Yes = nan(N,N);


for i = 1:N
    for j = 1:N
      
        if(i~=j)
        NYesAD=sum(~isnan(path_yesAD(:,i)) & ~isnan(path_yesAD(:,j)));
        NNoAD=sum(~isnan(path_noAD(:,i)) & ~isnan(path_noAD(:,j)));

            if(NYesAD > 3 && NNoAD >3)
                

            pathNamesRawSel{i}
            pathNamesRawSel{j}

            covMat_yesAD(i,j)=corr(path_yesAD(:,i),path_yesAD(:,j),'Rows','complete');
            covMat_noAD(i,j)=corr(path_noAD(:,i),path_noAD(:,j),'Rows','complete');
                
            cmpCovYes_gt_No(i,j) = TestCorrSig(covMat_noAD(i,j),covMat_yesAD(i,j),NNoAD,NYesAD)<.05;
            cmpCovNo_gt_Yes(i,j) = TestCorrSig(covMat_yesAD(i,j),covMat_noAD(i,j),NYesAD,NNoAD)<.05;

            corr(path_yesAD(:,i),path_yesAD(:,j),'Rows','complete')
            NYesAD
            corr(path_noAD(:,i),path_noAD(:,j),'Rows','complete')
            NNoAD
            [r z]=TestCorrSig(covMat_noAD(i,j),covMat_yesAD(i,j),NNoAD,NYesAD)

            [r z]=TestCorrSig(covMat_yesAD(i,j),covMat_noAD(i,j),NYesAD,NNoAD)
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

saveName='LBD'
minPath=min([path_yesAD; path_noAD],[],1,'omitnan');
maxPath=max([path_yesAD; path_noAD]-minPath+0.0015,[],1,'omitnan');
MakerVec=mean(path_yesAD,1,'omitnan');
MakerVec=2*(MakerVec-minPath)./maxPath;


H=figure(2)
imagesc(covMat_yesAD)
title('LBD-yesAD')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',pathNamesRawSel)
set(gca,'ytick',[1:N],'yticklabel',pathNamesRawSel)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ 'CovMat_yesAD.png']),'-dpng','-r400');

cRange = [0 1];
MakerVec=ones(1,SN);
colorVec=ones(1,SN);
LabelNames=pathNamesRawSel;
pSelect=1;
plotNetwork3(covMat_yesAD, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVec,panelAll,pSelect,colorVec,4)


H=figure(3)
imagesc(covMat_noAD)
title('LBD-noAD')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',pathNamesRawSel)
set(gca,'ytick',[1:N],'yticklabel',pathNamesRawSel)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ 'CovMat_noAD.png']),'-dpng','-r400');
pSelect=2;
plotNetwork3(covMat_noAD, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVec,panelAll,pSelect,colorVec,4)


H=figure(6)
imagesc(cmpCovYes_gt_No)
title('yesAD > noAD')
caxis([0 1])
colormap('jet')

set(gca,'xtick',[1:N],'xticklabel',pathNamesRawSel)
set(gca,'ytick',[1:N],'yticklabel',pathNamesRawSel)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TAUGTTDP.png']),'-dpng','-r400');
pSelect=3;
plotNetwork3(cmpCovYes_gt_No, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVec,panelAll,pSelect,colorVec,0)


H=figure(7)
imagesc(cmpCovNo_gt_Yes)
title('noAD > yesAD')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',pathNamesRawSel)
set(gca,'ytick',[1:N],'yticklabel',pathNamesRawSel)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_1TDPGTTAU.png']),'-dpng','-r400');
pSelect=4;
MakerVec=mean(path_noAD,1,'omitnan');
MakerVec=2*(MakerVec-minPath)/maxPath;
plotNetwork3(cmpCovNo_gt_Yes, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVec,panelAll,pSelect,colorVec,0)

print(H2,fullfile(baseSaveDir,[ saveName '_SigPlots.png']),'-dpng','-r400');