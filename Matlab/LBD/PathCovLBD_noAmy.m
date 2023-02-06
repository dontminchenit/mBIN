
%regNames=cell(1,SN);
    
%    for j = 1:SN
%        nameSplit = strrep(regNamesRaw{j},'_','-');
%        regNames{j} = nameSplit;
%    end



currCoM=pathCoM(:,:,1);


N = length(pathNamesRaw);
%thickCtrl = AllResults.Thickness.Mean(CohortsIndexs==1,currNetwork);
path_yesAD = pathDataGM(LBD_yesADIndx,:);
path_noAD = pathDataGM(LBD_noADIndx,:);



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
                
            covMat_yesAD(i,j)=corr(path_yesAD(:,i),path_yesAD(:,j),'Rows','complete');
            covMat_noAD(i,j)=corr(path_noAD(:,i),path_noAD(:,j),'Rows','complete');
                
            cmpCovYes_gt_No(i,j) = TestCorrSig(covMat_noAD(i,j),covMat_yesAD(i,j),NNoAD,NYesAD)<.05;
            cmpCovNo_gt_Yes(i,j) = TestCorrSig(covMat_yesAD(i,j),covMat_noAD(i,j),NYesAD,NNoAD)<.05;
            end
        end
        
    end
end

covMat_yesAD(2,:)=[];
covMat_noAD(2,:)=[];
cmpCovYes_gt_No(2,:)=[];
cmpCovNo_gt_Yes(2,:)=[];

covMat_yesAD(:,2)=[];
covMat_noAD(:,2)=[];
cmpCovYes_gt_No(:,2)=[];
cmpCovNo_gt_Yes(:,2)=[];
currCoM(2,:)=[]

pathNamesAdj=pathNamesRaw;
pathNamesAdj(2)=[];
H2=figure(10)
clf
panelAll = panel();
panelAll.pack(2,4);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

saveName='LBD'

path_yesAD_exp=path_yesAD;
path_noAD_exp=path_noAD;
minPath=min([path_yesAD_exp; path_noAD_exp],[],1,'omitnan');
maxPath=max([path_yesAD_exp; path_noAD_exp]-minPath+0.0015,[],1,'omitnan');
MakerVecYes=mean(path_yesAD_exp,1,'omitnan');
MakerVecYes=3*(MakerVecYes-minPath)./maxPath;
MakerVecNo=mean(path_noAD_exp,1,'omitnan');
MakerVecNo=3*(MakerVecNo-minPath)./maxPath;
MakerVecYes(2)=[]
MakerVecNo(2)=[]

H=figure(2)
imagesc(covMat_yesAD)
title('LBD-yesAD')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',pathNamesAdj)
set(gca,'ytick',[1:N],'yticklabel',pathNamesAdj)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ 'CovMat_yesAD.png']),'-dpng','-r400');

cRange = [0 1];
%MakerVec=ones(1,SN);
colorVec=ones(1,SN-1);
LabelNames=pathNamesAdj;
pSelectCol=1;
plotNetwork3(covMat_yesAD, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecYes,panelAll,pSelectCol,colorVec,4)


H=figure(3)
imagesc(covMat_noAD)
title('LBD-noAD')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',pathNamesAdj)
set(gca,'ytick',[1:N],'yticklabel',pathNamesAdj)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ 'CovMat_noAD.png']),'-dpng','-r400');
pSelectCol=2;
plotNetwork3(covMat_noAD, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecNo,panelAll,pSelectCol,colorVec,4)


H=figure(6)
imagesc(cmpCovYes_gt_No)
title('yesAD > noAD')
caxis([0 1])
colormap('jet')

set(gca,'xtick',[1:N],'xticklabel',pathNamesAdj)
set(gca,'ytick',[1:N],'yticklabel',pathNamesAdj)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_TAUGTTDP.png']),'-dpng','-r400');
pSelectCol=3;
plotNetwork3(cmpCovYes_gt_No, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecYes,panelAll,pSelectCol,colorVec,0)


H=figure(7)
imagesc(cmpCovNo_gt_Yes)
title('noAD > yesAD')
caxis([0 1])
colormap('jet')
set(gca,'xtick',[1:N],'xticklabel',pathNamesAdj)
set(gca,'ytick',[1:N],'yticklabel',pathNamesAdj)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_1TDPGTTAU.png']),'-dpng','-r400');
pSelectCol=4;
%MakerVecYes=mean(path_noAD,1,'omitnan');
%MakerVecYes=2*(MakerVecYes-minPath)/maxPath;
plotNetwork3(cmpCovNo_gt_Yes, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecNo,panelAll,pSelectCol,colorVec,0)

print(H2,fullfile(baseSaveDir,[ saveName '_SigPlots.png']),'-dpng','-r400');