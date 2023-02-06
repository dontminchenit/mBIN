thickcovMatCtrl = nan(N,N);
thickcovMatTau = nan(N,N);
thickcovMatTDP = nan(N,N);

cmpCovTau_gt_Ctrl = nan(N,N);
cmpCovTDP_gt_Ctrl = nan(N,N);
cmpCovTau_gt_Ctrl_z = nan(N,N);
cmpCovTDP_gt_Ctrl_z = nan(N,N);
cmpCovTau_gt_TDP = nan(N,N);
cmpCovTDP_gt_Tau = nan(N,N);

thickCtrl = matchedThicknessDataCat(CohortsIndexs==1,:);
thickTau = matchedThicknessDataCat(CohortsIndexs==2,:);
thickTDP = matchedThicknessDataCat(CohortsIndexs==3,:);


for i = 1:N
    for j = 1:N
      
        if(i~=j)
        NCtrl=sum(~isnan(thickCtrl(:,i)) & ~isnan(thickCtrl(:,j)));
        NTau=sum(~isnan(thickTau(:,i)) & ~isnan(thickTau(:,j)));
        NTDP=sum(~isnan(thickTDP(:,i)) & ~isnan(thickTDP(:,j)));

        if(NCtrl > 3)
            thickcovMatCtrl(i,j)=corr(thickCtrl(:,i),thickCtrl(:,j),'Rows','complete');
        end
        if(NTau > 3)
            thickcovMatTau(i,j)=corr(thickTau(:,i),thickTau(:,j),'Rows','complete');
        end
        if(NTDP > 3)
            thickcovMatTDP(i,j)=corr(thickTDP(:,i),thickTDP(:,j),'Rows','complete');
        end
        
        if(NCtrl > 3 && NTau > 3)
            [P z] = TestCorrSig(thickcovMatCtrl(i,j),thickcovMatTau(i,j),NCtrl,NTau);
            cmpCovTau_gt_Ctrl(i,j) = P<.05;
            cmpCovTau_gt_Ctrl_z(i,j) = z;
        end

        if(NCtrl > 3 && NTDP > 3)
            [P z] = TestCorrSig(thickcovMatCtrl(i,j),thickcovMatTDP(i,j),NCtrl,NTDP);
            cmpCovTDP_gt_Ctrl(i,j) = P<.05;
            cmpCovTDP_gt_Ctrl_z(i,j) = z;
        end


        end
        
    end
end

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
MakerVec=ones(1,N);
colorVec=ones(1,N);
LabelNames=regNames;
%pSelect=1;

plotNetwork3(cmpCovTau_gt_Ctrl, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVec,panelAll,2,colorVec,0)
plotNetwork3(cmpCovTDP_gt_Ctrl, schaefer400x7SurfaceMeshFromGII, currCoM,LabelNames,cRange,MakerVec,panelAll,1,colorVec,0)
print(H2,fullfile(baseSaveDir,[ saveName '_PathDiff_W.png']),'-dpng','-r400');


figure(10)
subplot(1,3,1)
imagesc(thickcovMatCtrl)
subplot(1,3,3)
imagesc(thickcovMatTau)
subplot(1,3,2)
imagesc(thickcovMatTDP)
