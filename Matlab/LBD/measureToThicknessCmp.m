networkNames={'DorsAttn','VentAttn','Limbic','SomMot','Default','Cont','_Vis_','All'};
measureName={'MMSE','ANIMALS', 'FLUENCY', 'BNT', 'PVLT', 'DF', 'DB', 'MMSEbaseline'};

baseSaveDir=fullfile(currBase,'Results_1_13_2023');
mkdir(baseSaveDir)


        H2=figure(10)
        %clf
        panelAll = panel();
        panelAll.pack(2,4);

        panelAll.de.margin = 1;
        panelAll.marginbottom = 1;
        panelAll.marginleft = -10;
        panelAll.margin = [0 0 0 0];
        
for ii=4:5
    if(ii==8)
        currNetwork=true(1,length(NetworkDataGeneral.Schaefer400x7.LUT.Label_Name2));
    else
        currNetwork=contains(NetworkDataGeneral.Schaefer400x7.LUT.Label_Name2,networkNames{ii})';
    end

    saveName=[networkNames{ii} ];

    regNamesRaw = NetworkDataGeneral.Schaefer400x7.LUT.Label_Name2(currNetwork);
    SN= length(regNamesRaw);
    regNames=cell(1,SN);

    for j = 1:SN
        nameSplit = strrep(regNamesRaw{j},'_','-');
        regNames{j} = nameSplit;
    end


    currCoM=NetworkDataGeneral.Schaefer400x7.CoM(currNetwork,:);
    N = sum(currNetwork);

    %currCoM=CoM;
    thickCtrl = CtrlResults.Thickness.Mean(:,currNetwork);

%    yesADSig=LBDsigThickness.YesAD < .05;
%    noADSig=LBDsigThickness.NoAD < .05;
%    thickyesAD = AllResults.Thickness.Mean(demoT.LBD0_LBDAD1 == 1,currNetwork & yesADSig);
%    thicknoAD = AllResults.Thickness.Mean(demoT.LBD0_LBDAD1 == 0,currNetwork & noADSig);

    thickyesAD = AllResults.Thickness.Mean(demoT.LBD0_LBDAD1 == 1,currNetwork);
    thicknoAD = AllResults.Thickness.Mean(demoT.LBD0_LBDAD1 == 0,currNetwork);


    Ctrl_Mean = mean(thickCtrl,1,'omitnan');
    Ctrl_SD = std(thickCtrl,1,'omitnan');

    YesAD_W=nan(size(thickyesAD));
    for i = 1:size(thickyesAD,1)
        YesAD_W(i,:)= (thickyesAD(i,:) - Ctrl_Mean)./Ctrl_SD;
    end
    
    %thickyesAD=YesAD_W;

    NoAD_W=nan(size(thicknoAD));
    for i = 1:size(thicknoAD,1)
        NoAD_W(i,:)= (thicknoAD(i,:) - Ctrl_Mean)./Ctrl_SD;
    end

    %thicknoAD=NoAD_W;

    Ctrl_W=nan(size(thicknoAD));
    for i = 1:size(thicknoAD,1)
        Ctrl_W(i,:)= (thickCtrl(i,:) - Ctrl_Mean)./Ctrl_SD;
    end
    %thickCtrl=Ctrl_W;

    mean_YesAD_w = mean(YesAD_W,1,'omitnan');
    mean_NoAD_w = mean(NoAD_W,1,'omitnan');

    covMatCtrl = nan(N,N);
    covMatyesAD = nan(N,N);
    covMatnoAD = nan(N,N);

    cmpCovYes_gt_Ctrl = nan(N,N);
    cmpCovNo_gt_Ctrl = nan(N,N);

    cmpCovYes_lt_Ctrl = nan(N,N);
    cmpCovNo_lt_Ctrl = nan(N,N);


    cmpCovYes_gt_No = nan(N,N);
    cmpCovNo_gt_Yes = nan(N,N);

    pthresh = .01;
    for i = 1:N
        for j = 1:N

            if(i~=j)
                NCtrl=sum(~isnan(thickCtrl(:,i)) & ~isnan(thickCtrl(:,j)));
                NYes=sum(~isnan(thickyesAD(:,i)) & ~isnan(thickyesAD(:,j)));
                NNo=sum(~isnan(thicknoAD(:,i)) & ~isnan(thicknoAD(:,j)));

                if(NCtrl > 3)
                    covMatCtrl(i,j)=corr(thickCtrl(:,i),thickCtrl(:,j),'Rows','complete');
                end
                if(NYes > 3)
                    covMatyesAD(i,j)=corr(thickyesAD(:,i),thickyesAD(:,j),'Rows','complete');
                end
                if(NNo > 3)
                    covMatnoAD(i,j)=corr(thicknoAD(:,i),thicknoAD(:,j),'Rows','complete');
                end
%
                if(NCtrl > 3 && NYes > 3)
                    cmpCovYes_gt_Ctrl(i,j) = TestCorrSig(covMatCtrl(i,j),covMatyesAD(i,j),NCtrl,NYes)<pthresh;
                    cmpCovYes_lt_Ctrl(i,j) = TestCorrSig(covMatyesAD(i,j),covMatCtrl(i,j),NYes,NCtrl)<pthresh;
                end
                if(NCtrl > 3 && NNo > 3)
                    cmpCovNo_gt_Ctrl(i,j) = TestCorrSig(covMatCtrl(i,j),covMatnoAD(i,j),NCtrl,NNo)<pthresh;
                    cmpCovNo_lt_Ctrl(i,j) = TestCorrSig(covMatnoAD(i,j),covMatCtrl(i,j),NNo,NCtrl)<pthresh;
                end
                if(NYes > 3 && NNo > 3)
                    cmpCovYes_gt_No(i,j) = TestCorrSig(covMatnoAD(i,j),covMatyesAD(i,j),NNo,NYes)<pthresh;
                    cmpCovNo_gt_Yes(i,j) = TestCorrSig(covMatyesAD(i,j),covMatnoAD(i,j),NYes,NNo)<pthresh;
                end
            end

        end
    end

efficiency_bin(cmpCovYes_gt_Ctrl)
efficiency_bin(cmpCovNo_gt_Ctrl)
efficiency_bin(covMatCtrl)
ROIyes=sum(cmpCovYes_gt_Ctrl,'omitnan')>-1;
ROIno=sum(cmpCovNo_gt_Ctrl,'omitnan')>-1;
key_thickness_yes=mean(YesAD_W(:,ROIyes),2);
key_thickness_no=mean(NoAD_W(:,ROIno),2);
for m = 8
switch m
    case 1
    measureData = measuresT_matched.MMSEANNCHANGE;
    case 2
    measureData = measuresT_matched.ANIMALSANNCHANGE;
    case 3
    measureData = measuresT_matched.FLUENCYANNCHANGE;
    case 4
    measureData = measuresT_matched.BNTANNCHANGE;
    case 5
    measureData = measuresT_matched.PVLTANNCHANGE;
    case 6
    measureData = measuresT_matched.DFANNCHANGE;
    case 7
    measureData = measuresT_matched.DBANNCHANGE;
    case 8
    measureData = measuresT_matched.MMSE_1;
end

measure_yes = measureData(demoT.LBD0_LBDAD1 == 1);
measure_no = measureData(demoT.LBD0_LBDAD1 == 0);
H3=figure(61)
scatter(key_thickness_yes,measure_yes)
ylabel([measureName{m} ' Annualized Change)'])
xlabel('Mean Thickness(W-score) at Significant Covariance ROIs')
[r p]=corr(key_thickness_yes,measure_yes,'rows','complete');
title(['LBD+AD, ' saveName ', r=' num2str(r) ' p=' num2str(p)]);
lsline
print(H3,fullfile(baseSaveDir,['measureCmp_YESAD_' saveName '_' measureName{m} '.tif']),'-dpng','-r400');

H3=figure(62)
scatter(key_thickness_no,measure_no)
ylabel([measureName{m} ' Annualized Change)'])
xlabel('Mean Thickness(W-score) at Significant Covariance ROIs')
[r p]=corr(key_thickness_no,measure_no,'rows','complete');
title(['LBD-AD, ' saveName ', r=' num2str(r) ' p=' num2str(p)]);
lsline
print(H3,fullfile(baseSaveDir,['measureCmp_NOAD_' saveName '_' measureName{m} '.tif']),'-dpng','-r400');

    %pthresh=.05;
    %cmpCovNo_gt_Ctrl=mtxAdjustMultCmp(cmpCovNo_gt_Ctrl,pthresh);
    %cmpCovNo_gt_Ctrl=mtxAdjustMultCmp(cmpCovNo_gt_Ctrl,pthresh);
    %cmpCovYes_gt_No=mtxAdjustMultCmp(cmpCovYes_gt_No,pthresh);
    %cmpCovNo_gt_Yes=mtxAdjustMultCmp(cmpCovNo_gt_Yes,pthresh);
end
end