networkNames={'DorsAttn','VentAttn','Limbic','SomMot','Default','Cont','_Vis_','All'};

H2=figure(10)
%clf
panelAll = panel();
panelAll.pack(2,4);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

NN = 1000;

NetworkNames=cell(NN,1);
ROIName1=cell(NN,1);
ROIName1_Path=cell(NN,1);
ROIName2=cell(NN,1);
ROIName2_Path=cell(NN,1);
r_YesAD=nan(NN,1);
r_NoAD=nan(NN,1);
N_YesAD=nan(NN,1);
N_NoAD=nan(NN,1);
z_YesAD_gt_NoAD=nan(NN,1);
z_NoAD_gt_YesAD=nan(NN,1);
p_YesAD_gt_NoAD=nan(NN,1);
p_NoAD_gt_YesAD=nan(NN,1);
SigEdgesT = table(NetworkNames,ROIName1,ROIName1_Path,ROIName2,ROIName2_Path,r_YesAD,r_NoAD,N_YesAD,N_NoAD,z_YesAD_gt_NoAD,z_NoAD_gt_YesAD,p_YesAD_gt_NoAD,p_NoAD_gt_YesAD);
SigEdgesCounter=0;

for ii=5
    if(ii==8)
        currNetwork=true(1,length(NetworkDataGeneral.Schaefer400x7.LUT.Label_Name2));
    else
        currNetwork=contains(NetworkDataGeneral.Schaefer400x7.LUT.Label_Name2,networkNames{ii})';
    end

    saveName=[networkNames{ii} ];

    pathLUTSorted=sortrows(pathLUT);
    pathNames=pathLUTSorted.PathName(currNetwork); 
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
                    if(covMatCtrl(i,j) < .1)
                        covMatCtrl(i,j)=nan;
                        NCtrl=0;
                    end
                end
                if(NYes > 3)
                    covMatyesAD(i,j)=corr(thickyesAD(:,i),thickyesAD(:,j),'Rows','complete');
                    if(covMatyesAD(i,j) < .1)
                        covMatyesAD(i,j)=nan;
                    end

                end
                if(NNo > 3)
                    covMatnoAD(i,j)=corr(thicknoAD(:,i),thicknoAD(:,j),'Rows','complete');
                    if(covMatnoAD(i,j) < .1)
                        covMatnoAD(i,j)=nan;
                    end
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

                if((cmpCovYes_gt_No(i,j) == 1 || cmpCovNo_gt_Yes(i,j) == 1) && i<j) %record significant edges
                    SigEdgesCounter=SigEdgesCounter+1;
                    ss=SigEdgesCounter;

                    SigEdgesT.NetworkNames{ss}=networkNames{ii};
                    SigEdgesT.ROIName1{ss}=regNames{i};
                    SigEdgesT.ROIName1_Path{ss}=pathNames{i};
                    SigEdgesT.ROIName2{ss}=regNames{j};
                    SigEdgesT.ROIName2_Path{ss}=pathNames{j};
                    SigEdgesT.r_YesAD(ss)=corr(thickyesAD(:,i),thickyesAD(:,j),'Rows','complete');
                    SigEdgesT.r_NoAD(ss)=corr(thicknoAD(:,i),thicknoAD(:,j),'Rows','complete');
                    SigEdgesT.N_YesAD(ss)=NYes;
                    SigEdgesT.N_NoAD(ss)=NNo;


                    [p z] = TestCorrSig(covMatnoAD(i,j),covMatyesAD(i,j),NNo,NYes);
                    SigEdgesT.p_YesAD_gt_NoAD(ss)=p;
                    SigEdgesT.z_YesAD_gt_NoAD(ss)=z;
                    [p z] = TestCorrSig(covMatyesAD(i,j),covMatnoAD(i,j),NYes,NNo);
                    SigEdgesT.p_NoAD_gt_YesAD(ss)=p;
                    SigEdgesT.z_NoAD_gt_YesAD(ss)=z;

                end

            end

        end
    end
    
    writetable(SigEdgesT,fullfile(baseSaveDir,'SigEdges_20230106.csv'))

    efficiency_bin(cmpCovYes_gt_Ctrl)
    efficiency_bin(cmpCovNo_gt_Ctrl)
    efficiency_bin(covMatCtrl)

    ROIyes=sum(cmpCovYes_gt_Ctrl,'omitnan')>0;
    ROIno=sum(cmpCovNo_gt_Ctrl,'omitnan')>0;

    key_thickness_yes=mean(thickyesAD(:,ROIyes),2);
    key_thickness_no=mean(thicknoAD(:,ROIno),2);

    measureData = measuresT_matched.BNTANNCHANGE;
    measure_yes = measureData(demoT.LBD0_LBDAD1 == 1);
    measure_no = measureData(demoT.LBD0_LBDAD1 == 0);


    H3=figure(61)
    scatter(key_thickness_yes,measure_yes)
    ylabel('MMSE Annualized Change)')
    xlabel('Mean Thickness at Significant Covariance ROIs')
    [r p]=corr(key_thickness_yes,measure_yes,'rows','complete');
    title(['LBD+AD, ' saveName ', r=' num2str(r) ' p=' num2str(p)]);
    lsline
    %print(H3,fullfile(baseSaveDir,['measureCmpMMSE_YES' saveName '_' measureName '.tif']),'-dpng','-r400');
    print(H3,fullfile(baseSaveDir,['measureCmpMMSE_YES' saveName '_' '.tif']),'-dpng','-r400'); % originally commented out
    

    H3=figure(62)
    scatter(key_thickness_no,measure_no)
    ylabel('MMSE Annualized Change)')
    xlabel('Mean Thickness at Significant Covariance ROIs')
    [r p]=corr(key_thickness_no,measure_no,'rows','complete');
    title(['LBD-AD, ' saveName ', r=' num2str(r) ' p=' num2str(p)]);
    lsline
    %print(H3,fullfile(baseSaveDir,['measureCmpMMSE_NO' saveName '_' measureName '.tif']),'-dpng','-r400');
    print(H3,fullfile(baseSaveDir,['measureCmpMMSE_NO' saveName '_' '.tif']),'-dpng','-r400'); % originally commented out

    %pthresh=.05;
    %cmpCovNo_gt_Ctrl=mtxAdjustMultCmp(cmpCovNo_gt_Ctrl,pthresh);
    %cmpCovNo_gt_Ctrl=mtxAdjustMultCmp(cmpCovNo_gt_Ctrl,pthresh);
    %cmpCovYes_gt_No=mtxAdjustMultCmp(cmpCovYes_gt_No,pthresh);
    %cmpCovNo_gt_Yes=mtxAdjustMultCmp(cmpCovNo_gt_Yes,pthresh);

    H3=figure(50)

    covMatyesADnz=covMatyesAD;
    covMatyesADnz(covMatyesADnz<0)=0;
    covMatnoADnz=covMatnoAD;
    covMatnoADnz(covMatnoADnz<0)=0;
    covMatCtrlnz=covMatCtrl;
    covMatCtrlnz(covMatCtrlnz<0)=0;

    degYes=sum(covMatyesADnz,'omitnan');
    degNo=sum(covMatnoADnz,'omitnan');
    degCtrl=sum(covMatCtrlnz,'omitnan');

    H3=figure(1)
    boxchart([degCtrl' degYes' degNo' ],'Notch','on');
    [H pval]=ttest2(degYes,degNo);
    title([saveName ' ,p=' num2str(pval)], 'Interpreter', 'none');
    ylabel('Degree')
    set(gca,'XTickLabel',{'Ctrl','Yes AD', 'No AD'})
    print(H3,fullfile(baseSaveDir,['nonZero_degCmp_' saveName '.tif']),'-dpng','-r400');

    % thickCtrl = CtrlResults.Thickness.Mean(:,currNetwork);
    % thickyesAD = AllResults.Thickness.Mean(demoT.LBD0_LBDAD1 == 1,currNetwork);
    % thicknoAD = AllResults.Thickness.Mean(demoT.LBD0_LBDAD1 == 0,currNetwork);
    
    H3=figure(2)
    subplot(1,2,1)
    scatter(degYes,mean(thickyesAD,1))
    [h p]=corrcoef(degYes,mean(thickyesAD,1))
    title([saveName ' ,r=' num2str(h(1,2)) ' ,p=' num2str(p(1,2))], 'Interpreter', 'none');
    lsline
    xlabel('Degree')
    ylabel('Thickness') 
    subplot(1,2,2)
    scatter(degNo,mean(thicknoAD,1))
    [h p]=corrcoef(degNo,mean(thicknoAD,1))
    title([saveName ' ,r=' num2str(h(1,2)) ' ,p=' num2str(p(1,2))], 'Interpreter', 'none');
    lsline
    xlabel('Degree')
    ylabel('Thickness')
    print(H3,fullfile(baseSaveDir,['nonZero_degCorr_' saveName '.tif']),'-dpng','-r400');
    
    plotON=1;
    if(plotON)

        plotMtx=0

        if(plotMtx)
            H=figure(1)
            imagesc(covMatCtrl)
            title('controls')
            caxis([0 1])
            colormap('jet')
            set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
            set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
            set(gca,'TickLabelInterpreter','none','FontSize',5)
            print(H,fullfile(baseSaveDir,[ saveName '_ctrl.png']),'-dpng','-r400');


            H=figure(2)
            imagesc(covMatyesAD)
            title('LBD-yesAD')
            caxis([0 1])
            colormap('jet')
            set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
            set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
            set(gca,'TickLabelInterpreter','none','FontSize',5)
            print(H,fullfile(baseSaveDir,[ saveName '_yesAD.png']),'-dpng','-r400');


            H=figure(3)
            imagesc(covMatnoAD)
            title('LBD-noAD')
            caxis([0 1])
            colormap('jet')
            set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
            set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
            set(gca,'TickLabelInterpreter','none','FontSize',5)
            print(H,fullfile(baseSaveDir,[ saveName '_noAD.png']),'-dpng','-r400');
            %

            H=figure(4)
            imagesc(cmpCovYes_gt_Ctrl)
            title('YesAD > Ctrl')
            caxis([0 1])
            colormap('jet')
            set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
            set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
            set(gca,'TickLabelInterpreter','none','FontSize',5)
            print(H,fullfile(baseSaveDir,[ saveName '_YesGTCTRL.png']),'-dpng','-r400');
        end

        cRange = [0 1];
        MakerVec=ones(1,sum(currNetwork));

        thick_yesAD_exp=YesAD_W;
        thick_noAD_exp=NoAD_W;

        minThick=min(min([thick_yesAD_exp; thick_noAD_exp],[],1,'omitnan'));
        maxThick=max(max([thick_yesAD_exp; thick_noAD_exp]-minThick+0.0015,[],1,'omitnan'));

        MakerVecYes=rescale(mean(thick_yesAD_exp,1,'omitnan'),.35,.9);
        MakerVecYes=1./MakerVecYes;
        %MakerVecYes=1./((MakerVecYes-minThick)./maxThick) -1;

        MakerVecNo=rescale(mean(thick_noAD_exp,1,'omitnan'),.35,.9);
        MakerVecNo=1./MakerVecNo;
        %MakerVecNo=1./((MakerVecNo-minThick)./maxThick) -1;

        colorVecYes=mean(thick_yesAD_exp,1,'omitnan');
        colorVecNo=mean(thick_noAD_exp,1,'omitnan');

        colorVecYes=ones(1,sum(currNetwork));
        colorVecNo=ones(1,sum(currNetwork));
        
        MakerVecYes=ones(1,sum(currNetwork));
        MakerVecNo=ones(1,sum(currNetwork));
        displayType=0;

        LabelNames=regNamesRaw;
        
        pSelectcol=1;
        plotNetwork3(cmpCovYes_gt_Ctrl, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecYes,panelAll,pSelectcol,colorVecYes,displayType)

        if(plotMtx)

            H=figure(5)
            imagesc(cmpCovNo_gt_Ctrl)
            title('NoAD > Ctrl')
            caxis([0 1])
            colormap('jet')
            set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
            set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
            set(gca,'TickLabelInterpreter','none','FontSize',5)
            print(H,fullfile(baseSaveDir,[ saveName '_NOGTCTRL.png']),'-dpng','-r400');
        end

        pSelectcol=2;
        plotNetwork3(cmpCovNo_gt_Ctrl, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecNo,panelAll,pSelectcol,colorVecNo,displayType)
        
        if(plotMtx)

            H=figure(6)
            imagesc(cmpCovYes_gt_No)
            title('YesAD > NoAD')
            caxis([0 1])
            colormap('jet')
            set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
            set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
            set(gca,'TickLabelInterpreter','none','FontSize',5)
            print(H,fullfile(baseSaveDir,[ saveName '_YesGTNO.png']),'-dpng','-r400');
        end
        pSelectcol=3;
        plotNetwork3(cmpCovYes_gt_No, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecYes,panelAll,pSelectcol,colorVecYes,displayType)


        if(plotMtx)

            H=figure(7)
            imagesc(cmpCovNo_gt_Yes)
            title('NoAD > YesAD')
            caxis([0 1])
            colormap('jet')
            set(gca,'xtick',[1:N],'xticklabel',regNamesRaw)
            set(gca,'ytick',[1:N],'yticklabel',regNamesRaw)
            set(gca,'TickLabelInterpreter','none','FontSize',5)
            print(H,fullfile(baseSaveDir,[ saveName '_NOGTYES.png']),'-dpng','-r400');
        end

        pSelectcol=4;
        plotNetwork3(cmpCovNo_gt_Yes, NetworkDataGeneral.Schaefer400x7.GII, currCoM,LabelNames,cRange,MakerVecNo,panelAll,pSelectcol,colorVecNo,displayType)

        print(H2,fullfile(baseSaveDir,[ saveName '_SigPlots.png']),'-dpng','-r400');
    end
end