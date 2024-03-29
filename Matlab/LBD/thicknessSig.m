%LoadData

        i=4;%thickness
        j=6;%mean
        
        H2=figure(2)
        clf
        panelAll = panel();
        panelAll.pack(4,6);

        panelAll.de.margin = 1;
        panelAll.marginbottom = 1;
        panelAll.marginleft = -10;
        panelAll.margin = [0 0 0 0];
        
        LBDsigThickness=struct;
        %measure = AllResults.(graphNames{i}).(measureNames{j});
        %if(j==1)%if betweenness, then we need to normalize
        %    measure = measure./((numLab-1)*(numLab-2));
        %end
        
        %we have the measure for each subgroup
        thickCtrl = CtrlResults.Thickness.Mean(:,currNetwork);
        thickyesAD = AllResults.Thickness.Mean(demoT.LBD0_LBDAD1 == 1,currNetwork);
        thicknoAD = AllResults.Thickness.Mean(demoT.LBD0_LBDAD1 == 0,currNetwork);
        GroupNames={'Ctrls','YesAD','NoAD'}

        GroupMeasures = cell(1,3);
        GroupMeasures{1}=thickCtrl;
        GroupMeasures{2}=thickyesAD;
        GroupMeasures{3}=thicknoAD;
        
        %saveDirSub = fullfile(saveDirBase,graphNames{i},measureNames{j});
        %mkdir(saveDir)
        %CRange=([0 prctile(measure(:),95)]);
        

        %We make these comparisons between groups
        testG1All=[2 3 2 3];
        testG2All=[1 1 3 2];
        %plotg1=[1 3 5 -1];
        %plotg2=[2 4 6 7];
        plotg1=[0 0 3 0];
        plotg2=[1 2 3 4];
            for g = 1:2
                G1 = testG1All(g);
                G2 = testG2All(g);
                testName=[GroupNames{G1} 'vs' GroupNames{G2}];
                TTestsAll_pVal = nan(1,numLab);
                TTestsAll_tstat = nan(1,numLab);
                %do t-test on all labels
                for t = 1:numLab
                    [h1 p ci stats]=ttest2(GroupMeasures{G1}(:,t),GroupMeasures{G2}(:,t),'Tail','left');
                    TTestsAll_pVal(t)=p;
                    TTestsAll_tstat(t) = stats.tstat;
                end
                
                pThresh=.05
                %correct for multiple comparison
                notNaN = ~isnan(TTestsAll_pVal);
                [padj_fdr,alpha_fdr] = multicmp(TTestsAll_pVal(notNaN),'fdr',pThresh);
                
                allPadj = {padj_fdr};
                allPadjName = {'FDR_Benjamini-Hochberg',...
                    'FWE_Hochberg-Bonferroni','FWE_Holm-Bonferroni'};
                
                for p = 1%:3
                    data = TTestsAll_pVal
                    data(~notNaN) = 1;
                    data(notNaN) = allPadj{p};
                    CRange=[0 1];
                    saveName = testName;

                    LBDsigThickness.(GroupNames{G1})=data;
                    %saveDir = fullfile(saveDirSub,testName,allPadjName{p});
                    %plotBrain2(labelmesh,data,LabelLUT,CRange,saveDir,saveName,0.001,panelAll,plotg2(g))
                    %plotBrain3(Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0.001,panelAll,plotg2(g))
                    %if(g~=3)
                    %plotBrain3Points(NetworkDataGeneral.Schaefer400x7.CoM,NetworkDataGeneral.Schaefer400x7.GII,data,pThresh,panelAll,plotg2(g),1,10)
                    %print(H2,fullfile(baseSaveDir,'thickCmp_p05.tif'),'-dpng','-r400');
                    %end
                end

                data = TTestsAll_tstat;
                CRange=[-4 4];
                
                if(g==1 || g==2) %if comparing to controls we only care about negative tstats (patients < HC)
                    data(data>0) = 0;
                end
                
                saveName = testName;
                %saveDir = fullfile(saveDirSub,testName,'tstat');
                %plotBrain2(labelmesh,data,LabelLUT,CRange,saveDir,saveName,0,panelAll,plotg1(g))
                %plotBrain3(Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,plotg1(g))
                if(g==3)
                    %plotBrain3Points(NetworkDataGeneral.Schaefer400x7.CoM,NetworkDataGeneral.Schaefer400x7.GII,data,0,panelAll,plotg1(g),2,10)
                    %print(H2,fullfile(baseSaveDir,'thickCmp2.tif'),'-dpng','-r400');
                end
            end
                
%            print('-dtiff','-r500',fnames)
%            saveas(H2,fullfile(saveDirBase,'fig1.tif'))
            
            %print(H2,fullfile(baseSaveDir,'thickCmp1.tif'),'-dpng','-r400');