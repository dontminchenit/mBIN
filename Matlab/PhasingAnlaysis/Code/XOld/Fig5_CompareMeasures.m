LoadDataSarah

for i = 2
    
plotg1=[1 2 3 4];
        H2=figure(2)
        clf
        panelAll = panel();
        panelAll.pack(4,6);

        panelAll.de.margin = 1;
        panelAll.marginbottom = 0;
        panelAll.margin = [1 1 1 1];
jj=0;
for j= [5 4 2 1]
        jj = jj+1;%index counter
    
        
        measure = AllResults.(graphNames{i}).(measureNames{j});
        if(j==1)%if betweenness, then we need to normalize
            measure = measure./((numLab-1)*(numLab-2));
        end
        
        %we have the measure for each subgroup
        GroupMeasures = cell(1,3);
        GroupMeasures{1}=measure(HCIND,:);
        GroupMeasures{2}=measure(TauIND,:);
        GroupMeasures{3}=measure(TDPIND,:);
        
        saveDirSub = fullfile(saveDirBase,graphNames{i},measureNames{j});
        %mkdir(saveDir)
        CRange=([0 prctile(measure(:),95)]);
        

        %We make these comparisons between groups
        testG1All=[3];
        testG2All=[1];
            for g = 1
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
                
                %correct for multiple comparison
                notNaN = ~isnan(TTestsAll_pVal);
                [padj_fdr,alpha_fdr] = multicmp(TTestsAll_pVal(notNaN),'fdr',0.001);
                
                allPadj = {padj_fdr};
                allPadjName = {'FDR_Benjamini-Hochberg',...
                    'FWE_Hochberg-Bonferroni','FWE_Holm-Bonferroni'};
                
                for p = 1%:3
                    data = TTestsAll_pVal;
                    data(~notNaN) = 1;
                    data(notNaN) = allPadj{p};
                    CRange=[0 1];
                    saveName = testName;
                    saveDir = fullfile(saveDirSub,testName,allPadjName{p});
                    %plotBrain2(labelmesh,data,LabelLUT,CRange,saveDir,saveName,0.001,panelAll,plotg2(g))
                end

                data = TTestsAll_tstat;
                CRange=[-4 4];
                saveName = testName;
                saveDir = fullfile(saveDirSub,testName,'tstat');
                %plotBrain2(labelmesh,data,LabelLUT,CRange,saveDir,saveName,0,panelAll,plotg1(jj))
                plotBrain3(Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,plotg1(jj))
            end
                
%            print('-dtiff','-r500',fnames)
%            saveas(H2,fullfile(saveDirBase,'fig1.tif'))
            
 
            
end
end

           print(H2,fullfile(saveDirBase,['fig5_'  testName '.tif']),'-dpng','-r400');