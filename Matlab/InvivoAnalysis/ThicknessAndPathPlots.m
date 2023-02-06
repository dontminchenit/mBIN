%LoadData

dorsAttnLabels = contains(LabelLUT.Label_Name2,'DorsAttn')';
ventAttnLabels = contains(LabelLUT.Label_Name2,'VentAttn')';
PickedLabels=ones(size(ventAttnLabels));
%LN=sum(CurrLabels);
%metrics
i=4;%thickness
j=6;%mean

%create pannels
H2=figure(2)
clf
panelAll = panel();
panelAll.pack(2,6);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

measure = AllResults.Thickness.Mean;
if(j==1)%if betweenness, then we need to normalize
    measure = measure./((numLab-1)*(numLab-2));
end

%we have the measure for each subgroup
GroupMeasures = cell(1,3);
GroupMeasures{1}=measure(CohortsIndexs==1,:);
GroupMeasures{2}=measure(CohortsIndexs==2,:);
GroupMeasures{3}=measure(CohortsIndexs==3,:);

pathmeasure=pathDataGM;
GroupPath = cell(1,3);
GroupPath{1}=pathmeasure(CohortsIndexs==1,:);
GroupPath{2}=pathmeasure(CohortsIndexs==2,:);
GroupPath{3}=pathmeasure(CohortsIndexs==3,:);

%saveDirSub = fullfile(saveDirBase,graphNames{i},measureNames{j});
%mkdir(saveDir)
CRange=([0 prctile(measure(:),95)]);

GroupNames={'HC','Tau','TDP'}

%We make these comparisons between groups
testG1All=[2 3 2 3];
testG2All=[1 1 3 2];
%plotg1=[1 3 5 -1];
%plotg2=[2 4 6 7];
plotg1=[0 0 3 4];
plotg2=[1 2 0 0];
plotg2=[1 2 3 4];%location of plot
for g = 1:2%tau or TDP
    G1 = testG1All(g);
    G2 = testG2All(g);
    testName=[GroupNames{G1} 'vs' GroupNames{G2}];
    TTestsAll_pVal = nan(1,numLab);
    TTestsAll_tstat = nan(1,numLab);
    zscore = nan(1,numLab);
    
    
    %do t-test on all labels
    for t = 1:numLab
        if(PickedLabels(t))
            
            currCtrlMean = mean(GroupMeasures{1}(:,t),'omitnan');
            currCtrlStd = std(GroupMeasures{1}(:,t),'omitnan');
            zscore(t)=mean((GroupMeasures{G1}(:,t)-currCtrlMean)/currCtrlStd,'omitnan');
            
            
            [h1 p ci stats]=ttest2(GroupMeasures{G1}(:,t),GroupMeasures{G2}(:,t),'Tail','left');
            TTestsAll_pVal(t)=p;
            TTestsAll_tstat(t) = stats.tstat;
        end
    end
    
    %correct for multiple comparison
    notNaN = ~isnan(TTestsAll_pVal);
    pThresh=0.001;
    [padj_fdr,alpha_fdr] = multicmp(TTestsAll_pVal(notNaN),'fdr',pThresh);
    
    allPadj = {padj_fdr};
    allPadjName = {'FDR_Benjamini-Hochberg',...
        'FWE_Hochberg-Bonferroni','FWE_Holm-Bonferroni'};
    
    %                for p = 1%:3
    p=1
    thicknessdata = TTestsAll_pVal
    thicknessdata(~notNaN) = 1;
    thicknessdata(notNaN) = allPadj{p};
    %                    CRange=[0 1];
    %                    saveName = testName;
    %                    saveDir = fullfile(saveDirSub,testName,allPadjName{p});
    %plotBrain2(labelmesh,data,LabelLUT,CRange,saveDir,saveName,0.001,panelAll,plotg2(g))
    %plotBrain3(Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0.001,panelAll,plotg2(g))
    %                    if(g~=3 )
    surfDisp=thicknessdata<pThresh;
    surfDisp=TTestsAll_tstat/3;
    surfDisp=zscore/3;
    pointDisp= mean(GroupPath{G1},1,'omitnan');
    pointDisp=pointDisp/max(pointDisp(:));
    plotBrain3Points(CoM,schaefer400x7SurfaceMeshFromGII,pointDisp,0,panelAll,plotg2(g),2,10,surfDisp)
    
    %                    end
    %                end
    
    %data = TTestsAll_tstat;
    %CRange=[-4 4];
    
    %if(g==1 || g==2) %if comparing to controls we only care about negative tstats (patients < HC)
    %    data(data>0) = 0;
    %end
    
    %saveName = testName;
    %saveDir = fullfile(baseSaveDir,testName,'tstat');
    %plotBrain2(labelmesh,data,LabelLUT,CRange,saveDir,saveName,0,panelAll,plotg1(g))
    %plotBrain3(Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,plotg1(g))
    %if(g==3 | g==4)
    %    plotBrain3Points(CoM,schaefer400x7SurfaceMeshFromGII,data,0.001,panelAll,plotg1(g),1,10)
    %end
end

%            print('-dtiff','-r500',fnames)
%            saveas(H2,fullfile(saveDirBase,'fig1.tif'))

print(H2,fullfile(baseSaveDir,'Thickness.tif'),'-dpng','-r400');