LoadDataSarah


NS=5;
H2=figure(2)
clf
panelAll = panel();
panelAll.pack(2,NS);

panelAll.de.margin = 1;
panelAll.marginbottom = 0;
panelAll.margin = [1 1 1 1];

compPairs = [2 3 2 3;...
    1 1 3 2]

FAByStage = nan(NS,2,100000);
MetricBySTage = nan(4,2,100000);

DistanceComp = nan(4,2,100000);
CountComp =  nan(4,2,2);
TotalCount =  nan(4,2,2);

StageDistanceComp = nan(5,2,100000);
StageCountComp =  nan(5,2,2);


WeakDistanceComp = nan(4,2,100000);
WeakCountComp =  nan(4,2,2);

HubWeakPercent = nan(NS,2);

inoutHubData = cell(NS,2,2);
inoutHubPtest = nan(NS,2);
inoutHubMeanDiff = nan(NS,2);

%WeakCount
for m=1%[1 2 4 5]
    for p = 1:2
        
        for s = 1:NS
            for g = 2
                p
                %Define Groups
                G1=compPairs(1,p);
                G2=compPairs(2,p);
                testName=[GroupNames{G1} 'vs' GroupNames{G2}];
                
                %find the mean FA for Group 1
                allmtx1=AllResults.(graphNames{g}).('WeightMtx')(GroupInd{G1});
                meanmtx1 = allmtx1{1};
                for k = 2:length(allmtx1)
                    meanmtx1 = meanmtx1+allmtx1{k};
                end
                meanmtx1=meanmtx1./length(allmtx1);
                caseEdgeVAls = meanmtx1;
                
                CaseMetric = AllResults.FAGraph.(measureNames{m})(GroupInd{G1},:);
                
                %find the mean FA for Group 2
                allmtx2=AllResults.(graphNames{g}).('WeightMtx')(GroupInd{G2});
                meanmtx2 = allmtx2{1};
                fullmtx2=zeros(size(allmtx2{1},1),size(allmtx2{1},2),length(allmtx2));
                fullmtx2(:,:,1)=allmtx2{1};
                for k = 2:length(allmtx2)
                    fullmtx2(:,:,k) = allmtx2{k};
                    meanmtx2 = meanmtx2+allmtx2{k};
                end
                meanmtx2=meanmtx2./length(allmtx2);
                ctrlEdgeVAls = meanmtx2;
                
                CtrlMetric = AllResults.FAGraph.(measureNames{m})(GroupInd{G2},:);
                
                %find the difference between the two group
                diffVals=meanmtx1-meanmtx2;%difference between group 1 and 2
                meanZscore = (meanmtx1-meanmtx2)./std(fullmtx2,0,3);
                
                valid = (meanmtx2~=0 & meanmtx1~=0);
                
                %set the 3 standard deviation threshold
                thresh=mean(diffVals(valid))-3*std(diffVals(valid)); %find bottom 3 stddev of groupwise differences
                
                meanmtx1(diffVals>thresh | diffVals>0 )=0;%remove all connections above 3 stddev
                
                saveDir = fullfile(saveDirBase,graphNames{g},'NetworkCmp',testName);
                mkdir(saveDir);
                
                val = CoMDistMatrix(meanmtx1>0) ;
                DistanceComp(s,p,1:length(val)) = val;
                CountComp(s,p,1) =  sum(meanmtx1(:)>0);
                
                %keepIndx = ones(size(CoM(:,1)));
                
                flagEpi0Adj1=0;
                
                if(flagEpi0Adj1)
                    stageA=s;
                else
                    stageA=1;
                end
                stageB=s+1;
                
                
                if(G1==2)
                    rmIndx=T_tauStagingLaus==stageA | T_tauStagingLaus==stageB;
                    rmIndx1=T_tauStagingLaus==stageA;
                    rmIndx2=T_tauStagingLaus==stageB;
                elseif(G1==3)
                    rmIndx=T_tdpStagingLaus==stageA | T_tdpStagingLaus==stageB;
                    rmIndx1=T_tdpStagingLaus==stageA;
                    rmIndx2=T_tdpStagingLaus==stageB;
                end
                keepIndx=rmIndx;
                %meanmtx1(rmIndx,~rmIndx)=0;
                %meanmtx1(~rmIndx,rmIndx)=0;
                keepmtx = ones(size(meanmtx1));
                keepmtx(~rmIndx,:)=0;
                keepmtx(:,~rmIndx)=0;
                keepmtx(rmIndx1,rmIndx1)=0;
                keepmtx(rmIndx2,rmIndx2)=0;
                
                
                stage1mtx = meanmtx1;
                stage1mtx(~rmIndx1,~rmIndx1)=0;
                stage2mtx = meanmtx1;
                stage2mtx(~rmIndx2,~rmIndx2)=0;
                
                val = CoMDistMatrix( stage1mtx>0) ;
                StageDistanceComp(stageA,p,1:length(val)) = val;
                StageCountComp(stageA,p,1) =  sum(stage1mtx(:)>0);
                
                
                val = CoMDistMatrix( stage2mtx>0) ;
                StageDistanceComp(stageB,p,1:length(val)) = val;
                StageCountComp(stageB,p,1) =  sum(stage2mtx(:)>0);
                
                
                meanmtx1(~keepmtx)=0;
                val = CoMDistMatrix(meanmtx1>0) ;
                WeakDistanceComp(s,p,1:length(val)) = val;
                WeakCountComp(s,p,1) =  sum(meanmtx1(:)>0);
                
                
                [ii,jj] = find(meanmtx1);
                
                weaknodes = unique([ii, jj]);
                
                stagecaseEdgeVAls=caseEdgeVAls;
                stagecaseEdgeVAls(~keepmtx)=0;
                
                stagectrlEdgeVAls=ctrlEdgeVAls;
                stagectrlEdgeVAls(~keepmtx)=0;
                
                val=stagecaseEdgeVAls(stagecaseEdgeVAls>0);
                val=stagecaseEdgeVAls(stagecaseEdgeVAls>0)./stagectrlEdgeVAls(stagecaseEdgeVAls>0);
                
                FAByStage(s,p,1:length(val)) = val;
                
                val = mean(CaseMetric(:,weaknodes),2);
                MetricBySTage(s,p,1:length(val(:))) = val(:);
                
                
                ColorMode = 2;
                
                switch ColorMode
                    
                    case 1
                        if(flagEpi0Adj1)
                            ColorVec=rmIndx1*s+rmIndx2*(s+1);
                        else
                            ColorVec=rmIndx1+rmIndx2*(s+1);
                        end
                        hubName='';
                        
                    case 2
                        ColorVec=rmIndx1+rmIndx2;
                        
                        switch 1
                            case 1
                                hubName='DiseaseHubs';
                                Hubs1 = rmIndx1 & Hubs.FAGraph.HCThresh.(GroupNames{G1})';
                                Hubs2 = rmIndx2 & Hubs.FAGraph.HCThresh.(GroupNames{G1})';
                                
                                
                            case 2
                                hubName='CtrlHubs';
                                Hubs1 = rmIndx1 & Hubs.FAGraph.HCThresh.(GroupNames{G2})';
                                Hubs2 = rmIndx2 & Hubs.FAGraph.HCThresh.(GroupNames{G2})';
                                
                            case 3
                                
                                hubName='LostHubs';
                                Hubs1 = rmIndx1 & ((Hubs.FAGraph.HCThresh.(GroupNames{G1})' - Hubs.FAGraph.HCThresh.(GroupNames{G2})')<0);
                                Hubs2 = rmIndx2 & ((Hubs.FAGraph.HCThresh.(GroupNames{G1})' - Hubs.FAGraph.HCThresh.(GroupNames{G2})')<0);
                        end
                        
                        nonHubs1 = rmIndx1 & ~Hubs1;
                        nonHubs2 = rmIndx2 & ~Hubs2;
                        
                        
                        ColorVec=ColorVec+Hubs1.*3+Hubs2;
                        
                        [f1 f2]=find(triu(meanmtx1));
                        
                        
                        hits = Hubs1(f1) | Hubs1(f2) | Hubs2(f1) | Hubs2(f2);
                        
                        HubWeakPercent(s,p)=sum(hits)/length(hits)
                        
                        diffVals = meanZscore;
                        
                        HubRelatedMatrix = zeros(size(diffVals));
                        HubRelatedMatrix(Hubs1 | Hubs2,:) = 1;
                        HubRelatedMatrix(:,Hubs1 | Hubs2) = 1;
                        InHub=logical(triu(HubRelatedMatrix));
                        
                        nonHubRelatedMatrix = zeros(size(diffVals));
                        nonHubRelatedMatrix(nonHubs1 | nonHubs2,:) = 1;
                        nonHubRelatedMatrix(:,nonHubs1 | nonHubs2) = 1;
                        OutHub=logical(triu(nonHubRelatedMatrix));
                        
                        InHubdiffVals = diffVals(InHub & valid);
                        outHubdiffVals = diffVals(OutHub & valid);
                        

                        inoutHubData{s,p,1}=InHubdiffVals;
                        inoutHubData{s,p,2}=outHubdiffVals;
                        
                        [h pval]=ttest2(InHubdiffVals,outHubdiffVals);
                        
                        inoutHubMeanDiff(s,p)=mean(InHubdiffVals) - mean(outHubdiffVals);
                        
                        inoutHubPtest(s,p) = pval;
                        
                end
                
                
                plotNetwork3(meanmtx1, Lausanne250SurfaceMeshFromGII, CoM,saveDir,testName,[],[0 max(meanmtx1(:))],1, keepIndx+.001,panelAll,s,ColorVec)
                
            end
        end
        if(flagEpi0Adj1)

            PlotComparisons
            title([GroupNames{G1} ' Using ' hubName])
            print(H2,fullfile(saveDirBase,['fig6_' GroupNames{G1} '_AdjStageHubAndEdges_' hubName '_' num2str(p)   '.tif']),'-dpng','-r400');

            %print(H2,fullfile(saveDirBase,['fig6_' GroupNames{G1} '_AdjStageTransition_' hubName '_' num2str(p)   '.tif']),'-dpng','-r400');
        else
            PlotComparisons
            title([GroupNames{G1} ' Using ' hubName])
            print(H2,fullfile(saveDirBase,['fig6_' GroupNames{G1} '_EpicenterHubAndEdges_' hubName '_' num2str(p)   '.tif']),'-dpng','-r400');

            %print(H2,fullfile(saveDirBase,['fig6_' GroupNames{G1} '_EpicenterStageTransition_' hubName '_' num2str(p)   '.tif']),'-dpng','-r400');
        end
    end
    
    % StagingPlot(DistanceComp,1,2,saveDirBase,'AllDist','Distance')
    % StagingPlot(CountComp,2,2,saveDirBase,'AllCount','Count')
    % StagingPlot(StageDistanceComp,1,2,saveDirBase,'StageDist','Distance')
    % StagingPlot(StageCountComp,2,2,saveDirBase,'StageCount','Count')
    % StagingPlot(WeakDistanceComp,1,2,saveDirBase,'EpiWeakAllDist','Distance')
    % StagingPlot(WeakCountComp,2,2,saveDirBase,'EpiWeakAllCount','Count')
    %Fig6PlotMetric
end