%LoadDataSarah
             
%this sets which mode of analysis we're doing
%Set to 1 for Figures 4a & 4b
%Set to 2 for Figures 4c & 4d (and then run Fig4cd_EdgeConnectedToHubs.m)
ModeOrig1_HubAnalysis2 = 2;

%this flag determines if we're looking at adjacent or
%epicenter relationships (Figure 4a/c vs 4b/d)
flagEpi0Adj1=1;


%Initialize
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

indvWeakEdgeCountNormalized= nan(200,NS,2);
indvWeakEdgeCount= nan(200,NS,2);

indvTauHubsEdgeCountNormalized= nan(200,NS,2);
indvTDPHubsEdgeCountNormalized= nan(200,NS,2);
indvTauHubsEdgeCount= nan(200,NS,2);
indvTDPHubsEdgeCount= nan(200,NS,2);
PrctHubsEdgeCount= nan(200,NS,2);
TotalHubsEdgeCount= nan(200,NS,2);

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
                fullmtx1=zeros(size(allmtx1{1},1),size(allmtx1{1},2),length(allmtx1));%This stores all the matrices for Group 1 (disease)
                for k = 1:length(allmtx1)
                    currMtx=allmtx1{k};
                    currMtx(currMtx==0)=nan;%set zero FA as NaN, since no measurements
                    fullmtx1(:,:,k) = currMtx;
                end
                                              
              
                %find the mean FA for Group 2
                allmtx2=AllResults.(graphNames{g}).('WeightMtx')(GroupInd{G2});
                fullmtx2=zeros(size(allmtx2{1},1),size(allmtx2{1},2),length(allmtx2));%This stores all the matrices for Group 2 (Control)
                for k = 1:length(allmtx2)
                    currMtx = allmtx2{k};
                    currMtx(currMtx==0)=nan;%set zero FA as NaN, since no measurements
                    fullmtx2(:,:,k)=currMtx;
                end
                
                %Calculate mean of both groups, and std of control group
                %(for each edge)
                meanmtx1 = mean(fullmtx1,3,'omitnan');
                meanmtx2 = mean(fullmtx2,3,'omitnan');
                stdmtx2 = std(fullmtx2,0,3,'omitnan');
                stdmtx2(stdmtx2==0)=nan;%if we have a stddev of zero, it probably means only a single data point, so exclude

                
                %find the zscore of disease mean relative to controls
                meanZscore=(meanmtx1-meanmtx2)./stdmtx2;
              

                
                %meanZscore=(meanmtx1-meanmtx2);
                
                valid = (meanmtx2~=0 & meanmtx1~=0 & ~isnan(meanmtx2) & ~isnan(meanmtx1) & ~isnan(stdmtx2));%valid edges
                meanZscore(~valid) = nan;
                
                %set the thresholds as below 3 standard deviation (-3
                %z-score)
                thresh=-3;
                %thresh=mean(meanZscore(valid)) - 3*std(meanZscore(valid));
                
                %keep valid connections below 3 stddev
                keepmtx = (meanZscore<=thresh) & valid;
                adjmtx1=meanmtx1;
                adjmtx1(~keepmtx)=nan;

                %set the phases
                if(flagEpi0Adj1)
                    stageA=s;
                    Epi0Adj1Name='Adj';
                else
                    stageA=1;
                    Epi0Adj1Name='Epi';

                end
                stageB=s+1;
                
                
                %depending on which cohort, we now pick stage regions.
                if(G1==2)
                    keepIndxBoth=T_tauStagingLaus==stageA | T_tauStagingLaus==stageB;
                    keepIndxA=T_tauStagingLaus==stageA;
                    keepIndxB=T_tauStagingLaus==stageB;
                elseif(G1==3)
                    keepIndxBoth=T_tdpStagingLaus==stageA | T_tdpStagingLaus==stageB;
                    keepIndxA=T_tdpStagingLaus==stageA;
                    keepIndxB=T_tdpStagingLaus==stageB;
                end

                %this matrix determines which edges we keep for current two
                %phases
                keepmtx = ones(size(adjmtx1)); 
                keepmtx(~keepIndxBoth,:)=0;
                keepmtx(:,~keepIndxBoth)=0;
                keepmtx(keepIndxA,keepIndxA)=0;
                keepmtx(keepIndxB,keepIndxB)=0;
                
                        
                %we remove the edges from the disease matrix not in phase
                adjmtx1(~keepmtx)=nan;
                
                %find the zscore of all subjects
                indvZscore = nan(size(fullmtx1));
                
                for k = 1:length(allmtx1)
                    indvZscore(:,:,k) = (fullmtx1(:,:,k)-meanmtx2)./stdmtx2;
                    
                    WeakEdge=double(indvZscore(:,:,k)<-3);
                    WeakEdge(WeakEdge==0)=nan;
                    WeakEdge(isnan(stdmtx2))=nan;
                    WeakEdge(~keepmtx)=nan;
                    WeakEdge(~valid)=nan;
                    
                    PossibleEdge=double(fullmtx1(:,:,k));
                    PossibleEdge(PossibleEdge==0)=nan;
                    %PossibleEdge(isnan(stdmtx2))=nan;
                    PossibleEdge(~keepmtx)=nan;
                    %PossibleEdge(~valid)=nan;
                    
                    indvWeakEdgeCountNormalized(k,s,p)=sum(~isnan(WeakEdge(:)))/(sum(keepIndxA)*sum(keepIndxB));
                    indvWeakEdgeCount(k,s,p)=sum(~isnan(WeakEdge(:)));
                end
                            
                
                
                
                %find all the locations where there are weak edges
                [ii,jj] = find(adjmtx1);
                weaknodes = unique([ii, jj]);
                
   
                switch ModeOrig1_HubAnalysis2
                    
                    case 1
                        if(flagEpi0Adj1)
                            ColorVec=keepIndxA*s+keepIndxB*(s+1);
                        else
                            ColorVec=keepIndxA+keepIndxB*(s+1);
                        end
                        ModeName = 'StageTransition_';
                        hubName='';
                        
                    case 2
                        ColorVec=keepIndxA+keepIndxB;
                        ModeName = 'StageHubAndEdges_';
                        switch 2
                            case 1
                                hubName='DiseaseHubs';
                                Hubs1 = keepIndxA & Hubs.FAGraph.HCThresh.(GroupNames{G1})';
                                Hubs2 = keepIndxB & Hubs.FAGraph.HCThresh.(GroupNames{G1})';
                                
                                
                            case 2
                                hubName='CtrlHubs';
                                Hubs1 = keepIndxA & Hubs.FAGraph.HCThresh.(GroupNames{G2})';
                                Hubs2 = keepIndxB & Hubs.FAGraph.HCThresh.(GroupNames{G2})';
                                
                            case 3
                                
                                hubName='LostHubs';
                                Hubs1 = keepIndxA & ((Hubs.FAGraph.HCThresh.(GroupNames{G1})' - Hubs.FAGraph.HCThresh.(GroupNames{G2})')<0);
                                Hubs2 = keepIndxB & ((Hubs.FAGraph.HCThresh.(GroupNames{G1})' - Hubs.FAGraph.HCThresh.(GroupNames{G2})')<0);
                        end
                        
                        %index of places that non hub in each group
                        nonHubs1 = keepIndxA & ~Hubs1;
                        nonHubs2 = keepIndxB & ~Hubs2;
                                                
                        ColorVec=ColorVec+2*Hubs1+Hubs2;
                        keepIndxBoth=keepIndxBoth+0.5*Hubs1+0.5*Hubs2;
                        
                        
                        %find weakened edges
                        [f1 f2]=find(triu(~isnan(adjmtx1)));
                        hits = Hubs1(f1) | Hubs1(f2) | Hubs2(f1) | Hubs2(f2);
                        
                        HubWeakPercent(s,p)=sum(hits)/length(hits);
                        

                        %using the meanZscores, we look at just where the
                        %weakened edges line up with the hubs/nonHubs
                        HubRelatedMatrix = zeros(size(meanZscore));
                        HubRelatedMatrix(Hubs1 | Hubs2,:) = 1;
                        HubRelatedMatrix(:,Hubs1 | Hubs2) = 1;
                        InHub=logical(triu(HubRelatedMatrix));
                        
                        nonHubRelatedMatrix = zeros(size(meanZscore));
                        nonHubRelatedMatrix(nonHubs1 | nonHubs2,:) = 1;
                        nonHubRelatedMatrix(:,nonHubs1 | nonHubs2) = 1;
                        nonHubRelatedMatrix(Hubs1 | Hubs2,:) = 0;
                        nonHubRelatedMatrix(:,Hubs1 | Hubs2) = 0;
                        OutHub=logical(triu(nonHubRelatedMatrix));
                        
                        InHubZscores = meanZscore(InHub & valid);
                        OutHubZscores = meanZscore(OutHub & valid);
                        

                        %save zscores associated with hubs/nonhubs
                        inoutHubData{s,p,1}=InHubZscores;
                        inoutHubData{s,p,2}=OutHubZscores;
                        
                        %now do it for individual level
                        
                        for kk= 1:size(indvZscore,3)
                            
                            curZscore = indvZscore(:,:,kk)<=-3;
                            currInHub = curZscore(InHub & valid);                            
                            currOutHub = curZscore(OutHub & valid);
                            
                            possibleIn = InHub & valid;
                            possibleOut = OutHub & valid;
                            
                            %countInhub = mean(currInHub(:),'omitnan');
                            %countOutHub = mean(currOutHub(:),'omitnan');
                        
                            countInhub = sum(currInHub(:))/(sum(currInHub(:))+sum(currOutHub(:)));
                            countOutHub = sum(currOutHub(:))/(sum(currInHub(:))+sum(currOutHub(:)));
                            
                            
                            PrctHubsEdgeCount(kk,s,p)= sum(currInHub(:))/(sum(currInHub(:))+sum(currOutHub(:)));
                            TotalHubsEdgeCount(kk,s,p)= sum(currInHub(:));
                        if (p==1)
                            indvTauHubsEdgeCount(kk,s,1)=countInhub;
                            indvTauHubsEdgeCount(kk,s,2)=countOutHub;
                            
                            indvTauHubsEdgeCountNormalized(kk,s,1)=countInhub/(sum(keepIndxA)*sum(keepIndxB));
                            indvTauHubsEdgeCountNormalized(kk,s,2)=countOutHub/(sum(keepIndxA)*sum(keepIndxB));
                        elseif (p==2)
                            indvTDPHubsEdgeCount(kk,s,1)=countInhub;
                            indvTDPHubsEdgeCount(kk,s,2)=countOutHub;
                            
                            indvTDPHubsEdgeCountNormalized(kk,s,1)=countInhub/(sum(keepIndxA)*sum(keepIndxB));
                            indvTDPHubsEdgeCountNormalized(kk,s,2)=countOutHub/(sum(keepIndxA)*sum(keepIndxB));
                        end
                        
                        end
                        
                        %do a ttest between the groups
                        [h pval]=ttest2(InHubZscores,OutHubZscores);
                        inoutHubPtest(s,p) = pval;
                        
                        %save the mean difference between the groups
                        inoutHubMeanDiff(s,p)=mean(InHubZscores,'omitnan') - mean(OutHubZscores,'omitnan');
                end
                
                %if not empty, use the max
                if(sum(~isnan(adjmtx1))>0)
                cRange=[0 max(adjmtx1(:))];
                else
                cRange=[0 1];
                end
                plotNetwork3(adjmtx1, Lausanne250SurfaceMeshFromGII, CoM,[],cRange,keepIndxBoth+.001,panelAll,s,ColorVec,0)
                
            end
        end
                
                
            print(H2,fullfile(saveDirBase,['fig5_' GroupNames{G1} '_' Epi0Adj1Name '_' ModeName '_' hubName '_' num2str(p)   '.tif']),'-dpng','-r400');
           % PlotComparisons
          %  print(H3,fullfile(saveDirBase,['fig8_' GroupNames{G1} '_' Epi0Adj1Name '_' ModeName '_' hubName '_' num2str(p)   '.tif']),'-dpng','-r400');
    end
    
    % StagingPlot(DistanceComp,1,2,saveDirBase,'AllDist','Distance')
    % StagingPlot(CountComp,2,2,saveDirBase,'AllCount','Count')
    % StagingPlot(StageDistanceComp,1,2,saveDirBase,'StageDist','Distance')
    % StagingPlot(StageCountComp,2,2,saveDirBase,'StageCount','Count')
    % StagingPlot(WeakDistanceComp,1,2,saveDirBase,'EpiWeakAllDist','Distance')
    % StagingPlot(WeakCountComp,2,2,saveDirBase,'EpiWeakAllCount','Count')
    %Fig6PlotMetric
end