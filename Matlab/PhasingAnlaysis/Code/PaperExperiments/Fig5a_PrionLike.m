%LoadDataSarah

%find center of mass, to show the location of nodes
meshlabels = labelmesh.attributes.attribute_array;
CoM=zeros(height(LabelLUT),3);
for i = 1:height(LabelLUT)
    currLabel = LabelLUT{i,1};
    idx = find(meshlabels==currLabel);
    CoM(i,:)=mean(labelmesh.points(idx,:),1);
end

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

%here we figure out what labels are adjacent to each other
adjmtx=zeros(numLab,numLab);

FN=size(Lausanne250SurfaceMeshFromGII.giiSurface_Both.faces,1);
vertexLabels = Lausanne250SurfaceMeshFromGII.cdata;
for f=1:FN
    
    currFace = Lausanne250SurfaceMeshFromGII.giiSurface_Both.faces(f,:);
    
    %convertFace vertex to unique labels
    FaceLabels = unique(vertexLabels(currFace));%unique face Labels
    
    if(length(FaceLabels)>1)
        
        for i = 1:length(FaceLabels)
            for j = 1:length(FaceLabels)
                if(i~=j && FaceLabels(i)~=0 && FaceLabels(j)~=0)
                    
                    labIndx1=find(LabelLUT.Label_ID==FaceLabels(i));
                    labIndx2=find(LabelLUT.Label_ID==FaceLabels(j));
                    
                    adjmtx(labIndx1,labIndx2)=adjmtx(labIndx1,labIndx2)+1;
                end
            end
        end
    end
    
end
figure(1)
clf
imshow(adjmtx)
colormap(jet)
caxis([0 max(adjmtx(:))])
colorbar
%WeakCount

adjRegcount = zeros(2,5,3);

HubWeakPercent = nan(NS,2);

for m=1%[1 2 4 5]
    
    prevKeepTau = ones(length(T_tauStagingLaus),1);
    prevKeepTDP= ones(length(T_tauStagingLaus),1);
    for p = 1:2
        
        for s = 1:NS
            for g = 2
                p
                G1=compPairs(1,p);
                G2=compPairs(2,p);
                testName=[GroupNames{G1} 'vs' GroupNames{G2}];
                
                
                %Pick Stages
                stageA=s;
                %stageA=1;
                stageB=s+1;
                
                
                
                %             if(G1==2)
                %                  rmIndx=(T_tauStagingLaus==stageA & prevKeepTau) | T_tauStagingLaus==stageB;
                %                  rmIndx1=(T_tauStagingLaus==stageA & prevKeepTau);
                %                  rmIndx2=T_tauStagingLaus==stageB;
                %             elseif(G1==3)
                %                  rmIndx=(T_tdpStagingLaus==stageA & prevKeepTDP) | T_tdpStagingLaus==stageB;
                %                  rmIndx1=(T_tdpStagingLaus==stageA & prevKeepTDP);
                %                  rmIndx2=T_tdpStagingLaus==stageB;
                %             end
                %
                if(G1==2)
                    rmIndx=(T_tauStagingLaus==stageA) | T_tauStagingLaus==stageB;
                    rmIndx1=(T_tauStagingLaus==stageA);
                    rmIndx2=T_tauStagingLaus==stageB;
                elseif(G1==3)
                    rmIndx=(T_tdpStagingLaus==stageA) | T_tdpStagingLaus==stageB;
                    rmIndx1=(T_tdpStagingLaus==stageA);
                    rmIndx2=T_tdpStagingLaus==stageB;
                end
                
                
                keepIndx=rmIndx;
                adjmtxCopy=adjmtx;
                adjmtxCopy(~rmIndx,:)=0;
                adjmtxCopy(:,~rmIndx)=0;
                adjmtxCopy(rmIndx1,rmIndx1)=0;
                adjmtxCopy(rmIndx2,rmIndx2)=0;
                
                
                if(G1==2)
                    prevKeepTau = sum(adjmtxCopy,2) > 0;
                elseif(G1==3)
                    prevKeepTDP = sum(adjmtxCopy,2) > 0;
                end
                
                
                adjRegcount(p,s,1)=sum(adjmtxCopy(:)>0)/2;
                
                keepmtx = ones(size(adjmtx));
                keepmtx(~rmIndx,:)=0;
                keepmtx(:,~rmIndx)=0;
                keepmtx(rmIndx1,rmIndx1)=0;
                keepmtx(rmIndx2,rmIndx2)=0;
                
                adjRegcount(p,s,2)=sum(rmIndx1>0);
                adjRegcount(p,s,3)=sum(rmIndx2>0);
                
                %ColorVec=rmIndx1*s+rmIndx2*(s+1);
                %ColorVec=rmIndx1+rmIndx2*(s+1);
                
                ColorMode = 1;
                
                switch ColorMode
                    
                    case 1
                        ColorVec=rmIndx1*s+rmIndx2*(s+1);
                        %             ColorVec=rmIndx1+rmIndx2*(s+1);
                        
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
                        
                        [f1 f2]=find(triu(adjmtxCopy));
                        
                        hits = Hubs1(f1) | Hubs1(f2) | Hubs2(f1) | Hubs2(f2);
                        
                        HubWeakPercent(s,p)=sum(hits)/length(hits)
                        
                        HubRelatedMatrix = zeros(length(Hubs1),length(Hubs1));
                        HubRelatedMatrix(Hubs1 | Hubs2,:) = 1;
                        HubRelatedMatrix(:,Hubs1 | Hubs2) = 1;
                        InHub=logical(triu(HubRelatedMatrix));
                        
                        nonHubRelatedMatrix = zeros(length(Hubs1),length(Hubs1));
                        nonHubRelatedMatrix(nonHubs1 | nonHubs2,:) = 1;
                        nonHubRelatedMatrix(:,nonHubs1 | nonHubs2) = 1;
                        OutHub=logical(triu(nonHubRelatedMatrix));
                        
                        InHubdiffVals = diffVals(InHub & valid);
                        outHubdiffVals = diffVals(OutHub & valid);
                        
                        
                        inoutHubData{s,p,1}=InHubdiffVals;
                        inoutHubData{s,p,2}=outHubdiffVals;
                        
                        [h pval]=ttest2(InHubdiffVals,outHubdiffVals);
                        
                        inoutHubPtest(s,p) = pval;
                        
                        
                        
                        
                        
                        
                end
                
                
                plotNetwork3(adjmtxCopy, Lausanne250SurfaceMeshFromGII, CoM,[],[0 max(adjmtx(:))], .75*keepIndx+.001,panelAll,s,ColorVec,2)
                
            end
        end
        
        %    print(H2,fullfile(saveDirBase,['fig8_' GroupNames{G1} '_AdjStageConnectionsEpiRegion_' hubName '_' num2str(p)   '.tif']),'-dpng','-r400');
        print(H2,fullfile(saveDirBase,['fig8_' GroupNames{G1} '_AdjStageTransition' num2str(p)   '.tif']),'-dpng','-r400');
    end
    
    % StagingPlot(DistanceComp,1,2,saveDirBase,'AllDist','Distance')
    % StagingPlot(CountComp,2,2,saveDirBase,'AllCount','Count')
    % StagingPlot(StageDistanceComp,1,2,saveDirBase,'StageDist','Distance')
    % StagingPlot(StageCountComp,2,2,saveDirBase,'StageCount','Count')
    % StagingPlot(WeakDistanceComp,1,2,saveDirBase,'EpiWeakAllDist','Distance')
    % StagingPlot(WeakCountComp,2,2,saveDirBase,'EpiWeakAllCount','Count')
    %Fig6PlotMetric
end