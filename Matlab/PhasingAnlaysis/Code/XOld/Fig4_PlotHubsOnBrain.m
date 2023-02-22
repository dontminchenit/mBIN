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
    
    
    %findHubs
    betweenness = (AllResults.(graphNames{i}).(measureNames{1}))./((numLab-1)*(numLab-2));
    degree =  AllResults.(graphNames{i}).(measureNames{5});
    
    %We calculate thresholds
    for g=1:3
        betweennessMean.(GroupNames{g}) = nanmean(betweenness(GroupInd{g},:),1);
        degreeMean.(GroupNames{g}) = nanmean(degree(GroupInd{g},:),1);

        a=betweennessMean.(GroupNames{g});
        HubsBetweennessThresh.(GroupNames{g}) = nanmean(a) + std(a,'omitnan');
        a=degreeMean.(GroupNames{g});
        HubsDegThresh.(GroupNames{g}) = nanmean(a) + std(a,'omitnan');
    end

    % Find hubs
    for g1=1:3
        for g2=1:3
            a = betweennessMean.(GroupNames{g1}) > HubsBetweennessThresh.(GroupNames{g2});
            b = degreeMean.(GroupNames{g1}) > HubsDegThresh.(GroupNames{g2});
            Hubs.([(GroupNames{g2}) 'Thresh']).(GroupNames{g1}) = a | b;

            if (g1==1 && g2==1)
            %plot Hubs
            data = Hubs.([(GroupNames{g2}) 'Thresh']).(GroupNames{g1});
            saveName = [GroupNames{g1} '_Hubs'];
            saveDir = fullfile(saveDirBase,graphNames{i},'Hubs',[(GroupNames{g2}) 'Thresh'],GroupNames{g1});
            CRange = [0 1];
            
            plotBrain3(Lausanne250SurfaceMeshFromGII,~data,CRange,saveDir,saveName,0.5,panelAll,plotg1(1))

            end
        end
    end
    
    %Find Hub Differences
    compGroup1 = [2 3 2];
    compGroup2 = [1 1 3];
    for jj = 1:3
 
        g1 = compGroup1(jj);
        g2 = compGroup2(jj);
        
            Hub1 = Hubs.([(GroupNames{1}) 'Thresh']).(GroupNames{g1});
            Hub2 = Hubs.([(GroupNames{1}) 'Thresh']).(GroupNames{g2});

            %plot Hub differences
            data = Hub1 - Hub2;
            saveName = [GroupNames{g1} '_wo_' GroupNames{g2}];
            saveDir = fullfile(saveDirBase,graphNames{i},'Hubs',[(GroupNames{1}) 'Thresh'],[GroupNames{g1} 'without' GroupNames{g2}]);
            CRange = [-1 1];
            if jj == 3
                plotBrain3(Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,plotg1((jj)+1),3) %force a new colomap for the last comparison
            else
                plotBrain3(Lausanne250SurfaceMeshFromGII,data,CRange,saveDir,saveName,0,panelAll,plotg1((jj)+1))
            end

    end
    
           print(H2,fullfile(saveDirBase,'fig4.tif'),'-dpng','-r400'); 

end