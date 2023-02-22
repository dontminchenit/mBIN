for i = 1:2
    
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
        for g_thresh=1:3
            a = betweennessMean.(GroupNames{g1}) > HubsBetweennessThresh.(GroupNames{g_thresh});
            b = degreeMean.(GroupNames{g1}) > HubsDegThresh.(GroupNames{g_thresh});
            Hubs.(graphNames{i}).([(GroupNames{g_thresh}) 'Thresh']).(GroupNames{g1}) = a | b;
            
            a = betweenness > HubsBetweennessThresh.(GroupNames{g_thresh});
            b = degree > HubsDegThresh.(GroupNames{g_thresh});
            IndivHubs.(graphNames{i}).([(GroupNames{g_thresh}) 'Thresh']).(GroupNames{g1}) = a | b;

            %if (g_thresh==1)
                %plot Hubs
            %    data = Hubs.(graphNames{i}).([(GroupNames{g_thresh}) 'Thresh']).(GroupNames{g1});
            %    saveName = [GroupNames{g1} '_Hubs'];
            %    saveDir = fullfile(saveDirBase,graphNames{i},'Hubs',[(GroupNames{g_thresh}) 'Thresh'],GroupNames{g1});
            %    CRange = [0 1];
                %plotBrain(labelmesh,data,LabelLUT,CRange,saveDir,saveName,0)
            %end
        end
    end
    
end