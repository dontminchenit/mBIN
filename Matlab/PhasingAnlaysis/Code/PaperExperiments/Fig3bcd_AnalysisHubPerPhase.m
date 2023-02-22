%for m = 1:5
%cross reference with metrics 
allData = nan(6,3);

mName = 'Number of Hubs'

for s=1:6
for j = 2:3
    
    if j == 2
        StagingSel=T_tauStagingLaus;
    else 
        StagingSel=T_tdpStagingLaus;
    end
    
    
    vals = Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{j})' & StagingSel==s;
    allData(s,j,1) = sum(vals)/sum(StagingSel==s);%each subject/region value is a datapoint
end
end


H=figure(1)
%bar(allData(:,2:3))
bar(allData(:,2:3)*100)
legend({'FTLD-tau', 'FTLD-TDP'})
xlabel('Disease Phases')
ylabel('% Nodes Are Disease Hubs')
ytickformat('percentage')
set(gca,'FontSize',20)
print(H,fullfile(saveDirBase,['fig3_' mName '.tif']),'-dpng','-r400');



for s=1:6
for j = 2:3
    
    if j == 2
        StagingSel=T_tauStagingLaus;
    else 
        StagingSel=T_tdpStagingLaus;
    end
    
    
    vals = Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{1})' & StagingSel==s;
    sum(vals)
    sum(StagingSel==s)
    allData(s,j,1) = sum(vals)/sum(StagingSel==s);%each subject/region value is a datapoint
end
end


H=figure(2)
%bar(allData(:,2:3))
bar(allData(:,2:3)*100)
legend({'FTLD-tau', 'FTLD-TDP'})
xlabel('Disease Phases')
ylabel('% Nodes Are HC Hubs')
ytickformat('percentage')
set(gca,'FontSize',20)
print(H,fullfile(saveDirBase,['fig3_HC_' mName '.tif']),'-dpng','-r400');


%for m = 1:5
%cross reference with metrics 
allData = nan(6,3);

mName = 'Number of Hubs'

for s=1:6
for j = 2:3
    
    if j == 2
        StagingSel=T_tauStagingLaus;
    else 
        StagingSel=T_tdpStagingLaus;
    end
    
    
    diffHub = Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{j})' - Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{1})';
    vals = diffHub<0 & StagingSel==s;
    allData(s,j,1) = sum(vals)/sum(Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{1})' & StagingSel==s);%each subject/region value is a datapoint
end
end


H=figure(3)
bar(allData(:,2:3)*100)
%legend({'FTLD-tau', 'FTLD-TDP'})
xlabel('Disease Phases')
ylabel('% HC Hubs Lost in Disease')
ytickformat('percentage')
set(gca,'FontSize',20)
print(H,fullfile(saveDirBase,['fig3_Lostprct_' mName '.tif']),'-dpng','-r400');



for s=1:6
for j = 2:3
    
    if j == 2
        StagingSel=T_tauStagingLaus;
    else 
        StagingSel=T_tdpStagingLaus;
    end
    
    
        
    diffHub = Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{j})' - Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{1})';
    vals = diffHub>0 & StagingSel==s;
    %sum(vals)
    %sum(Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{1})' & StagingSel==s)
    allData(s,j,1) = sum(vals)/sum(Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{1})' & StagingSel==s);%each subject/region value is a datapoint
end
end


H=figure(4)
bar(allData(:,2:3)*100)
%legend({'FTLD-tau', 'FTLD-TDP'}, 'Location', 'nw')
xlabel('Disease Phases')
ylabel('% New Hubs in Disease')
ytickformat('percentage')
set(gca,'FontSize',20)
print(H,fullfile(saveDirBase,['fig3_Gainedprct_' mName '.tif']),'-dpng','-r400');



