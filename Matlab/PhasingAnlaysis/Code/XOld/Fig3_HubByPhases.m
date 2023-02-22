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


H=figure(3)
%bar(allData(:,2:3))
bar(allData(:,2:3))
legend(GroupNames{2:3})
xlabel('Disease Phases')
ylabel('Fraction of Nodes Associated with Disease Hubs')

print(H,fullfile(saveDirBase,['fig3_' mName '.tif']),'-dpng','-r400');



for s=1:6
for j = 2:3
    
    if j == 2
        StagingSel=T_tauStagingLaus;
    else 
        StagingSel=T_tdpStagingLaus;
    end
    
    
    vals = Hubs.FAGraph.([(GroupNames{1}) 'Thresh']).(GroupNames{1})' & StagingSel==s;
    allData(s,j,1) = sum(vals)/sum(StagingSel==s);%each subject/region value is a datapoint
end
end


H=figure(4)
%bar(allData(:,2:3))
bar(allData(:,2:3))
legend(GroupNames{2:3})
xlabel('Disease Phases')
ylabel('Fraction of Nodes Associated with Control Hubs')

print(H,fullfile(saveDirBase,['fig3_HC_' mName '.tif']),'-dpng','-r400');

