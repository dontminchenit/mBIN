

H=figure(2)
clf
panelAll = panel();
panelAll.pack(2,3);

panelAll.de.margin = 15;
panelAll.marginbottom = 15;
panelAll.margin = [15 15 15 15];


allMetricData = cell(1,5);
for m = 1:5
%cross reference with metrics 
SN=6; % of stages
allData = nan(SN,2,10000);
allDataSubMean = nan(SN,2,10000);
allDataRegionMean = nan(SN,2,10000);

mName = measureNames{m}

for s=1:6
for j = 2:3
    
    if j == 2
        StagingSel=T_tauStagingLaus;
    else 
        StagingSel=T_tdpStagingLaus;
    end
    
    normaMean = nanmean(AllResults.FAGraph.(mName)(GroupInd{1},StagingSel==s),1); 
    vals = AllResults.FAGraph.(mName)(GroupInd{j},StagingSel==s)-normaMean;
    allData(s,j-1,1:length(vals(:))) = vals(:);%each subject/region value is a datapoint
    
    valsSubMean = nanmean(vals,2);%average region, each subject is datapoint
    valsRegionMean = nanmean(vals,1);%average subjects, each region is a datapoint
    
    allDataSubMean(s,j-1,1:length(valsSubMean(:))) = valsSubMean(:);
    allDataRegionMean(s,j-1,1:length(valsRegionMean(:))) = valsRegionMean(:);
end
end


allMetricData{m}=allDataSubMean;

for k = 2%:3

x = 1:SN;
switch k
    case 1
        y=allData;        
        savename = 'All Data';
    case 2
        y=allDataSubMean;     
        savename = 'Subject Data (Regions Averaged)';
    case 3
        y=allDataRegionMean;
        savename = 'Region Data (Subjects Averaged)';
        
end
% Plot boxplots

if(m<=3)
    panelAll(1,m).select();
elseif(m>3)
    panelAll(2,m-3).select();
end

h = boxplot2(y,x,'Notch','off');

cmap = get(0, 'defaultaxescolororder');
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), h);
end
set([h.lwhis h.uwhis], 'linestyle', '-');
set(h.out, 'marker', '.');


% title(savename)
 ylabel(mName)
 xlabel('Disease Phases')
 set(gca,'xtick',1:SN)
 xtickangle(45)
 if(m==1)
legend({'FTLD-Tau','FTLD-TDP'},'Location','northeast');
 end

end
end
print(H,fullfile(saveDirBase,['fig7b_' mName '.tif']),'-dpng','-r600');