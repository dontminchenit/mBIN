for m = 1:5
%cross reference with metrics 
allData = nan(5,2,10000);
allDataSubMean = nan(5,2,10000);
allDataRegionMean = nan(5,2,10000);

mName = measureNames{m}

for s=1:5
for j = 2:3
    
    if j == 2
        StagingSel=T_tauStagingLaus;
    else 
        StagingSel=T_tdpStagingLaus;
    end
    
    vals = AllResults.FAGraph.(mName)(GroupInd{j},StagingSel==s);
    allData(s,j-1,1:length(vals(:))) = vals(:);%each subject/region value is a datapoint
    
    valsSubMean = nanmean(vals,2);%average region, each subject is datapoint
    valsRegionMean = nanmean(vals,1);%average subjects, each region is a datapoint
    
    allDataSubMean(s,j-1,1:length(valsSubMean(:))) = valsSubMean(:);
    allDataRegionMean(s,j-1,1:length(valsRegionMean(:))) = valsRegionMean(:);
end
end


H=figure(2)
clf
panelAll = panel();
panelAll.pack(1,3);

panelAll.de.margin = 15;
panelAll.marginbottom = 15;
panelAll.margin = [15 15 15 15];


for k = 1:3

x = 1:5;
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

panelAll(1,k).select();

h = boxplot2(y,x,'Notch','on');

cmap = get(0, 'defaultaxescolororder');
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), h);
end
set([h.lwhis h.uwhis], 'linestyle', '-');
set(h.out, 'marker', '.');


 title(savename)
 ylabel(mName)
 xlabel('Stages')
 set(gca,'xtick',1:5)
 xtickangle(45)
legend({'tau','TDP'},'Location','northeast');


end
print(H,fullfile(saveDirBase,['fig7b_' mName '.tif']),'-dpng','-r400');
end