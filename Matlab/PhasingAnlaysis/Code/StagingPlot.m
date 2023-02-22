function StagingPlot(tableIn,plotType,transType,saveDirBase,saveName,MeasureName)
H=figure(25)
clf
x = 1:size(tableIn,1);
y=tableIn;

% Plot boxplots

h = boxplot2(y,x,'notch','on');

% Alter linestyle and color

cmap = get(0, 'defaultaxescolororder');
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), h);
end
set([h.lwhis h.uwhis], 'linestyle', '-');
set(h.out, 'marker', '.');


 ylabel(MeasureName)
 
 if(size(tableIn,1)==5)
     set(gca,'xtick',1:5, 'xticklabel',{'1','2','3','4','5'})
 else
 if(transType == 1)
     set(gca,'xtick',1:4, 'xticklabel',{'1-2','2-3','3-4','4-5'})
 else
     set(gca,'xtick',1:4, 'xticklabel',{'1-2','1-3','1-4','1-5'})
 end
 end
 xtickangle(45)
legend({'tau','TDP'},'Location','northeast');

print(H,fullfile(saveDirBase,[saveName  '_Transition.tif']),'-dpng','-r400'); 