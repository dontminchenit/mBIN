
H=figure(1)
clf
x = 1:4;
y=FAByStage;
y=MetricBySTage;

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


 ylabel('FA')
 set(gca,'xtick',1:4, 'xticklabel',{'1-2','1-3','1-4','1-5'})
 xtickangle(45)
legend({'tau','TDP'},'Location','northeast');

print(H,fullfile(saveDirBase,'fig6FAEpiTransition.tif'),'-dpng','-r400'); 