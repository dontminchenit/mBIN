%Must run Fig4ab_EpiAdjStageTrans_zscore.m first
%Set Epi0Adj1Name flag in Fig4ab_EpiAdjStageTrans_zscore.m to change between 
%Epicenter (Figure 4c) or Adjacent (Figure 4d) analysis
%Make sure to Set ModeOrig1_HubAnalysis2=2 for this analysis
SN=5; % of stages
%allData = indvWeakEdgeCountNormalized;
%allData = indvWeakEdgeCount;
labelNames={'Tau','TDP'}

allData =PrctHubsEdgeCount;
%labelNames={'In','Out'}
%allData(:,1,:)=sum(allData,2);


x = 1:SN;
y=allData;
H3=figure(3)
clf
if(flagEpi0Adj1)
    names={'1->2','2->3','3->4','4->5','5->6'};
else
    names={'1->2','1->3','1->4','1->5','1->6'};
end
for s=1:SN
    [H pval]=ttest2(allData(:,s,1),allData(:,s,2))
    
    %pval = ranksum(allData(:,s,1),allData(:,s,2))
    if(pval<0.05)
        names{s} = [names{s} '*'];
    end
end

x = names;
y = allData*100;
iosr.statistics.boxPlot(x,y, ...
    'style','hierarchy',...
    'xSeparator',true,...
    'symbolColor','k',...
    'symbolMarker','.',...
    'medianColor','k',...
    'boxcolor','auto',...
    'showScatter',true,...
    'scatterColor', 'r', ...
    'scatterMarker', '.', ...
    'groupLabels', labelNames, ...
    'groupLabelFontSize', 10 ...
);
set(gca,'FontSize',20)
ytickformat('percentage')
print(H3,fullfile(saveDirBase,['fig6_WeakEdgeCmp_' Epi0Adj1Name '_' hubName '.tif']),'-dpng','-r400');

% cmap = get(0, 'defaultaxescolororder');
% for ii = 1:2
%     structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
%         'markeredgecolor', cmap(ii,:)), h);
% end
% set([h.lwhis h.uwhis], 'linestyle', '-');
% set(h.out, 'marker', '.');
% legend({'Tau','TDP'},'Location','northeast');
% set(gca,'XTick',1:SN,'XTickLabel',names)
% title([GroupNames{G1} ' Using ' hubName])