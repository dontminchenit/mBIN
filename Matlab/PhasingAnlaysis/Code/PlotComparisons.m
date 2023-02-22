

SN=5; % of stages
allData = nan(SN,2,30000);



for s=1:SN
    for q=1:2
    
    vals=squeeze(inoutHubData{s,p,q});
    allData(s,q,1:length(vals(:))) = vals(:);%each subject/region value is a datapoint
    
    end
end


x = 1:SN;
y=allData;
H3=figure(3)
clf
if(flagEpi0Adj1)
    names={'1->2','2->3','3->4','4->5','5->6'}
else
    names={'1->2','1->3','1->4','1->5','1->6'}
end
for s=1:SN
    
    if(inoutHubPtest(s,p)<0.05)
        names{s} = [names{s} '*'];
    end
end

h = boxplot2(y,x,'Notch','off');


cmap = get(0, 'defaultaxescolororder');
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), h);
end
set([h.lwhis h.uwhis], 'linestyle', '-');
set(h.out, 'marker', '.');
legend({'Hubs','NonHubs'},'Location','northeast');
set(gca,'XTick',1:SN,'XTickLabel',names)
title([GroupNames{G1} ' Using ' hubName])