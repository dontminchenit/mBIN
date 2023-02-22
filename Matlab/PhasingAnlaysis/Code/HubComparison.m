HubsTauLost = ((Hubs.FAGraph.HCThresh.(GroupNames{2})' - Hubs.FAGraph.HCThresh.(GroupNames{1})')<0);
HubsTDPLost = ((Hubs.FAGraph.HCThresh.(GroupNames{3})' - Hubs.FAGraph.HCThresh.(GroupNames{1})')<0);


TauMeanGMPathHub = mean(TauGMPath(:,HubsTauLost),1,'omitnan');
TDPMeanGMPathHub = mean(TDPGMPath(:,HubsTDPLost),1,'omitnan');

TauMeanWMPathHub = mean(TauWMPath(:,HubsTauLost),1,'omitnan');
TDPMeanWMPathHub = mean(TDPWMPath(:,HubsTDPLost),1,'omitnan');

indata=cell(1,2);
indata{1}=TauMeanGMPathHub;
indata{2}=TauMeanWMPathHub;

boxPlotAdv(indata)
[h, p] = ttest2(TauMeanGMPathHub,TauMeanWMPathHub)