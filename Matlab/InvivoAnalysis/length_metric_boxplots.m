metricTau=betweenness_wei(thickcovMatTau_delta);
metricTDP=betweenness_wei(thickcovMatTDP_delta);

metricTau=strengths_und(thickcovMatTau_delta)';
metricTDP=strengths_und(thickcovMatTDP_delta)';

metricTau=CoMDistMatrix(logical(triu(thickcovMatTau_delta)));
metricTDP=CoMDistMatrix(logical(triu(thickcovMatTDP_delta)));

metricTau(metricTau==0)=nan;
metricTDP(metricTDP==0)=nan;

X=[metricTau;metricTDP]
G=[ones(size(metricTau));2*ones(size(metricTDP))]


[h p]=ttest2(metricTau,metricTDP)
figure(12)
boxplot(X,G,'Labels',{'FTLD-tau','FTLD-TDP'})