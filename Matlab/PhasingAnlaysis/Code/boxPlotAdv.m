function boxPlotAdv(inData)

c_1=inData{1,1};
c_2=inData{1,2};
C = [c_1; c_2];
grp = [zeros(1,length(c_1)),ones(1,length(c_2))];
names={'HubFA', 'NonHubFA'};
names={'Tau', 'TDP'};
boxplot(C,grp,'labels',names,'whisker',2)
%ylabel('Average Z-score of Edges Connected to Node')
ylabel('Sum # of Weak Edges Connected to Node')


end