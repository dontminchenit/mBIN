UniqueScans=unique(thicknessAllraw(:,{'id', 'date'}));

%UniqueIDs=unique(UniqueScans.id);

UniqueIDs=thickIDs(demoT.LBD0_LBDAD1 == 0)%|(demoT.LBD0_LBDAD1 == 0));

numScans=zeros(height(UniqueIDs),1);

for i = 1:height(UniqueIDs)
currID=UniqueIDs(i);
numScans(i)=sum(UniqueScans.id==currID);

end

figure(20)
hist(numScans)
xlabel('Number of Scans')
ylabel('Number of Subjects With')
