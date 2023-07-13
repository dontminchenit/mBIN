
%AllCasesInd = [ GroupInd{2}; GroupInd{3}]
%AllCaseMtx = [AllResults.FAGraph.WeightMtxNormal(GroupInd{2}); AllResults.FAGraph.WeightMtxNormal(GroupInd{3})]

baseData=SubVolsThick;

data1 = baseData(strcmpi(subDemo.Group,'tau'),:);
data2 = baseData(strcmpi(subDemo.Group,'tdp'),:);
dataAll= [data1; data2];

N=height(dataAll);
allCorr = nan(N,N);
for i = 1:N
    for j = 1:N

        if(i ~= j)

        %figure(2)
        %clf
        %scatter(a(:),b(:))
        %pause
        a=dataAll(i,:);
        b=dataAll(j,:);
        
        %[RHO,PVAL] = corr(a(:),b(:));
        RHO = icc(a,b);
        allCorr(i,j)=RHO;
        end
    end
end

figure(1)
clf
imagesc(allCorr)