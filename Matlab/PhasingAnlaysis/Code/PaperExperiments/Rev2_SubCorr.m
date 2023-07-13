
AllCasesInd = [ GroupInd{2}; GroupInd{3}]

AllCaseMtx = [AllResults.FAGraph.WeightMtx(GroupInd{2}); AllResults.FAGraph.WeightMtx(GroupInd{3})]
N=length(AllCaseMtx)

allCorr = nan(N,N);
for i = 1:N
    for j = 1:N

        if(i ~= j)

        a=AllCaseMtx{i};
        a(a==0)=nan;
        b=AllCaseMtx{j};
        b(b==0)=nan;
%         figure(2)
%         clf
%         scatter(a(:),b(:))
%         pause
        [RHO,PVAL] = corr(a(:),b(:),'rows','Complete');
        allCorr(i,j)=RHO;
        end
    end
end

figure(1)
clf
imagesc(allCorr)