LoadData


globMeasureNames={'GlobalEff','MeanStr','CharPathLen'};

LTests = cell(1,4);

LTests{1}=T.BNT_BL;
LTests{2}=T.FWordsCorrect_BL;
LTests{3}=T.AnimalsTotal_BL;
LTests{4}=T.CookieFluency_BL;

corrMtx = nan(2,4,3);%Dim 1-Measure, 2-Test, 3-Group)
corrMtxP = nan(2,4,3);%Dim 1-Measure, 2-Test, 3-Group)

Group = cell(1,3);
Group{1}=HCIND;
Group{2}=svIND;
Group{3}=naIND;

TestNames={'BNT','FWordsCorrect','AnimalsTotal','CookieFluency'};


saveDirBaseSub=fullfile(saveDirBase,'fig3');
mkdir(saveDirBaseSub);
for q = 1:5

switch q
    case 1
        CurrResults = AllResults;
        NodeName='AllNodes';
        gRange = 1:3;
    case 2
        CurrResults = SubnetworkResultsSV;
        NodeName='svAtrophied';
        gRange = 2;
    case 3
        CurrResults = SubnetworkResultsNA;
        NodeName='naAtrophied';
        gRange = 3;
    case 4
        CurrResults = nonSubnetworkResultsSV;
        NodeName='svNonAtrophied';
        gRange = 2;
    case 5
        CurrResults = nonSubnetworkResultsNA;
        NodeName='naNonAtrophied';
        gRange = 3;
end


%x1=SubnetworkResultsSV.FAGraph.(mName)(svIND);
%x2=SubnetworkResultsNA.FAGraph.(mName)(naIND);
%x3=nonSubnetworkResultsSV.FAGraph.(mName)(svIND);
%x4=nonSubnetworkResultsNA.FAGraph.(mName)(naIND);


for i = 2
for g = gRange
    for j = 1:2
        measure = CurrResults.(graphNames{i}).(globMeasureNames{j});
        measure(isinf(measure))=nan;
        %we have the measure for each subgroup
        for t = 1:4
            test = LTests{t};
            [r p]=corr(measure(Group{g}),test(Group{g}),'Rows','complete');
            corrMtx(j,t,g)= r;
            corrMtxP(j,t,g)= p;
        end
    end
    H=figure(1)
        imagesc(corrMtx(:,:,g))
        caxis([-1 1])
        colorbar
        set(gca,'xtick',[1:4],'xticklabel',TestNames)
        set(gca,'ytick',[1:3],'yticklabel',globMeasureNames)
        saveas(H,fullfile(saveDirBaseSub,[NodeName '_' GroupNames{g} '_corrMtx.png']));
    H=figure(2)
        imagesc(corrMtxP(:,:,g)<0.05)
        set(gca,'xtick',[1:4],'xticklabel',TestNames)
        set(gca,'ytick',[1:3],'yticklabel',globMeasureNames)
        saveas(H,fullfile(saveDirBaseSub,[NodeName '_' GroupNames{g} '_corrMtxP.png']));
end
end
end


