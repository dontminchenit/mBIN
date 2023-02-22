%LoadDataSarah

%comparison Against GlobalEFf

pvals = nan(2,4);

for i = 1:2
    
    H=figure(i);
    clf
    if i ==1
        mName = 'GlobalEff';
        labelName = 'Global Efficiency';
    else
        mName = 'MeanStr';
        labelName = 'Mean Nodal Strength';
    end
    
    Nmax = max([sum(HCIND) sum(TauIND) sum(TDPIND)]);
    allData = nan(1,3,Nmax);
    
    cmp1=[2 3 2 3];
    cmp2=[1 1 3 2];
    
    for j = 1:4
        vals = AllResults.FAGraph.(mName)(GroupInd{cmp1(j)});
        vals2 = AllResults.FAGraph.(mName)(GroupInd{cmp2(j)});
        
        [h,p,ci,stats] = ttest2(vals,vals2,'Tail','left');
        pvals(i,j) = p;
    end
    
    
    
    
    for j = 1:3
        vals = AllResults.FAGraph.(mName)(GroupInd{j});
        allData(1,j,1:length(vals)) = vals;
        
        %     vals = SubnetworkResultsSV.FAGraph.(mName)(GroupInd{j});
        %     allData(2,j,1:length(vals)) = vals;
        %
        %     vals =SubnetworkResultsNA.FAGraph.(mName)(GroupInd{j});
        %     allData(3,j,1:length(vals)) = vals;
        %
        %     vals =nonSubnetworkResultsSV.FAGraph.(mName)(GroupInd{j});
        %     allData(4,j,1:length(vals)) = vals;
        %
        %     vals =nonSubnetworkResultsNA.FAGraph.(mName)(GroupInd{j});
        %     allData(5,j,1:length(vals)) = vals;
    end
    
    
    
    
    x = 1;
    y=allData;
    % Plot boxplots
    
    h = boxplot2(y,x);
    
    % Alter linestyle and color
    
    cmap = get(0, 'defaultaxescolororder');
    for ii = 1:3
        structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
            'markeredgecolor', cmap(ii,:)), h);
    end
    set([h.lwhis h.uwhis], 'linestyle', '-');
    set(h.out, 'marker', '.');
    
    
    ylabel(labelName)
    set(gca,'xtick',[0.75 1 1.25], 'xticklabel',{'HC','FTLD-Tau','FTLD-TDP'})
    xtickangle(45)
    %legend({'HC','Tau','TDP'},'Location','northeast');
    
    print(H,fullfile(saveDirBase,['fig2_' num2str(i) '.tif']),'-dpng','-r400');
end




%
%
%