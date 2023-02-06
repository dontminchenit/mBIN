function pvals = PlotMetricComparisons(TauMetrics,TDPMetrics)


        
        MetricNames = {'ClusterCoeff','LocalEff','EigenCen'};
        %MetricNames = {'LocalEff'}%'EigenCen'};
        
        
        Nmax = 20;

        MN = length(MetricNames);
        allData = nan(MN,2,Nmax);

        pvals=ones(1,MN);

        for i = 1:MN
           
            val = TauMetrics.(MetricNames{i});
            allData(i,1,1:length(val)) = val;
            
            val2 = TDPMetrics.(MetricNames{i});
            allData(i,2,1:length(val2)) = val2;
            
            [H,P]=ttest(val,val2);
            pvals(i)=P;
            
        end
        
        
        
        
        
        x = 1:MN;
        y=allData;
        % Plot boxplots

        h = boxplot2(y,x,'notch','on');

        % Alter linestyle and color

        cmap = get(0, 'defaultaxescolororder');
        for ii = 1:2
            structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
                'markeredgecolor', cmap(ii,:)), h);
        end
        set([h.lwhis h.uwhis], 'linestyle', '-');
        set(h.out, 'marker', '.');


         ylabel('Measure')
       %  ylim([0 .9])
         set(gca,'xtick',1:MN, 'xticklabel',MetricNames)
         xtickangle(45)
         legend({'Tau','TDP43'},'Location','southeast');
        
        

end