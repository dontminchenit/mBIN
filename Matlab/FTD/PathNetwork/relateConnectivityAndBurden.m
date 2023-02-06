function NetworkBurdenCorr=relateConnectivityAndBurden(AdjMtx,Burden,MarkerNames)


AllMetrics = calcNetworkMetrics(AdjMtx);
%AllMetrics = calcNetworkMetrics_bin(AdjMtx);



mName = {'Degree','ClusterCoeff','BetweenCen',};
mFullName = {'Degree', 'Clustering\nCoefficient', 'Betweenness\nCentrality'};
%H=figure(10*PlotID)
NetworkBurdenCorr = struct;
NetworkBurdenCorr.('mName')=mName;
NetworkBurdenCorr.('mFullName')=mFullName;
NetworkBurdenCorr.RegionNames=MarkerNames;

X=Burden;
NetworkBurdenCorr.Burden = Burden;

for i = 1:3
%    subplot(3,3,ColID + 3*(i-1));
    
    Y=AllMetrics.(mName{i});
    valid=~isnan(X') & ~isnan(Y);
    
%    scatter(X,Y,25)
    
    %            regPlot(Burden',AllMetrics.(mName{i}))
%    if (ColID==1)
%        ylabel(sprintf(mFullName{i}))
%    end
    
%     if(i==3)
%         switch ColID
%             
%             case 1
%                 xlabel('GM %AO')
%                 
%             case 2
%                 xlabel('WM %AO')
%                 
%             case 3
%                 xlabel('GM & WM %AO')
%         end
%     end
%     lsline
    [R P LB UB]=corrcoef(X(valid),Y(valid));
    [rho,pval] = corr(X(valid)',Y(valid))
    NetworkBurdenCorr.(mName{i})=Y;
    NetworkBurdenCorr.([mName{i} '_corr'])=rho%R(1,2);
    NetworkBurdenCorr.([mName{i} '_corrPval'])=pval%P(1,2);
    NetworkBurdenCorr.([mName{i} '_corrLB'])=0%LB(1,2);
    NetworkBurdenCorr.([mName{i} '_corrUB'])=0%UB(1,2);
    NetworkBurdenCorr.([mName{i} '_corrN'])=sum(valid(:));
    
    %[b,bint,r,rint,stats] = regress(X(valid)',[Y(valid) ones(length(Y(valid)),1)])
    [b,bint,r,rint,stats] = regress(Y(valid),[X(valid)' ones(length(X(valid)),1)],.05)
    
    NetworkBurdenCorr.([mName{i} '_regcoeff'])=b;
    NetworkBurdenCorr.([mName{i} '_regpval'])=stats(3);
    
    
    
%     str=sprintf('r= %1.2f',tmp(1,2));
%     T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
%     set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
end

end


