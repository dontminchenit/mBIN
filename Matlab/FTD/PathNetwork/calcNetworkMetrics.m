function AllResults = calcNetworkMetrics(currAdjMtx)

AllResults=struct;


        currAdjMtx=triu(currAdjMtx)+triu(currAdjMtx,1)';
        
        currAdjMtx(currAdjMtx<0) = 0;
        currAdjMtx(isnan(currAdjMtx)) = 0;
        
        W=currAdjMtx+currAdjMtx';
        W_normalized = weight_conversion(W, 'normalize');
        L = weight_conversion(W, 'lengths');
        
        %calc the local measures
        AllResults.('WeightMtx') = W;
        AllResults.('WeightMtxNormal') = W_normalized;
        AllResults.('LengthMtxNormal') = L;
        AllResults.('BetweenCen') = betweenness_wei(L);
        AllResults.('ClusterCoeff') = clustering_coef_wu(W_normalized);
        AllResults.('LocalEff') = efficiency_wei(W_normalized, 2);
        %AllResults.('EigenCen') = eigenvector_centrality_und(W);
        degree = sum(W,2);
        AllResults.('Degree') = degree;