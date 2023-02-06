function AllResults = calcNetworkMetrics_bin(currAdjMtx)

        currAdjMtx=triu(currAdjMtx)+triu(currAdjMtx,1)';
        currAdjMtx(currAdjMtx<0) = 0;
        currAdjMtx(isnan(currAdjMtx)) = 0;


        AllResults=struct;
        thresh=.7;
        currAdjMtx = currAdjMtx > thresh;
        
        W=currAdjMtx+currAdjMtx';
        W_normalized = weight_conversion(W, 'normalize');
        L = weight_conversion(W, 'lengths');
        
        %calc the local measures
        AllResults.('WeightMtx') = W;
        AllResults.('WeightMtxNormal') = W_normalized;
        AllResults.('LengthMtxNormal') = L;
        AllResults.('BetweenCen') = betweenness_bin(L);
        AllResults.('ClusterCoeff') = clustering_coef_bu(W_normalized);
        AllResults.('LocalEff') = efficiency_bin(W_normalized, 2);
        AllResults.('EigenCen') = eigenvector_centrality_und(W);
        degree = sum(W,2);
        AllResults.('Degree') = degree;