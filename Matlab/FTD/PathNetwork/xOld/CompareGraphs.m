



TauSigLower=double(TauMtx<TDPLBmtx);

TauSigLower(isnan(TauMtx))=nan;

TDPSigLower=double(TDPMtx<TauLBmtx);

TDPSigLower(isnan(TDPMtx))=nan;



for k = 1:2
    if k==1
plotMtx = TauSigLower;
saveName = 'tauLower_comp';    
    else
plotMtx = TDPSigLower;
saveName = 'tdpLower_comp';    
    end

plotGraphs
end