function inMtx=mtxAdjustMultCmp(inMtx,pthresh)


    notNaN = ~isnan(inMtx);
    [padj_fdr,alpha_fdr] = multicmp(inMtx(notNaN),'fdr',pthresh);
    inMtx(notNaN)=padj_fdr<pthresh;
    inMtx(~notNaN)=0;
