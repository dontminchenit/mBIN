function compareTauTDP(tau,tdp,titleName)

    [H,P]=ttest2(tau,tdp)
    x = tau;
    y = tdp;
    z = [x; y];
    g = [zeros(length(x),1); ones(length(y),1)];
    figure(10)
    clf
    boxplot(z,g,'labels',{'tau','TDP'});
    title(titleName);

