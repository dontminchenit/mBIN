

H=figure(21)
scatter(mean_Tau_w,Tau_Degree)
[RHO,PVAL]=corr(mean_Tau_w',Tau_Degree')
xlabel('Mean W-score')
ylabel('Degree')
lsline
title([saveName ', Tau r=' num2str(RHO) ' p=' num2str(PVAL)])
print(H,fullfile(baseSaveDir,[ saveName '_TAU_wscoreVsDEg.png']),'-dpng','-r400');

H=figure(22)
scatter(mean_TDP_w,TDP_Degree)
[RHO,PVAL]=corr(mean_TDP_w',TDP_Degree')
xlabel('Mean W-score')
ylabel('Degree')
lsline
title([saveName ', TDP r=' num2str(RHO) ' p=' num2str(PVAL)])
print(H,fullfile(baseSaveDir,[ saveName '_TDP_wscoreVsDEg.png']),'-dpng','-r400');