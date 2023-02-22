macBase = fullfile('/Users','Min','Dropbox (Personal)','Research','Projects','FTD','PhasingAnalysis');
currBase = 'D:\Min\Dropbox (Personal)\Research\Publications\StagingConnectivity_2021\BrainSubmission\R1';

inInfo = fullfile(currBase,'mc_revisions_bvols_mmse_cdr.csv');
T = readtable(inInfo);

TauIND2 = strcmpi(T.Group,'tau');
TDPIND2 = strcmpi(T.Group,'tdp');

a=T.MMSETotal(TauIND2);
a(a<=1)=nan;
b=T.MMSETotal(TDPIND2);
b(b<=1)=nan;
[H,P]=ttest2(a,b,'tail','both')

boxplot(T.MMSETotal,T.Group,'notch','on')