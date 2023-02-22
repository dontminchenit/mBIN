inData = '/Users/min/Dropbox (Personal)/Research/Projects/PhasingAnalysis/Data/DTI cohort ordinal data.xlsx'

T=readtable(inData);

tauIdx = strcmpi(T.Group,'tau')
tdpIdx = strcmpi(T.Group,'tdp')



tau=T.MFTau(tauIdx);
tdp=T.MFTDP43(tdpIdx);
compareTauTDP(tau,tdp,'MFC');

tau=T.SMTTau(tauIdx);
tdp=T.SMTTDP43(tdpIdx);
compareTauTDP(tau,tdp,'SMT');