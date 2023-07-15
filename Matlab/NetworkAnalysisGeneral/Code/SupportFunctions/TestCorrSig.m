function [Probs, zObs]= TestCorrSig(Corr1,Corr2,NCorr1,NCorr2)

zCorr1=atanh(Corr1);
zCorr2=atanh(Corr2);
zObs = (zCorr1 - zCorr2)./ ((1 ./ (NCorr1 - 3)) + (1 ./ (NCorr2 - 3))).^(.5);
Probs = normcdf(zObs);
