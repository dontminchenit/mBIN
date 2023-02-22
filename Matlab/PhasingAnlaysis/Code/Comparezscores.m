
NTau=size(TAUIndvZscore,3);
TauZScore=nan(NTau,1);

for i = 1:NTau
    CurrWeakEdge = TAUIndvZscore(:,:,i);
    TauZScore(i) = mean(CurrWeakEdge(:),'omitnan');
end

NTDP=size(TDPIndvZscore,3);
TDPZscores=nan(NTDP,1);

for i = 1:NTDP
    CurrWeakEdge = TDPIndvZscore(:,:,i);
    TDPZscores(i) = mean(CurrWeakEdge(:),'omitnan');
end

indata=cell(1,2);
indata{1}=TauZScore;
indata{2}=TDPZscores;

boxPlotAdv(indata)
[h, p] = ttest2(TauZScore,TDPZscores)