
NTau=size(TAUWeakEdge,3);
TauEdgeCount=nan(NTau,1);

for i = 1:NTau
    CurrWeakEdge = TAUWeakEdge(:,:,i)==1;
    TauEdgeCount(i) = sum(CurrWeakEdge(:));
end

NTDP=size(TDPWeakEdge,3);
TDPEdgeCount=nan(NTDP,1);

for i = 1:NTDP
    CurrWeakEdge = TDPWeakEdge(:,:,i)==1;
    TDPEdgeCount(i) = sum(CurrWeakEdge(:));
end

indata=cell(1,2);
indata{1}=TauEdgeCount;
indata{2}=TDPEdgeCount;

boxPlotAdv(indata)
[h, p] = ttest2(TauEdgeCount,TDPEdgeCount)