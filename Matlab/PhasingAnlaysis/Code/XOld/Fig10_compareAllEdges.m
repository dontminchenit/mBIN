%find the mean FA for Group 1
HCmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{1});

%find the mean FA for Group 2
TAUmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{2});

%find the mean FA for Group 3
TDPmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{3});

prctVal=0;

thresh = prctile(CoMDistMatrix(:),prctVal);

threshMatrix = double(CoMDistMatrix > thresh);
threshMatrix(threshMatrix~=1)=nan;

%figure(5)
%imagesc(threshMatrix)

[TAUWeakEdge, TauWeakEdgeReg, TAUIndvZscore, TAUMeanZscore, TAUIndvZscoreReg] = calcWeakEdges(TAUmtxs,HCmtxs,threshMatrix);
[TDPWeakEdge, TDPWeakEdgeReg, TDPIndvZscore, TDPMeanZscore, TDPIndvZscoreReg] = calcWeakEdges(TDPmtxs,HCmtxs,threshMatrix);

figure(9)

%indata{1}=mean(TAUIndvZscoreReg,1,'omitnan')';
%indata{2}=mean(TDPIndvZscoreReg,1,'omitnan')';


%indata{1}=mean(TAUIndvZscoreReg,2,'omitnan');
%indata{2}=mean(TDPIndvZscoreReg,2,'omitnan');


indata{1}=sum(TauWeakEdgeReg,1)';
indata{2}=sum(TDPWeakEdgeReg,1)';


%indata{1}=sum(TauWeakEdgeReg,2);
%indata{2}=sum(TDPWeakEdgeReg,2);

boxPlotAdv(indata)
%histogram(indata{1}-indata{2},100);
%xlabel('Avg Z-score Difference (Tau - TDP43) at Node')
%ylabel('Count')
%title('Histogram of 
[h, p] = ttest2(indata{1},indata{2})