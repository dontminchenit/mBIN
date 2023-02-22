function [PathWeakEdge,PathWeakEdgeReg, PathIndvZscore, PathMeanZscore, PathIndvZscoreReg] = calcWeakEdges(Pathmtxs,HCmtxs,threshMatrix)

%find the mean FA for Group 1
%HCmtxs=AllResults.(graphNames{2}).('WeightMtx')(GroupInd{1});
fullHCmtx=zeros(size(HCmtxs{1},1),size(HCmtxs{1},2),length(HCmtxs));%This stores all the matrices for Group 1 (HC)
for k = 1:length(HCmtxs)
    currMtx=HCmtxs{k};
    currMtx(currMtx==0)=nan;%set zero FA as NaN, since no measurements
    fullHCmtx(:,:,k) = currMtx;
end

%find the mean FA for Group 2
%Pathmtxs=AllResults.(graphNames{g}).('WeightMtx')(GroupInd{2});
fullPathmtx=zeros(size(Pathmtxs{1},1),size(Pathmtxs{1},2),length(Pathmtxs));%This stores all the matrices for Group 2 (Tau)
for k = 1:length(Pathmtxs)
    currMtx = Pathmtxs{k};
    currMtx(currMtx==0)=nan;%set zero FA as NaN, since no measurements
    fullPathmtx(:,:,k)=currMtx;
end

%Calculate mean of both groups, and std of control group
%(for each edge)
meanHCmtx = mean(fullHCmtx,3,'omitnan');
ctrlCount = sum(~isnan(fullHCmtx),3);

meanPathmtx = mean(fullPathmtx,3,'omitnan');
stdHCmtx = std(fullHCmtx,0,3,'omitnan');
stdHCmtx(stdHCmtx==0)=nan;%if we have a stddev of zero, it probably means only a single data point, so exclude
stdHCmtx(ctrlCount<5)=nan;

%find the zscore of disease mean relative to controls
PathMeanZscore=(meanPathmtx-meanHCmtx)./stdHCmtx;

%find the zscore of all subjects
PathIndvZscore = nan(size(fullPathmtx));
PathWeakEdge = nan(size(fullPathmtx));
PathWeakEdgeReg = nan(size(fullPathmtx,2),size(fullPathmtx,3)); 
PathIndvZscoreReg = nan(size(fullPathmtx,2),size(fullPathmtx,3)); 

for k = 1:size(fullPathmtx,3)
    PathIndvZscore(:,:,k) = (fullPathmtx(:,:,k)-meanHCmtx)./stdHCmtx;

    PathWeakEdgeCurr=double(PathIndvZscore(:,:,k)<-3);
    PathWeakEdgeCurr(PathWeakEdgeCurr==0)=nan;
    PathWeakEdgeCurr(isnan(stdHCmtx))=nan;
    
    PathWeakEdge(:,:,k)=PathWeakEdgeCurr;
    
    PathWeakEdgeReg(:,k) = sum(PathWeakEdgeCurr,1,'omitnan');
    PathIndvZscoreReg(:,k)=mean(PathIndvZscore(:,:,k).*threshMatrix,1,'omitnan');
end

