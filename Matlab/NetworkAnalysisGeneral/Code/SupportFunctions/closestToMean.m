function closestToCoM=closestToMean(pointsList)

CoM = mean(pointsList,1);

[D I] = pdist2(pointsList,CoM,'euclidean','smallest',1);

closestToCoM = pointsList(I,:);
