function [PathCenterVal]=findPathCenterCount(PathVal, PathLabel)

PathCenterVal = nan(size(PathLabel));

for k = 1:size(PathLabel,1)
   
    currSubPath = PathLabel(k,:);
    unqLabel= unique(currSubPath);
    unqLabel(isnan(unqLabel))=[];
    
    for i = 1:length(unqLabel)
        loc = ceil(median(find(currSubPath==unqLabel(i))));
        currLabVal = mean(double(PathVal(k,currSubPath==unqLabel(i))),'omitnan');
        double(PathVal(k,currSubPath==unqLabel(i)))
        PathCenterVal(k,loc) = currLabVal;
        
    end
    
    
end