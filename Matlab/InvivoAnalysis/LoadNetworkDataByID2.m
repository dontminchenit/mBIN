function [AllResults] = LoadNetworkDataByID(ID,allThickVals,outMatSaveFile,atlasName,numLab)

N = length(ID);
AllResults = struct;
 AllResults.Thickness.Mean = nan(N,numLab);
 AllResults.Thickness.Median = nan(N,numLab);

 AllResults.Volume.Total = nan(N,numLab);
 AllResults.Volume.ICV = nan(N,numLab);
 
for i = 1:N
    i
    currID = ID(i);
    currDates=allThickVals{allThickVals.id == currID ,'date'};
    uniqDates = unique(currDates);

    num2str(currID)
    
    strcmpi(allThickVals.date,uniqDates{1});
 
    currVolume = allThickVals{allThickVals.id == currID & strcmpi(allThickVals.date,uniqDates{1}) & strcmpi(allThickVals.system,atlasName),'volume'};
    currICV  = allThickVals{allThickVals.id == currID & strcmpi(allThickVals.date,uniqDates{1}) & strcmpi(allThickVals.system,atlasName),'icv'}; 
    currLabelsVol = allThickVals{allThickVals.id == currID & strcmpi(allThickVals.date,uniqDates{1}) & strcmpi(allThickVals.system,atlasName),'label'};
    
    for L = 1:numLab
        if(sum(currLabelsVol==L)==1)
            if(~isempty(currVolume))
            AllResults.Volume.Total(i,L) = currVolume(currLabelsVol==L);
            AllResults.Volume.ICV(i,L) = currICV(currLabelsVol==L);
            end
        else
            disp(['Missing Volume Label ' num2str(L) ' ID ' num2str(currID)])
        end
    end
    
end

save(outMatSaveFile,'AllResults');
