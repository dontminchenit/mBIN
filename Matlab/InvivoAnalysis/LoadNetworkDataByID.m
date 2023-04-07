function [AllResults] = LoadNetworkDataByID(ID,allThickVals,outMatSaveFile,atlasName,numLab)

N = length(ID);
AllResults = struct;
 AllResults.Thickness.Mean = nan(N,numLab);
 AllResults.Thickness.Median = nan(N,numLab);
 AllResults.Volume.Total = nan(N,numLab);

for i = 1:N
    i
    currID = ID(i);
    currDates=allThickVals{allThickVals.id == currID ,'date'};
    uniqDates = unique(currDates);

    num2str(currID)
    
    strcmpi(allThickVals.date,uniqDates{1});

    currThickness = allThickVals{allThickVals.id == currID & strcmpi(allThickVals.date,uniqDates{1}) & strcmpi(allThickVals.system,atlasName) & strcmpi(allThickVals.measure,'thickness') & strcmpi(allThickVals.metric,'mean'),'value'};
    currVolume = allThickVals{allThickVals.id == currID & strcmpi(allThickVals.date,uniqDates{1}) & strcmpi(allThickVals.system,atlasName) & strcmpi(allThickVals.measure,'volume') & strcmpi(allThickVals.metric,'numeric'),'value'};
    currLabelsThick = allThickVals{allThickVals.id == currID & strcmpi(allThickVals.date,uniqDates{1}) & strcmpi(allThickVals.system,atlasName) & strcmpi(allThickVals.measure,'thickness') & strcmpi(allThickVals.metric,'mean'),'label'};
    currLabelsVol = allThickVals{allThickVals.id == currID & strcmpi(allThickVals.date,uniqDates{1}) & strcmpi(allThickVals.system,atlasName) & strcmpi(allThickVals.measure,'volume') & strcmpi(allThickVals.metric,'numeric'),'label'};
    
    for L = 1:numLab
        if(sum(currLabelsThick==L)==1)
            if(~isempty(currThickness))
            AllResults.Thickness.Mean(i,L) = currThickness(currLabelsThick==L);
            end
        else
            disp(['Missing Thickness Label ' num2str(L) ' ID ' num2str(currID)])
        end

        if(sum(currLabelsVol==L)==1)
            if(~isempty(currVolume))
            AllResults.Volume.Total(i,L) = currVolume(currLabelsVol==L);
            end
        else
            disp(['Missing Volume Label ' num2str(L) ' ID ' num2str(currID)])
        end
    end
   
end

save(outMatSaveFile,'AllResults');
