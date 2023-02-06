function [matchedPathDataAtlas matchedPathData TotalSubsMatched TotalRegionMatched TotalSumbRegionMatched] = matchToPath(AllIDs,pathT,pathLUT,LabelLUT,AtlasToPathLUT,GM0WM1)



NS=size(AllIDs,1);
NL=size(LabelLUT,1);
NPL=height(AtlasToPathLUT);


matchedPathDataAtlas=nan(NS,NL);
matchedPathData=nan(NS,NPL,2);

if (GM0WM1==0)
    Tissuetype='GM';
elseif(GM0WM1==1)
    Tissuetype='WM';
end
    
TotalSubsMatched=0;
TotalRegionMatched=0;
TotalSumbRegionMatched=zeros(NS,1);

for s = 1:NS
    currID = AllIDs(s);
    
    pathDataSubIndex = find(pathT.INDDID==currID);
    if(~isempty(pathDataSubIndex))

        %match subject to path labels
        
        for p = 1:NPL
            
            AOName=AtlasToPathLUT.PathSpreadSheetNames{p};
            pathDataName = ['L_' Tissuetype '_' AOName];
            pathDataRegionIdx = find(matches(pathT.Properties.VariableNames,pathDataName,'IgnoreCase',true));
            matchedPathData(s,p,1)=pathT{pathDataSubIndex,pathDataRegionIdx};
            
            pathDataName = ['R_' Tissuetype '_' AOName];
            pathDataRegionIdx = find(matches(pathT.Properties.VariableNames,pathDataName,'IgnoreCase',true));
            matchedPathData(s,p,2)=pathT{pathDataSubIndex,pathDataRegionIdx};
        end
        
        %match to Atlas
        
        TotalSubsMatched=TotalSubsMatched+1;
        for l = 1:NL

            currLabelName=LabelLUT.Label_Name{l};
            pathLabelIdx = find(matches(pathLUT.Label_Name400x7,currLabelName));
            pathLabelName = pathLUT.PathName{pathLabelIdx};
            pathHemi = pathLUT.Label_Hemi{pathLabelIdx};

            AONameIdx = find(matches(AtlasToPathLUT.Matched,pathLabelName));
            if(~isempty(AONameIdx))
                AOName = AtlasToPathLUT.PathSpreadSheetNames{AONameIdx};
                %AOName
            else
                continue
            end
            pathDataName = [pathHemi '_' Tissuetype '_' AOName];
            %pathDataName = ['Core_' TissueType '_' pathHemi '_' 

            pathDataRegionIdx = find(matches(pathT.Properties.VariableNames,pathDataName,'IgnoreCase',true));

            %pathT.Properties.VariableNames(pathDataRegionIdx)
            if(~isempty(pathDataRegionIdx))
                %pathT{pathDataSubIndex,pathDataRegionIdx}
                matchedPathDataAtlas(s,l) = pathT{pathDataSubIndex,pathDataRegionIdx};
                TotalRegionMatched=TotalRegionMatched+1;
                TotalSumbRegionMatched(s)=TotalSumbRegionMatched(s)+1;
            end
        end
    end
end