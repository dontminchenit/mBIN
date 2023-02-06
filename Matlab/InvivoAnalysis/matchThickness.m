function [matchedThicknessData, matchedPathData, pathCoM] = matchThickness(thicknessVals,pathVals,pathLUT,AtlasToPathLUT,AtlasCoM)



NS=size(thicknessVals,1);
NL=height(AtlasToPathLUT);

matchedThicknessData=nan(NS,NL,2);
matchedPathData=nan(NS,NL,2);
pathCoM=zeros(NL,3,2);

hemi = {'L','R'}

        for l = 1:NL
            for r=1:2
            AOName = AtlasToPathLUT.Matched{l};
            
            pathLabelIdx = find(matches(pathLUT.PathName,AOName) & matches(pathLUT.Label_Hemi,hemi{r}));
            if(~isempty(pathLabelIdx))
                atlasIdx = pathLUT.Label_ID400x7(pathLabelIdx);
                for s = 1:NS
                    matchedThicknessData(s,l,r) = mean(thicknessVals(s,atlasIdx),2,'omitnan');
                    matchedPathData(s,l,r) = mean(pathVals(s,atlasIdx),2,'omitnan');
                    
                    pathCoM(l,:,r)=mean(AtlasCoM(atlasIdx,:),1,'omitnan');
                end

            end
            end
        end
    
