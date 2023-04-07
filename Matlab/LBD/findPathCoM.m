function [pathCoM, pathToAtlasIndex] = findPathCoM(pathLUT,AtlasToPathLUT,AtlasCoM)



%NS=size(thicknessVals,1);
NL=height(AtlasToPathLUT);

%matchedThicknessData=nan(NS,NL,2);
%matchedPathData=nan(NS,NL,2);
pathCoM=zeros(NL,3,2);
pathToAtlasIndex = cell(NL,2); 
hemi = {'L','R'}

        for l = 1:NL
            for r=1:2
            AOName = AtlasToPathLUT.Matched{l};
            
            pathLabelIdx = find(matches(pathLUT.PathName,AOName) & matches(pathLUT.Label_Hemi,hemi{r}));
            if(~isempty(pathLabelIdx))
                atlasIdx = pathLUT.Label_ID400x7(pathLabelIdx);
             
                pathCoM(l,:,r)=mean(AtlasCoM(atlasIdx,:),1,'omitnan');
                pathToAtlasIndex{l,r}=atlasIdx;
                    
                %giiSurf=GII.GII.giiSurface_Both;
                %valDisp=100*(ismember(GII.GII.cdata,atlasIdx));
                %figure(1)
                %clf
                %t=trisurf(giiSurf.faces,giiSurf.vertices(:,1),giiSurf.vertices(:,2),giiSurf.vertices(:,3),valDisp(:),'EdgeColor','none');
                %alpha(t,.05)
 %               for s = 1:NS
                    %matchedThicknessData(s,l,r) = mean(thicknessVals(s,atlasIdx),2,'omitnan');
                    %matchedPathData(s,l,r) = mean(pathVals(s,atlasIdx),2,'omitnan');
                    
                    %hold on;
                    %scatter3(meanloc(1))
                    %scatter3(AtlasCoM(atlasIdx,1),AtlasCoM(atlasIdx,2),AtlasCoM(atlasIdx,3))
                    %scatter3(pathCoM(l,1,r),pathCoM(l,2,r),pathCoM(l,3,r),50)
                    %hold off;
  %              end
            end
            end
        end
    
