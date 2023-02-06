function [Path2LaussIdx, PathCoM, OffsetV1, OffsetV2] = PathLabels2Laussane(PathLabelName,LaussLabelLUT,LaussCoM)

NP=length(PathLabelName);
NL=height(LaussLabelLUT);

Path2LaussIdx=nan(NP,NL);

PathCoM=zeros(NP,3);
offsetVal = -200;
GM0WM1=false(NP,1);

for p = 1:NP
    
    currPName = PathLabelName{p};
    
    if(contains(currPName,'_WM_'))
        GM0WM1(p)=true;
    end
    
    currPName = strrep(currPName,'KeyCore','L');
    currPName = strrep(currPName,'Core','L');
    currPName = strrep(currPName,'_GM_',' ');
    compName = strrep(currPName,'_WM_',' ');
    c = 0;
    for l = 1:NL
        
        currLName = LaussLabelLUT.Label_AbbrevName{l};
        match=strcmpi(currLName,compName);
        if(match)
           c=c+1;
            Path2LaussIdx(p,c) = l;
        end
    end
    
    
    if (c==0)
       disp([currPName ' cannot find match']);
    else
        matchedIdx = Path2LaussIdx(p,~isnan(Path2LaussIdx(p,:)));
        PathCoM(p,:) = mean(LaussCoM(matchedIdx,:),1);
    end
end

%set offset for WM values
OffsetV1=zeros(NP,3);
OffsetV2=zeros(NP,3);

OffsetV1(GM0WM1,3) = offsetVal;
OffsetV2(GM0WM1,2) = offsetVal;

end