function StagingLaus = matchToLaussaune(Staging,LabelLUT)


StagingLaus = zeros(height(LabelLUT),1);

for i = 1:height(LabelLUT)
   
   StagingLaus(i) = Staging.Stage(strcmpi(LabelLUT.Label_Name{i},Staging.Region));
    
end

StagingLaus(isnan(StagingLaus))=0;