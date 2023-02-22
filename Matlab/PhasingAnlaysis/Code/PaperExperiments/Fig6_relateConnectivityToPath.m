%Load Data
%LoadDataSarah

PathDataIn= fullfile(currBase,'PhaseData','FTLD Library 05-06-21 with DTI.xlsx');
PathToLaussLUT=fullfile(currBase,'PhaseData','copath_labels_lausanne250_20191104_V2.lut.csv');
PathInfo = fullfile(currBase,'PhaseData','dti_subs_5subRemoved.csv');


PathInfo_T=readtable(PathInfo);%List of subjects in cohort
PathDataIn_T=readtable(PathDataIn);%Path data
PathToLaussLUT=readtable(PathToLaussLUT);%LUT for path region to lussane atlas

LabelN = height(LabelLUT); %#number of labels in atlas
SubN = height(PathInfo_T);%number of subjects
GMPath = nan(SubN,LabelN);%save each subjects GM Path data for each lassuane label
WMPath = nan(SubN,LabelN);%save each subjects WM Path data for each lassuane label

GMPathLabel = nan(SubN,LabelN);%save the Path region label for each lassuane label
WMPathLabel = nan(SubN,LabelN);%save the Path region label for each lassuane label

GMCount=0;
WMCount=0;

%parse region names from path info
for s = 1:SubN
    
    indSub=find(PathDataIn_T.INDDID==PathInfo_T.INDDID(s));%find current subject in path data
    
    %cycle through all laussane labels
    for i = 1:LabelN
        
       currRegionName = PathToLaussLUT.Label_AbbrevName{i};%current name of label
       if (~isempty(currRegionName))
           %currRegionLR=PathToLaussLUT.Label_Hemisphere{i}(1);

           currReg = strrep(currRegionName,' ','_');%modify label name to match path region name
           
           srchPatGM=[currReg '_GM'];
           srchPatWM=[currReg '_WM'];

           %find region name in path data
           indGM = find(contains(lower(PathDataIn_T.Properties.VariableNames),lower(srchPatGM)));
           indWM = find(contains(lower(PathDataIn_T.Properties.VariableNames),lower(srchPatWM)));
           
      
           %pull GM/WM value for this region for the subject
           pathValGM = PathDataIn_T{indSub,indGM};
           pathValWM = PathDataIn_T{indSub,indWM};
           
           %save it to the data structure
           if(~isempty(pathValGM))
               srchPatGM
               if(iscell(pathValGM))
                 GMPath(s,i)= str2double(pathValGM{1});
               else
                 GMPath(s,i)= pathValGM;
               end
               GMPathLabel(s,i)=indGM;
           end
           
           if(~isempty(pathValWM))
               if(iscell(pathValWM))
                 WMPath(s,i)= str2double(pathValWM{1});
               else
                 WMPath(s,i)= pathValWM;
               end
               WMPathLabel(s,i)=indWM; 
           end
       end
       
    end
end
