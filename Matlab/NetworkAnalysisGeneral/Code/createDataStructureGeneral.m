winBase = fullfile('D:','Min','Dropbox (Personal)','Research','Projects','FTD','NetworkAnalysisGeneral');
macBase = fullfile('/Users','Min','Dropbox (Personal)','Research','Projects','FTD','NetworkAnalysisGeneral');
currBase = macBase;
%Load Laussane Atlas
Lausanne250LabelLUT = readtable(fullfile(currBase,'Data','Lausanne_Scale250.csv'));
Lausanne250numLab=height(Lausanne250LabelLUT);
%labelmeshFile=fullfile(currBase,'Code','NetworkViewer','lausanne250_mesh.mat');
%MeshFile = load(labelmeshFile);
%labelmesh = MeshFile.mesh;
load(fullfile(currBase,'Data','ProcessedResults','Lausanne250SurfaceMeshFromGII.mat'))

%find center of mass, to show the location of nodes
meshlabels = Lausanne250SurfaceMeshFromGII.cdata;
Lausanne250CoM=zeros(height(Lausanne250LabelLUT),3);
for i = 1:height(Lausanne250LabelLUT)
    currLabel = Lausanne250LabelLUT{i,1};
    idx = find(meshlabels==currLabel);
    Lausanne250CoM(i,:)=closestToMean(Lausanne250SurfaceMeshFromGII.giiSurface_Both.vertices(idx,:));
end


Lausanne250CoMDistMatrix = zeros(Lausanne250numLab,Lausanne250numLab);
for i = 1:Lausanne250numLab
    for j = 1:Lausanne250numLab
        Lausanne250CoMDistMatrix(i,j) = norm(Lausanne250CoM(i,:) - Lausanne250CoM(j,:));
    end
end

%Load Sarah's Phase Analysis on Laussane 
tauStaging = fullfile(currBase,'Data','connectivity_autopsy_study','sublists','tc_tau.csv');
tdpStaging = fullfile(currBase,'Data','connectivity_autopsy_study','sublists','tc_tdp.csv');
T_tauStaging = readtable(tauStaging);
T_tdpStaging = readtable(tdpStaging);

T_tauStagingLaus = matchToLaussaune(T_tauStaging,Lausanne250LabelLUT);
T_tauStagingLaus(T_tauStagingLaus==0)=6;

T_tdpStagingLaus = matchToLaussaune(T_tdpStaging,Lausanne250LabelLUT);
T_tdpStagingLaus(T_tdpStagingLaus==0)=6;


SENSAASlangNet=readtable(fullfile(currBase,'Data','sensaas_masked_lausanne250_percentages33.csv'));
SENSAASlangNetMask = zeros(Lausanne250numLab,1);
for i = 1:height(SENSAASlangNet)
    
    currLab = SENSAASlangNet{i,1};
    idx = find(currLab==Lausanne250LabelLUT.Label_ID);
    SENSAASlangNetMask(idx)=1;
end

LausannLeftMask=zeros(Lausanne250numLab,1);
LausannLeftMask(224:448)=1;
LausannRightMask=zeros(Lausanne250numLab,1);
LausannRightMask(1:223)=1;
%Create Data Structure
NetworkDataGeneral=struct;
NetworkDataGeneral.Lausanne250.LUT=Lausanne250LabelLUT;
NetworkDataGeneral.Lausanne250.GII=Lausanne250SurfaceMeshFromGII;
NetworkDataGeneral.Lausanne250.CoM=Lausanne250CoM;
NetworkDataGeneral.Lausanne250.NumLab=Lausanne250numLab;
NetworkDataGeneral.Lausanne250.EucDist=Lausanne250CoMDistMatrix;

NetworkDataGeneral.Lausanne250.masks.TauPhases=T_tauStagingLaus;
NetworkDataGeneral.Lausanne250.masks.TDPPhases=T_tdpStagingLaus;
NetworkDataGeneral.Lausanne250.masks.SENSAASLangNet=SENSAASlangNetMask;
NetworkDataGeneral.Lausanne250.masks.LeftHemi=LausannLeftMask;
NetworkDataGeneral.Lausanne250.masks.RightHemi=LausannRightMask;


%load Schaefer400 Atlas
Schaefer400x7LabelLUT=readtable(fullfile(currBase,'Data','schaefer400x7_lut.csv'));
load(fullfile(currBase,'Data','ProcessedResults','schaefer400x7SurfaceMeshFromGII.mat'))

%find center of mass, to show the location of nodes
meshlabels = schaefer400x7SurfaceMeshFromGII.cdata;
Schaefer400x7CoM=zeros(height(Schaefer400x7LabelLUT),3);
for i = 1:height(Schaefer400x7LabelLUT)
    i
    currLabel = Schaefer400x7LabelLUT.Label_ID(i);
    idx = find(meshlabels==currLabel);
    Schaefer400x7CoM(i,:)=closestToMean(schaefer400x7SurfaceMeshFromGII.giiSurface_Both.vertices(idx,:));
end
Schaefer400x7numLab=height(Schaefer400x7LabelLUT);
Schaefer400x7CoMDistMatrix = zeros(Schaefer400x7numLab,Schaefer400x7numLab);

for i = 1:Schaefer400x7numLab
    for j = 1:Schaefer400x7numLab
        Schaefer400x7CoMDistMatrix(i,j) = norm(Schaefer400x7CoM(i,:) - Schaefer400x7CoM(j,:));
    end
end

NetworkDataGeneral.Schaefer400x7.LUT=Schaefer400x7LabelLUT;
NetworkDataGeneral.Schaefer400x7.GII=schaefer400x7SurfaceMeshFromGII;
NetworkDataGeneral.Schaefer400x7.CoM=Schaefer400x7CoM;
NetworkDataGeneral.Schaefer400x7.NumLab=Schaefer400x7numLab;
NetworkDataGeneral.Schaefer400x7.EucDist=Schaefer400x7CoMDistMatrix;

save(fullfile(currBase,'FTDGeneralData_20221114.mat'),'NetworkDataGeneral')