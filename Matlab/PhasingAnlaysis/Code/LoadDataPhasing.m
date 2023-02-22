%winBase = fullfile('D:','Min','Dropbox (Personal)','Research','Projects','FTD','PhasingAnalysis');
%macBase = fullfile('/Users','Min','Dropbox (Personal)','Research','Projects','FTD','PhasingAnalysis');
%currBase = winBase;

currBase=dataDir;

inInfo = fullfile(currBase,'PhaseData','subj95_baseline.xlsx');
T = readtable(inInfo);

inInfo = fullfile(currBase,'PhaseData','dti_subs_5subRemoved.csv');
TPhases = readtable(inInfo);

TauIND = strcmpi(TPhases.Group,'tau');
TDPIND = strcmpi(TPhases.Group,'tdp');

LabelLUT = readtable(fullfile(currBase,'PhaseData','Lausanne_Scale250.csv'));
numLab=height(LabelLUT);

labelmeshFile=fullfile(currBase,'PhaseData','lausanne250_mesh.mat');
MeshFile = load(labelmeshFile);
labelmesh = MeshFile.mesh;
load(fullfile(currBase,'PhaseData','Lausanne250SurfaceMeshFromGII.mat'))
%meshlabelsLauss = labelmesh.attributes.attribute_array;

%saveDirBase = fullfile(currBase,'RevisionResults_9_14_2022');
saveDirBase = fullfile(outputDir,'RevisionResults_9_14_2022');
mkdir(saveDirBase)

savedResults = fullfile(currBase,'PhaseData','ProcessedData','allMeasures_AllGraphs_11_1_2020.mat');
PPAAllResults = load(savedResults);

savedResults = fullfile(currBase,'PhaseData','ProcessedData','allMeasures_AllGraphs_Sarah_7_29_2020.mat');
load(savedResults)

HCIND_Only =strcmpi(T.StudyGroup,'eld');

S = AllResults.FAGraph;
T = PPAAllResults.AllResults.FAGraph;
Z = cell2struct(cellfun(@vertcat,struct2cell(S),struct2cell(T),'uni',0),fieldnames(S),1);

S1 = AllResults.SCGraph;
T1 = PPAAllResults.AllResults.SCGraph;
Z1 = cell2struct(cellfun(@vertcat,struct2cell(S1),struct2cell(T1),'uni',0),fieldnames(S1),1);



AllResults = struct;
AllResults.FAGraph=Z;
AllResults.SCGraph=Z1;


HCIND = logical([zeros(length(TauIND),1); HCIND_Only]);
TauIND =logical([TauIND; zeros(length(HCIND_Only),1)]);
TDPIND =logical([TDPIND; zeros(length(HCIND_Only),1)]);

graphNames={'SCGraph', 'FAGraph', 'TractLenGraph', 'Thickness'};
measureNames={'BetweenCen','ClusterCoeff','LocalEff','EigenCen','Degree','Mean','Median'};

% LaussCoM=zeros(height(LabelLUT),3);
% for i = 1:height(LabelLUT)
%     currLabel = LabelLUT{i,1};
%     idx = find(meshlabelsLauss==currLabel);
%     LaussCoM(i,:)=mean(labelmesh.points(idx,:),1);
% end


%find center of mass, to show the location of nodes
meshlabels = Lausanne250SurfaceMeshFromGII.cdata;
CoM=zeros(height(LabelLUT),3);
for i = 1:height(LabelLUT)
    currLabel = LabelLUT{i,1};
    idx = find(meshlabels==currLabel);
    CoM(i,:)=closestToMean(Lausanne250SurfaceMeshFromGII.giiSurface_Both.vertices(idx,:));
end

CoMDistMatrix = zeros(numLab,numLab);

for i = 1:numLab
    for j = 1:numLab
        CoMDistMatrix(i,j) = norm(CoM(i,:) - CoM(j,:));
    end
end


saveMeans = 1;
saveComp = 1;

GroupNames = {'HC','TAU','TDP'};
GroupInd = cell(1,3);
GroupInd{1}=HCIND;
GroupInd{2}=TauIND;
GroupInd{3}=TDPIND;

%define hubs and atrophied network nodes.
defineHubs

tauStaging = fullfile(currBase,'PhaseData','connectivity_autopsy_study','connectivity_autopsy_study','sublists','tc_tau.csv');
tdpStaging = fullfile(currBase,'PhaseData','connectivity_autopsy_study','connectivity_autopsy_study','sublists','tc_tdp.csv');
T_tauStaging = readtable(tauStaging);
T_tdpStaging = readtable(tdpStaging);

T_tauStagingLaus = matchToLaussaune(T_tauStaging,LabelLUT);
T_tauStagingLaus(T_tauStagingLaus==0)=6;

T_tdpStagingLaus = matchToLaussaune(T_tdpStaging,LabelLUT);
T_tdpStagingLaus(T_tdpStagingLaus==0)=6;
