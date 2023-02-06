winBase = fullfile('D:','Min','Dropbox (Personal)','Research','Projects','FTD','InvivoPath');
macBase = fullfile('/Users','Min','Dropbox (Personal)','Research','Projects','FTD','InvivoPath');
currBase = winBase;

baseSaveDir=fullfile(currBase,'Results_9_16_2022');
mkdir(baseSaveDir)
demoT = readtable(fullfile(currBase,'Data','AUTPOSY MRI DEMOGRAPHICS(2022.03.25 13.28).xlsx'));

TauIDs=demoT.inddid(demoT.BVFTD==1 & demoT.TDP1Tau2==2);
TDPIDS=demoT.inddid(demoT.BVFTD==1 & demoT.TDP1Tau2==1);

thicknessCtrl = readtable(fullfile(currBase,'Data','ctrl_4007.csv'));
thicknessTauTDP = readtable(fullfile(currBase,'Data','tautdp_4007.csv'));

thicknessAll=[thicknessCtrl; thicknessTauTDP];

ctrlIDs = unique(thicknessCtrl.id);

AllIDs = [ctrlIDs; TauIDs; TDPIDS];
CohortsIndexs = [ones(length(ctrlIDs),1);2*ones(length(TauIDs),1);3*ones(length(TDPIDS),1)];
numLab=400;


thicknessSaveFile=fullfile(currBase,'proccessedResults','thicknessVales.mat');

%[AllResults] = LoadNetworkDataByID(AllIDs,thicknessAll,thicknessSaveFile)
load(thicknessSaveFile)

LabelLUT=readtable(fullfile(currBase,'Data','schaefer400x7_lut.csv'));

load(fullfile(currBase,'Data','schaefer400x7SurfaceMeshFromGII.mat'))

%load Path data
pathT = readtable(fullfile(currBase,'Data','Final_Network_Path_Sheet_5_7_21_wCore.xlsx'));

pathLUT = readtable(fullfile(currBase,'Data','schaefer_path_20210719_20220328.csv'));

AtlasToPathLUT = readtable(fullfile(currBase,'Data','PathToAtlasLUT_6_10_2022.xlsx'));



[pathDataGM subMatchedPathDataGM  TotalSubsMatchedGM TotalRegionMatchedGM TotalSumbRegionMatchedGM] = matchToPath(AllIDs,pathT,pathLUT,LabelLUT,AtlasToPathLUT,0); 
[pathDataWM subMatchedPathDataWM TotalSubsMatchedWM TotalRegionMatchedWM TotalSumbRegionMatchedWM] = matchToPath(AllIDs,pathT,pathLUT,LabelLUT,AtlasToPathLUT,1); 
[matchedThicknessData matchedPathDataGM pathCoM] = matchThickness(AllResults.Thickness.Mean,pathDataGM,pathLUT,AtlasToPathLUT,CoM)
