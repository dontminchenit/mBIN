winBase = fullfile('D:','Min','Dropbox (Personal)','Research','Projects','FTD');
macBase = fullfile('/Users','Min','Dropbox (Personal)','Research','Projects','FTD');
currBaseGen=winBase;
currBase = fullfile(currBaseGen,'LBD');

load(fullfile(currBaseGen,'NetworkAnalysisGeneral','FTDGeneralData_20221114.mat'))

baseSaveDir=fullfile(currBase,'Results_12_9_2022');
mkdir(baseSaveDir)
demoT = readtable(fullfile(currBase,'Data','MRI Schaeffer Demographics classification.xlsx'));
pathT = readtable(fullfile(currBase,'Data','LBD_AO_aSYN.csv'));

LBD_yesADIndx=strcmpi(pathT.Group,'AD+LBD') & ~strcmpi(pathT.CDx,'AD');
LBD_noADIndx=strcmpi(pathT.Group,'LBD');

pathDataGM=log(.01*pathT{:,2:7}+0.00015);
pathNamesRaw = pathT.Properties.VariableNames(2:7);
SN= length(pathNamesRaw);

thicknessCtrlraw = readtable(fullfile(currBase,'Data','ctrl_4007.csv'));
ctrlIDs = unique(thicknessCtrlraw.id);
%thicknessTauTDP = readtable(fullfile(currBase,'Data','tautdp_4007.csv'));

thickIDs=demoT.INDDID;
thicknessAllraw=readtable(fullfile(currBase,'Data','ftdc_dlb_ctx_volume.csv'));
volAllraw = readtable(fullfile(currBase,'Data','ftdc_dlb_tian_s1.csv'));



%ctrlIDs = unique(thicknessCtrl.id);

AllIDs = pathT.INDDID;
%CohortsIndexs = [ones(length(ctrlIDs),1);2*ones(length(TauIDs),1);3*ones(length(TDPIDS),1)];
numLab=400;

mkdir(fullfile(currBase,'proccessedResults'));
thicknessSaveFile=fullfile(currBase,'proccessedResults','thicknessVales.mat');
subVolSaveFile=fullfile(currBase,'proccessedResults','subVolSaveFile.mat');
thicknessSaveFileCtrl=fullfile(currBase,'proccessedResults','thicknessValesCtrl.mat');
subVolSaveFileCtrl=fullfile(currBase,'proccessedResults','subVolSaveFileCtrl.mat');

[AllResults, numTimePoints, YrFromBaseline, baselineYear] = LoadNetworkDataByID_Longitudinal(thickIDs,thicknessAllraw,thicknessSaveFile,'schaefer400x7v1',400);
%[SubCorticalVol] = LoadNetworkDataByID2(thickIDs,volAllraw,subVolSaveFile,'Tian_Subcortex_S1_3T',16);

%[CtrlResults] = LoadNetworkDataByID(ctrlIDs,thicknessCtrlraw,thicknessSaveFileCtrl,'schaefer400x7',400);
%[SubCorticalVolCtrl] = LoadNetworkDataByID(ctrlIDs,subcorticalVolCtrlraw,subVolSaveFileCtrl,'Tiavn_Subcortex_S1_3T',16);

%load(thicknessSaveFileCtrl)


%load Path data
atlasBase=fullfile(currBaseGen,'InvivoPath','Data');
pathLUT = readtable(fullfile(atlasBase,'schaefer_path_20210719_20220328.csv'));

AtlasToPathLUT = readtable(fullfile(currBase,'Data','PathToAtlasLUT_10_7_2022.xlsx'));

[pathCoMunordered] = findPathCoM(pathLUT,AtlasToPathLUT,NetworkDataGeneral.Schaefer400x7.CoM);
pathCoM=nan(SN,3,2);
for s=1:SN

    idx = find(strcmpi(AtlasToPathLUT.PathSpreadSheetNames,pathNamesRaw(s)));
    pathCoM(s,:,:)=pathCoMunordered(idx,:,:);

end

%[pathDataGM subMatchedPathDataGM  TotalSubsMatchedGM TotalRegionMatchedGM TotalSumbRegionMatchedGM] = matchToPath(AllIDs,pathT,pathLUT,LabelLUT,AtlasToPathLUT,0); 
%[pathDataWM subMatchedPathDataWM TotalSubsMatchedWM TotalRegionMatchedWM TotalSumbRegionMatchedWM] = matchToPath(AllIDs,pathT,pathLUT,LabelLUT,AtlasToPathLUT,1); 


