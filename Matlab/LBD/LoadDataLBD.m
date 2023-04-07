%winBase = fullfile('D:','Min','Dropbox (Personal)','Research','Projects','FTD');
%macBase = fullfile('/Users','Min','Dropbox (Personal)','Research','Projects','FTD');
%currBase = fullfile(winBase,'LBD');
%load(fullfile(winBase,'NetworkAnalysisGeneral','FTDGeneralData_20221114.mat'))
%baseSaveDir=fullfile(currBase,'Results_1_6_2023');

currBase=dataDir;
baseSaveDir=outputDir;

mkdir(baseSaveDir)

demoT = readtable(fullfile(currBase,'LBDData','MRI Schaeffer Demographics classification.xlsx'));
pathT = readtable(fullfile(currBase,'LBDData','LBD_AO_aSYN.csv'));
measuresT=readtable(fullfile(currBase,'LBDData','LBD neuropsych annualized change no formulas.xlsx'));

%pathT = readtable(fullfile(currBase,'Data','LBD %AO Tau LBD 11152022.csv')); %tau %AO


LBD_yesADIndx=strcmpi(pathT.Group,'AD+LBD') & ~strcmpi(pathT.CDx,'AD');
LBD_noADIndx=strcmpi(pathT.Group,'LBD');

pathDataGM=log(.01*pathT{:,2:7}+0.00015);
pathNamesRaw = pathT.Properties.VariableNames(2:7);
SN= length(pathNamesRaw);

thicknessCtrlraw = readtable(fullfile(currBase,'LBDData','ctrl_4007.csv'));
ctrlIDs = unique(thicknessCtrlraw.id);
%thicknessTauTDP = readtable(fullfile(currBase,'Data','tautdp_4007.csv'));

thickIDs=demoT.INDDID;

thicknessAllraw=readtable(fullfile(currBase,'LBDData','ftdc_dlb_ctx_volume.csv'));
volAllraw = readtable(fullfile(currBase,'LBDData','ftdc_dlb_tian_s1.csv'));

[idx idx2] = ismember(thickIDs, measuresT.INDDID);
measuresT_matched=measuresT(idx2(idx),:);

%ctrlIDs = unique(thicknessCtrl.id);

AllIDs = pathT.INDDID;
%CohortsIndexs = [ones(length(ctrlIDs),1);2*ones(length(TauIDs),1);3*ones(length(TDPIDS),1)];
numLab=400;

mkdir(fullfile(currBase,'proccessedResultsLBD'));

thicknessSaveFile=fullfile(currBase,'proccessedResultsLBD','thicknessVales.mat');
subVolSaveFile=fullfile(currBase,'proccessedResultsLBD','subVolSaveFile.mat');
thicknessSaveFileCtrl=fullfile(currBase,'proccessedResultsLBD','thicknessValesCtrl.mat');
subVolSaveFileCtrl=fullfile(currBase,'proccessedResultsLBD','subVolSaveFileCtrl.mat');

[AllResults] = LoadNetworkDataByID(thickIDs,thicknessAllraw,thicknessSaveFile,'schaefer400x7v1',400);
[SubCorticalVol] = LoadNetworkDataByID2(thickIDs,volAllraw,subVolSaveFile,'Tian_Subcortex_S1_3T',16);
[CtrlResults] = LoadNetworkDataByID(ctrlIDs,thicknessCtrlraw,thicknessSaveFileCtrl,'schaefer400x7',400);

%[SubCorticalVolCtrl] = LoadNetworkDataByID(ctrlIDs,subcorticalVolCtrlraw,subVolSaveFileCtrl,'Tiavn_Subcortex_S1_3T',16);
%load(thicknessSaveFileCtrl)


%load Path data
%atlasBase=fullfile(winBase,'InvivoPath','Data');
atlasBase=currBase;  

pathLUT = readtable(fullfile(atlasBase,'schaefer_path_20210719_20220328.csv'));
AtlasToPathLUT = readtable(fullfile(currBase,'LBDData','PathToAtlasLUT_10_7_2022.xlsx'));

[pathCoMunordered, pathToAtlasIndexunordered] = findPathCoM(pathLUT,AtlasToPathLUT,NetworkDataGeneral.Schaefer400x7.CoM);
pathCoM=nan(SN,3,2);
pathToAtlasIndex = cell(SN,2);
for s=1:SN

    idx = find(strcmpi(AtlasToPathLUT.PathSpreadSheetNames,pathNamesRaw(s)));
    pathCoM(s,:,:)=pathCoMunordered(idx,:,:);
    pathToAtlasIndex(s,:)=pathToAtlasIndexunordered(idx,:);
end

N=size(AllResults.Thickness.Mean,1);
PN = size(pathDataGM,2);

AllthicknessAtPath=nan(N,PN,2);
CtrlthicknessAtPath=nan(N,PN,2);


for n = 1:N
    for p = 1:PN
       for r=1:2
           curIdx=pathToAtlasIndex{p,r};
            AllthicknessAtPath(n,p,r)= mean(AllResults.Thickness.Mean(n,curIdx));
            CtrlthicknessAtPath(n,p,r)= mean(CtrlResults.Thickness.Mean(n,curIdx));
       end
    end
end


%[pathDataGM subMatchedPathDataGM  TotalSubsMatchedGM TotalRegionMatchedGM TotalSumbRegionMatchedGM] = matchToPath(AllIDs,pathT,pathLUT,LabelLUT,AtlasToPathLUT,0); 
%[pathDataWM subMatchedPathDataWM TotalSubsMatchedWM TotalRegionMatchedWM TotalSumbRegionMatchedWM] = matchToPath(AllIDs,pathT,pathLUT,LabelLUT,AtlasToPathLUT,1); 


