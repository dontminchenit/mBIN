macBase='/Users/min/Dropbox (Personal)/Research/Projects/FTD/PhasingAnalysis/';
winBase='D:\Min\Dropbox (Personal)\Research\Projects\FTD\PhasingAnalysis\';
currBase = macBase;
%loaddata
ctrlVols = readtable(fullfile(currBase,'Data','w-score','ctrl_vl250.csv'));
ctrldemo = readtable(fullfile(currBase,'Data','w-score','ctrl_demos.csv'));
subDemo = readtable(fullfile(currBase,'Data','w-score','mc_revisions_bvols_mmse_cdr.csv'));

NCtrl = height(ctrldemo);
CtrlVolsData = nan(NCtrl,numLab);
CtrlTotalVol = ctrldemo.BVOL;
CtrlAge = ctrldemo.AgeatMRI;
CtrlSex = strcmpi(ctrldemo.Sex,'Female');

for i = 1:NCtrl
    currLabels = ctrlVols.label(ctrlVols.id==ctrldemo.INDDID(i));
    currVols =  ctrlVols.value(ctrlVols.id==ctrldemo.INDDID(i));
    CtrlVolsData(i,:)=currVols;
end

NSub = height(subDemo);
SubVolsData = nan(NSub,numLab);
SubVolsThick = nan(NSub,numLab);
SubTotalVol = subDemo.BVOL;
SubAge = subDemo.AgeatMRI;
SubSex = strcmpi(subDemo.Sex,'Female');

dataDir=fullfile(currBase,'Data','pathcsvs20210125','csvsSent');
for i = 1:NSub
    currID = subDemo.INDDID(i);

    matchedFile = dir(fullfile(dataDir,[num2str(currID) '*lausanne.csv']));
    currAllData =readtable(fullfile(dataDir,matchedFile.name));
    volIdx=strcmp(currAllData.measure,'volume') & strcmp(currAllData.system,'lausanne250');
    currLabesl = currAllData.label(volIdx);
    currVols = currAllData.value(volIdx);
    SubVolsData(i,:)=currVols;

    currThickIdx =strcmp(currAllData.metric,'mean') & strcmp(currAllData.measure,'thickness') & strcmp(currAllData.system,'lausanne250');
    currThick = currAllData.value(currThickIdx);
    SubVolsThick(i,:)=currThick;
end
%calculate control coefficients

HCCoeff=nan(numLab,3); %1-Intercept, 2-Age, 3-TotalVol
HCSigma=nan(numLab,1); %Sigma of Residual
HCModels = cell(numLab,1); %saves the regression model

%X=[CtrlAge CtrlTotalVol];
X=[CtrlAge CtrlTotalVol CtrlSex];
for i = 1:numLab

    y=CtrlVolsData(:,i);
    mdl = fitlm(X,y);

    HCModels{i}=mdl;
    HCCoeff(i,1)=mdl.Coefficients{1,1};
    HCCoeff(i,2)=mdl.Coefficients{2,1};
    HCCoeff(i,3)=mdl.Coefficients{3,1};
    HCSigma(i) = std(mdl.Residuals.Raw);


end

%calculate w-score using HC coeff
%X=[SubAge SubTotalVol];
X=[SubAge SubTotalVol SubSex];
SubVolsWscore  = nan(NSub,numLab);
for k = 1:numLab
    currMdl = HCModels{k};
    labVolPredict = predict(currMdl,X);
    labCurrVols=SubVolsData(:,k);
    SubVolsWscore(:,k) = (labCurrVols-labVolPredict)./HCSigma(k);
end
