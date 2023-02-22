H2=figure(2)
clf
panelAll = panel();
panelAll.pack(2,6);
panelAll.de.margin = 1;
panelAll.marginbottom = 0;
panelAll.margin = [1 1 1 1];


GII=NetworkDataGeneral.Lausanne250.GII;
CoM=NetworkDataGeneral.Lausanne250.CoM;
numLab=NetworkDataGeneral.Lausanne250.NumLab;

alphaThresh=[];
pointSize=50;
setViews=1:6;
ViewLoc=1:6;
fillFlag=3;
surfDisp=[]
cmapTypeForce=2;

data = mean(SubVolsWscore(strcmpi(subDemo.Group,'tau'),:))/2;
data(data>0)=0;
plotBrain3Points(CoM,GII,data,alphaThresh,panelAll,1,cmapTypeForce,pointSize,surfDisp,setViews,ViewLoc)

data = mean(SubVolsWscore(strcmpi(subDemo.Group,'tdp'),:))/2;
data(data>0)=0;
plotBrain3Points(CoM,GII,data,alphaThresh,panelAll,2,cmapTypeForce,pointSize,surfDisp,setViews,ViewLoc)

saveDirBase='D:\Min\Dropbox (Personal)\Research\Projects\FTD\PhasingAnalysis\revResults';
mkdir(saveDirBase);
print(H2,fullfile(saveDirBase,'wscore_TauTop_TDPBot.tif'),'-dpng','-r400');
