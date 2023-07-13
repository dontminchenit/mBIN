
%AllCasesInd = [ GroupInd{2}; GroupInd{3}]
%AllCaseMtx = [AllResults.FAGraph.WeightMtxNormal(GroupInd{2}); AllResults.FAGraph.WeightMtxNormal(GroupInd{3})]

baseData=SubVolsWscore;
%baseData(baseData>0)=0;

data1 = baseData(strcmpi(subDemo.Group,'tau'),:);
data2 = baseData(strcmpi(subDemo.Group,'tdp'),:);
dataAll= [data1; data2];

str1='tau'
str2='TDP'

N=size(baseData,2);
allRegR = nan(1,N);
allRegL = nan(1,N);
for i = 1:N

%        if(i ~= j)

        %figure(2)
        %clf
        %scatter(a(:),b(:))
        %pause
        
        currData = baseData(:,i);
        %[p,tbl] = anova1(currData,subDemo.Group,'off')
        [h p]=ttest2(currData(strcmpi(subDemo.Group,'tau')),currData(strcmpi(subDemo.Group,'tdp')),Tail='right');
        allRegR(i)=p;

        [h p]=ttest2(currData(strcmpi(subDemo.Group,'tau')),currData(strcmpi(subDemo.Group,'tdp')),Tail='left');
        allRegL(i)=p;

        %[RHO,PVAL] = corr(a(:),b(:));
  %      RHO = icc(a,b);
   %     allCorr(i,j)=RHO;
 %       end
end

%figure(1)
%clf
%imagesc(allCorr)

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
pointSize=15;
setViews=1:6;
ViewLoc=1:6;
fillFlag=3;
surfDisp=[]
cmapTypeForce=3;


dataR=allRegR<.05;
dataL=allRegL<.05;
data = double(dataR) - double(dataL);
plotBrain3Points(CoM,GII,data,alphaThresh,panelAll,1,cmapTypeForce,pointSize,surfDisp,setViews,ViewLoc)
print(H2,fullfile(saveDirBase,'wscore_sigdiff.tif'),'-dpng','-r400');