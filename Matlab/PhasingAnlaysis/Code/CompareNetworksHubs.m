%function plotBrain(mesh,data,LabelLUT,CRange,saveDir,saveName,alpha)
%saveDirBase = 'D:\Research\ResearchResults\PPA\NetworkResultsZeroThresh';
mkdir(saveDirBase);
savedResults = 'D:\Research\ResearchResults\PPA\allMeasures_AllGraphs_11_1_2020.mat';
%load(savedResults)

LabelLUT = readtable('../NetworkViewer/Lausanne_Scale250.csv');
loadMesh = load('../NetworkViewer/lausanne250_mesh.mat');
Mesh = loadMesh.mesh;
inInfo = '../NetworkViewer/subj95_baseline.xlsx';
T = readtable(inInfo);

HCIND =strcmpi(T.StudyGroup,'eld');
svIND =strcmpi(T.StudyGroup,'svPPA');
naIND =strcmpi(T.StudyGroup,'naPPA');

meshlabels = Mesh.attributes.attribute_array;
%find center of mass, to show the location of nodes
CoM=zeros(height(LabelLUT),3);
for i = 1:height(LabelLUT)
    currLabel = LabelLUT{i,1};
    idx = find(meshlabels==currLabel);
    CoM(i,:)=mean(Mesh.points(idx,:),1);
end

%groups we'll be comparing against
GroupIND = {HCIND,svIND,naIND};
GroupNames = {'HC','SV','NA'};
GraphType={'FAGraph','SCGraph'};
compPairs = [2 3;...
             1 1]
            
for g = 1
    for p = 1:2
        p
        G1=compPairs(1,p);
        G2=compPairs(2,p);
        testName=[GroupNames{G1} 'vs' GroupNames{G2}];
        
        %allmtx=allResults.(app.GraphType).('WeightMtx')(app.Group1Index);
        allmtx1=AllResults.(GraphType{g}).('WeightMtx')(GroupIND{G1});
        meanmtx1 = allmtx1{1};
        for k = 2:length(allmtx1)
            meanmtx1 = meanmtx1+allmtx1{k};
        end
        meanmtx1=meanmtx1./length(allmtx1);
        
        allmtx2=AllResults.(GraphType{g}).('WeightMtx')(GroupIND{G2});
        meanmtx2 = allmtx2{1};
        for k = 2:length(allmtx2)
            meanmtx2 = meanmtx2+allmtx2{k};
        end
        meanmtx2=meanmtx2./length(allmtx2);

%         for i = 1:size(meanmtx1,1)
%             for j = 1:size(meanmtx1,2)
%                 if(i~=j)
%                 diff = CoM(i,:)-CoM(j,:);
%                 dist = sqrt(diff(1)*diff(1)+diff(2)*diff(2)+diff(3)*diff(3));
%                 meanmtx1(i,j) = meanmtx1(i,j)/dist;
%                 meanmtx2(i,j) = meanmtx2(i,j)/dist;
%                 end
%             end
%         end

        
        diffVals=meanmtx1-meanmtx2;%difference between group 1 and 2
        thresh=mean(diffVals(diffVals~=0))-3*std(diffVals(diffVals~=0)); %find bottom 3 stddev of groupwise differences
        
        meanmtx1(diffVals>thresh | diffVals>0 )=0;%remove all connections above 3 stddev
        saveDir = fullfile(saveDirBase,GraphType{g},'NetworkCmpInOutThinningtoHubs',testName);
        mkdir(saveDir);
        
        
        %keepIndx = ones(1,numLab);
        %keepIndx=graphSubNetworks.(GraphType{g}).(GroupNames{G1}).p05;
        %rmIndx = ~keepIndx;
        %meanmtx1(rmIndx,:)=0;
        %meanmtx1(:,rmIndx)=0;
        
        rmIndx=Hubs.SCGraph.HCThresh.HC | graphSubNetworks.('Thickness').(GroupNames{G1}).p001;
        keepIndx=rmIndx;
        meanmtx1(rmIndx,~rmIndx)=0;
        meanmtx1(~rmIndx,rmIndx)=0;
        meanmtx1(~rmIndx2,~rmIndx)=0;
        
        plotNetwork(meanmtx1, Mesh, CoM,saveDir,testName,[],[0 max(meanmtx1(:))],1,keepIndx+.001)
    end
end