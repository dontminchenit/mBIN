

IDS = TPhases.INDDID;

wscoreT= readtable(fullfile(currBase,'PhaseData','wscore_merged_clinphen_tdptype.csv'));

adjPhaseFilters=cell(5,2);

%adjmtx
for p = 1:2
    if(p==1)
        currPathStage=T_tauStagingLaus;
        pathname = 'tau'
    elseif(p==2)
        currPathStage=T_tdpStagingLaus;
        pathname = 'tdp';
        
    end
    
    for s=1:5
        
        PhaseA=s
        PhaseB=s+1
        
        nodesA=currPathStage==PhaseA;
        nodesB=currPathStage==PhaseB;
        
        adjmtx2=adjmtx>0;
        
        %remove nodes not in either phase
        adjmtx2(~(nodesA | nodesB),:)=0;
        adjmtx2(:,~(nodesA | nodesB))=0;
        
        %remove nodes internal to each phase
        adjmtx2(nodesA,nodesA)=0;
        adjmtx2(nodesB,nodesB)=0;
        
        adjPhaseFilters{s,p}=adjmtx2;
        
        %
        %         H2=figure(1)
        %         clf
        %         imshow(adjmtx2)
        %         colormap(jet)
        %         caxis([0 max(adjmtx2(:))])
        %         colorbar
        %         title([pathname ' - Phase ' num2str(PhaseA) ' to ' num2str(PhaseB)])
        %         print(H2,fullfile(saveDirBase,[pathname '_prion_stage' num2str(s)   '.tif']),'-dpng','-r400');
        %
        
    end
    
end

countData=zeros(length(IDS),5);

for i = 1:length(IDS)
    
    currIDs = IDS(i);
    
    Tau1Tdp0=strcmpi(TPhases.Group(i),'tau');
    
    for s=1:5
        
        if(Tau1Tdp0)
            currPathStage=T_tauStagingLaus;
            adjfilter= adjPhaseFilters{s,1};
        else
            currPathStage=T_tdpStagingLaus;
            adjfilter= adjPhaseFilters{s,2};
        end
        
        
        indx= find(wscoreT.INDDID==currIDs);
        
        currWscore=wscoreT{indx,20:467}<-0.5;
        
        adjmtx3=zeros(size(adjmtx));
        adjmtx3(currWscore,currWscore)=1;
        adjmtx3=adjmtx3.*adjfilter;
        adjmtx3=triu(adjmtx3);
        %countData(i,s)=sum(adjmtx3(:))/(sum(currPathStage==s));
        countData(i,s)=sum(adjmtx3(:))/(sum(currPathStage==s)*sum(currPathStage==(s+1)));
        %countData(i,s)=sum(adjmtx3(:))/(sum(adjfilter(:))/2);
    end
    %     figure(2)
    %     clf
    %     imshow(adjmtx3)
    %     colormap(jet)
    %     caxis([0 max(adjmtx3(:))])
    %     colorbar
    %
    %     title(['Subject - atrophy plot phase 1-2'])
    %print(H2,fullfile(saveDirBase,[pathname '_prion_stage' num2str(s)   '.tif']),'-dpng','-r400');
end

tauIndx=strcmpi(TPhases.Group,'tau');
tauData=countData(tauIndx,:);
TDPData=countData(~tauIndx,:);

allData=nan(length(IDS),5,2);

allData(1:length(tauData),:,1)=tauData;
allData(1:length(TDPData),:,2)=TDPData;
labelNames={'Tau','TDP'}
names={'1->2','2->3','3->4','4->5','5->6'};
x = names;
y = allData;

H2=figure(2);
clf

iosr.statistics.boxPlot(x,y, ...
    'style','hierarchy',...
    'xSeparator',true,...
    'symbolColor','k',...
    'symbolMarker','.',...
    'medianColor','k',...
    'boxcolor','auto',...
    'showScatter',true,...
    'scatterColor', 'r', ...
    'scatterMarker', '.', ...
    'groupLabels', labelNames, ...
    'groupLabelFontSize', 10 ...
    );

pvalues = nan(5,1);
print(H2,fullfile(saveDirBase,'fig8b_prioncount.tif'),'-dpng','-r400');

for s=1:5
    [H pval ci stats]= ttest2(allData(:,s,1),allData(:,s,2))
    pvalues(s)=p;
end