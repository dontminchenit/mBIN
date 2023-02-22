%LoadDataSarah

%find center of mass, to show the location of nodes
meshlabels = labelmesh.attributes.attribute_array;
CoM=zeros(height(LabelLUT),3);
for i = 1:height(LabelLUT)
    currLabel = LabelLUT{i,1};
    idx = find(meshlabels==currLabel);
    CoM(i,:)=mean(labelmesh.points(idx,:),1);
end


H2=figure(2)
clf
panelAll = panel();
panelAll.pack(2,5);

panelAll.de.margin = 1;
panelAll.marginbottom = 0;
panelAll.margin = [1 1 1 1];

compPairs = [2 3 2 3;...
    1 1 3 2];

FAByStage = nan(5,2,100000);

for q = 1:2
%for s = 1:5
    for g = 2
        for s = 1:5            

            p
            G1=compPairs(1,q);
            G2=compPairs(2,q);
            testName=[GroupNames{G1} 'vs' GroupNames{G2}];
            
            allmtx1=AllResults.(graphNames{g}).('WeightMtx')(GroupInd{G1});
            meanmtx1 = allmtx1{1};
            for k = 2:length(allmtx1)
                meanmtx1 = meanmtx1+allmtx1{k};
            end
            meanmtx1=meanmtx1./length(allmtx1);
            caseEdgeVAls = meanmtx1;
            
            
            allmtx2=AllResults.(graphNames{g}).('WeightMtx')(GroupInd{G2});
            meanmtx2 = allmtx2{1};
            for k = 2:length(allmtx2)
                meanmtx2 = meanmtx2+allmtx2{k};
            end
            meanmtx2=meanmtx2./length(allmtx2);
            ctrlEdgeVAls = meanmtx1;
            
            diffVals=meanmtx1-meanmtx2;%difference between group 1 and 2
            thresh=mean(diffVals(diffVals~=0))-3*std(diffVals(diffVals~=0)); %find bottom 3 stddev of groupwise differences
            
            meanmtx1(diffVals>thresh | diffVals>0 )=0;%remove all connections above 3 stddev
            saveDir = fullfile(saveDirBase,graphNames{g},'NetworkCmp',testName);
            mkdir(saveDir);
            
            %keepIndx = ones(size(CoM(:,1)));
            
            
            
            if(q==1)
                 rmIndx=T_tauStagingLaus==s; 
                 pathName = 'Tau';
            else
                 rmIndx=T_tdpStagingLaus==s;
                 pathName = 'TDP';
            end
            keepIndx=rmIndx;
            %meanmtx1(rmIndx,~rmIndx)=0;
            %meanmtx1(~rmIndx,rmIndx)=0;
            %meanmtx1(~rmIndx,:)=0;
            meanmtx1(~rmIndx,~rmIndx)=0;
            
            stagecaseEdgeVAls=caseEdgeVAls;
            stagecaseEdgeVAls(~rmIndx,~rmIndx)=nan;
            
            stagectrlEdgeVAls=ctrlEdgeVAls;
            stagectrlEdgeVAls(~rmIndx,~rmIndx)=nan;
            
            val=stagecaseEdgeVAls(~isnan(stagecaseEdgeVAls));
            
            FAByStage(s,q,1:length(val)) = val;
            
            colorVec=keepIndx*s;% + (keepIndx & Hubs.FAGraph.HCThresh.HC');
            
            plotNetwork3(meanmtx1, labelmesh, CoM,saveDir,testName,[],[0 max(meanmtx1(:))],1, keepIndx+.001,panelAll,s,colorVec)
            
        end
    end
    
    print(H2,fullfile(saveDirBase,['fig6a_' pathName 'Stages.tif']),'-dpng','-r400');
end
%end