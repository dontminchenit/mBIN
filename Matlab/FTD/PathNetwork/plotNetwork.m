function plotNetwork(adjMtx, Mesh, CoM,saveDir,saveName,LabelNames,cRange,figNum,MakerVec,CoMoffsetV1, CoMoffsetV2,GMWMOffset,Tau1TDP2)

%make graph
adjMtx=triu(adjMtx);
G = graph(adjMtx,'upper');
plotNames={'L-Lateral','Dorsal'};

H=figure(figNum);
clf
CoMAdj = CoM;
meshOffset = 0;
for v = 1%:2
    %subplot(1,2,v)
    switch v
        case 1
            viewA=-90;
            viewE=15;
            if(GMWMOffset)
                CoMAdj = CoM+CoMoffsetV1;
                meshOffset = min(CoMoffsetV1,[],1);
            end
        case 2
            viewA=0;
            viewE=90;
            if(GMWMOffset)
                CoMAdj = CoM+CoMoffsetV2;
                meshOffset = min(CoMoffsetV2,[],1);
            end
    end
    
    if(Tau1TDP2==1)
        NodeColor = [1 .5 0];
    else
        NodeColor = [205 56, 234]./255;
    end
    
    if(GMWMOffset)
       MarkerSize = 5*MakerVec;
       FontSize = 10;
       LWidths=1;
    else
       MarkerSize = 15*MakerVec; 
       FontSize = 30;
       LWidths=2;
    end
    
    
    h=plot(G,'XData',CoMAdj(:,1),'YData',CoMAdj(:,2),'ZData',CoMAdj(:,3),'interpreter','none','linewidth',1,'MarkerSize',MarkerSize,'Marker','o','NodeColor',NodeColor,'LineWidth',LWidths);
    highlight(h,find(strcmpi(LabelNames,'HIP')),'NodeColor',[0 .4 0])
    for i = 1:length(LabelNames)
        %labelnode(h,i,LabelNames{i})
    end
    h.NodeLabel = {};
    h.NodeFontSize = FontSize;
    axis equal
    weight=nonzeros(adjMtx);
    h.EdgeCData=weight;
    %        caxis([0  5000])
    caxis(cRange);
    colormap(jet)
    hold on
    t=trisurf(Mesh.triangles,Mesh.points(:,1),Mesh.points(:,2),Mesh.points(:,3),'FaceColor', [.9 .9 .9],'EdgeColor','none');
    alpha(t,.1)
    if(GMWMOffset)
        t2=trisurf(Mesh.triangles,Mesh.points(:,1)+meshOffset(1),Mesh.points(:,2)+meshOffset(2),Mesh.points(:,3)+meshOffset(3),'FaceColor', [.9 .9 .9],'EdgeColor','none');
        alpha(t2,.1)
    end
    
    hold off
    view(viewA,viewE)
    setTightFig
    set(gca,'visible','off')
    set(gcf, 'PaperPosition', [-2.5313    2.2083   13.5625    6.5833]);
    
end
%colorbar
%        saveas(H,fullfile(saveDir,[ saveName '.png']))
print(H,fullfile(saveDir,[ saveName '.png']),'-dpng','-r400');
end
