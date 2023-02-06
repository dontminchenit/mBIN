function plotNetwork3_double(adjMtx, GIImesh, CoM,LabelNames,cRange,MakerVec,MakerVec2,panelAll,pSelect,colorVec,colorVec2,blue0orng1brwn2,surfDisp)
%,UnderlayData)

surfLabels = GIImesh.cdata;
surfval=GIImesh.cdata;
if(exist('surfDisp','var') && ~isempty(surfDisp))
    for i = 1:N
        currLabel = GIImesh.LabelLUT{i,1};
        currVal = surfDisp(i);
        surfval(surfLabels==currLabel) = currVal;
    end
    surfvalDisp = surfval;
else
    surfvalDisp = zeros(length(surfval),1);
end

%make graph
adjMtx(isnan(adjMtx))=0;
adjMtx=triu(adjMtx);
G = graph(adjMtx,'upper');
plotNames={'L-Lateral','Dorsal'};

%H=figure(figNum);
%clf
for v = 1:2
    
    
    panelAll(v,pSelect).select();
    
    switch v
        case 1
            viewA=-90;
            viewE=15;
        case 2
            viewA=0;
            viewE=90;
    end
    
    
    ColorVRGB = zeros(length(colorVec),3);
    

     cmap=colormap('hot');

     c1=cmap(round(256*1/6),:);
     c2=cmap(round(256*2/6),:);
     c3=cmap(round(256*3/6),:);
     c4=cmap(round(256*4/6),:);
     c5=cmap(round(256*5/6),:);
     c6=cmap(round(256*6/6),:);

     ColorVRGB(colorVec==0,:) = repmat(c6, sum(colorVec==0), 1);%c1
     ColorVRGB(colorVec==1,:) = repmat(c1, sum(colorVec==1), 1);%c1
     ColorVRGB(colorVec==2,:) = repmat(c2, sum(colorVec==2), 1);%c1
     ColorVRGB(colorVec==3,:) = repmat(c3, sum(colorVec==3), 1);%c1
     ColorVRGB(colorVec==4,:) = repmat(c4, sum(colorVec==4), 1);%c1
     ColorVRGB(colorVec==5,:) = repmat(c5, sum(colorVec==5), 1);%c1
     ColorVRGB(colorVec==6,:) = repmat(c6, sum(colorVec==6), 1);%c1
     ColorVRGB(colorVec==7,:) = repmat(c6, sum(colorVec==7), 1);%c1


    
    ColorVRGB(colorVec==0,:) = repmat([1 1 1], sum(colorVec==0), 1);%c1
    ColorVRGB(colorVec==1,:) = repmat([1 0 0], sum(colorVec==1), 1);%c1
    ColorVRGB(colorVec==2,:) = repmat([0.3281 0.0859 0.7031], sum(colorVec==2), 1);%c1
    ColorVRGB(colorVec==3,:) = repmat([0 .5 0], sum(colorVec==3), 1);%c1
    ColorVRGB(colorVec==4,:) = repmat([1 .5 0], sum(colorVec==4), 1);%c1
    ColorVRGB(colorVec==5,:) = repmat([1 1 0], sum(colorVec==5), 1);%c1
    ColorVRGB(colorVec==6,:) = repmat([240, 163, 255]/255, sum(colorVec==6), 1);%c1
    ColorVRGB(colorVec==7,:) = repmat([240,163,255]/255, sum(colorVec==7), 1);%c1
    
    ColorVRGB2 = zeros(length(colorVec),3);

    ColorVRGB2(colorVec2==0,:) = repmat([1 1 1], sum(colorVec2==0), 1);%c1
    ColorVRGB2(colorVec2==1,:) = repmat([1 0 0], sum(colorVec2==1), 1);%c1
    ColorVRGB2(colorVec2==2,:) = repmat([0.3281 0.0859 0.7031], sum(colorVec2==2), 1);%c1
    ColorVRGB2(colorVec2==3,:) = repmat([0 .5 0], sum(colorVec2==3), 1);%c1
    ColorVRGB2(colorVec2==4,:) = repmat([1 .5 0], sum(colorVec2==4), 1);%c1
    ColorVRGB2(colorVec2==5,:) = repmat([1 1 0], sum(colorVec2==5), 1);%c1
    ColorVRGB2(colorVec2==6,:) = repmat([240, 163, 255]/255, sum(colorVec2==6), 1);%c1
    ColorVRGB2(colorVec2==7,:) = repmat([240,163,255]/255, sum(colorVec2==7), 1);%c1
    
    
    if(blue0orng1brwn2==1)
        %Orange lines
        h=plot(G,'XData',CoM(:,1),'YData',CoM(:,2),'ZData',CoM(:,3),'interpreter','none','EdgeColor',[251, 97, 26]/255,'linewidth',3,'MarkerSize',4.25*MakerVec,'Marker','o','NodeColor',ColorVRGB);
    elseif(blue0orng1brwn2==2)
        %Orange lines
        h=plot(G,'XData',CoM(:,1),'YData',CoM(:,2),'ZData',CoM(:,3),'interpreter','none','EdgeColor',[61, 43, 31]/255,'linewidth',3,'MarkerSize',4.25*MakerVec,'Marker','o','NodeColor',ColorVRGB);
    else
        %blue lines
        h1=scatter3(CoM(:,1),CoM(:,2),CoM(:,3),25*MakerVec.*MakerVec,ColorVRGB,'o','filled','MarkerFaceAlpha',.5);
        hold on
        h2=scatter3(CoM(:,1),CoM(:,2),CoM(:,3),15*MakerVec2.*MakerVec2,ColorVRGB2,'o','LineWidth', 2,'MarkerFaceAlpha',.5);
%        h2=scatter3('XData',CoM(:,1),'YData',CoM(:,2),'ZData',CoM(:,3),'interpreter','none','linewidth',3,'MarkerSize',4.25*MakerVec2,'Marker','o','NodeColor',ColorVRGB2,'MarkerFaceAlpha',.1);
        hold off
    end
    
    showLabels=0;
    for i = 1:length(LabelNames)
        if(showLabels)
            labelnode(h,i,LabelNames{i});
        else
            %labelnode(h1,i,'');
            %labelnode(h2,i,'');
        end
    end
    axis equal
    weight=nonzeros(adjMtx);
    %h.EdgeCData=weight;
            caxis([0  .15])
    %cRange
    %caxis(cRange);
    %colormap(jet)
    %colormap(gca,brewermap(100,'PiYG'))
    hold on
    
    t=plotgiiSurf(GIImesh.giiSurface_Both, surfvalDisp, viewA,viewE,[0 .15],0,4,1);
    
    alpha(t,.05)
    hold off
    view(viewA,viewE)
    setTightFig
    set(gca,'visible','off')
    set(gcf, 'PaperPosition', [-2.5313    2.2083   13.5625    6.5833]);
    
end
% colorbar
%        saveas(H,fullfile(saveDir,[ saveName '.png']))
end
