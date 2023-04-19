function plotNetwork3(adjMtx, GIImesh, CoM,LabelNames,cRange,MakerVec,panelAll,pSelectCol,colorVec,blue0orng1brwn2jet3,surfDisp)

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

adjMtx(isnan(adjMtx))=0;
adjMtx=triu(adjMtx);
G = graph(adjMtx,'upper');
plotNames={'L-Lateral','Dorsal'};

for v = 1:2
    panelAll(v,pSelectCol).select();
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

     c1
     c2
     c3
     c4
     c5
     c6
        
     ColorVRGB

     ColorVRGB(colorVec==0,:) = repmat(c6, sum(colorVec==0), 1);%c1
     ColorVRGB(colorVec==1,:) = repmat(c1, sum(colorVec==1), 1);%c1
     ColorVRGB(colorVec==2,:) = repmat(c2, sum(colorVec==2), 1);%c1
     ColorVRGB(colorVec==3,:) = repmat(c3, sum(colorVec==3), 1);%c1
     ColorVRGB(colorVec==4,:) = repmat(c4, sum(colorVec==4), 1);%c1
     ColorVRGB(colorVec==5,:) = repmat(c5, sum(colorVec==5), 1);%c1
     ColorVRGB(colorVec==6,:) = repmat(c6, sum(colorVec==6), 1);%c1
     ColorVRGB(colorVec==7,:) = repmat(c6, sum(colorVec==7), 1);%c1    

     ColorVRGB

    if(blue0orng1brwn2jet3==1)
        h=plot(G,'XData',CoM(:,1),'YData',CoM(:,2),'ZData',CoM(:,3),'interpreter','none','EdgeColor',[251, 97, 26]/255,'linewidth',3,'MarkerSize',4.25*MakerVec,'Marker','o','NodeColor',ColorVRGB);
    elseif(blue0orng1brwn2jet3==2)
        h=plot(G,'XData',CoM(:,1),'YData',CoM(:,2),'ZData',CoM(:,3),'interpreter','none','EdgeColor',[61, 43, 31]/255,'linewidth',3,'MarkerSize',4.25*MakerVec,'Marker','o','NodeColor',ColorVRGB);
    elseif(blue0orng1brwn2jet3==0)
        h=plot(G,'XData',CoM(:,1),'YData',CoM(:,2),'ZData',CoM(:,3),'interpreter','none','linewidth',3,'MarkerSize',4.25*MakerVec,'Marker','o','NodeColor',ColorVRGB);
    elseif(blue0orng1brwn2jet3==4)


        Cdata=colorVec';
        Cdata(Cdata>2)=2;
        Cdata(Cdata<-2)=-2;
        cmin = -2;
        cmax = 2;
        cmap=colormap('jet');
        m = length(cmap);
        index = fix((Cdata-cmin)/(cmax-cmin)*m)+1; 

        RGB = squeeze(ind2rgb(index,cmap));

        h=plot(G,'XData',CoM(:,1),'YData',CoM(:,2),'ZData',CoM(:,3),'interpreter','none','linewidth',3,'MarkerSize',4.25*MakerVec,'Marker','o','NodeColor',RGB);
    else
        h=plot(G,'XData',CoM(:,1),'YData',CoM(:,2),'ZData',CoM(:,3),'interpreter','none','linewidth',3,'MarkerSize',4.25*MakerVec,'Marker','o','NodeColor',ColorVRGB);
        h.EdgeCData=G.Edges.Weight;

        clim(cRange);
        colormap(jet);
    end
    
    showLabels=1;
    for i = 1:length(LabelNames)
        if(showLabels)
            labelnode(h,i,LabelNames{i});
        else
            labelnode(h,i,'');
        end
    end
    axis equal
    
    hold on
    
    % Creates the 3D Brain surface
t=plotgiiSurf(GIImesh.giiSurface_Both, surfvalDisp, viewA,viewE,cRange,0,5,1);
    
    alpha(t,.05)
    hold off
    view(viewA,viewE)
    setTightFig
    set(gca,'visible','off')
    set(gcf, 'PaperPosition', [-2.5313    2.2083   13.5625    6.5833]);
    
end
end
