function plotBrain3Points(CoM,GIImesh,data,alphaThresh,panelAll,pSelect_row,cmapTypeForce,pointSize,surfDisp,setViews,ViewLoc,fillFlag)
%function plotBrain3Points(CoM,GIImesh,data,CRange,saveDir,saveName,alphaThresh,panelAll,pSelect,cmapTypeForce,pointSize,surfDisp)


N = length(data);
LeftIdx = GIImesh.R0_L1_Index;
RightIdx = ~GIImesh.R0_L1_Index;
surfLabels = GIImesh.cdata;
surfval=GIImesh.cdata;

XMid=0;
ZMid=15;

if(exist('surfDisp','var') && ~isempty(surfDisp))
for i = 1:N
    currLabel = GIImesh.LabelLUT{i,1};
    currVal = surfDisp(i);
    surfval(surfLabels==currLabel) = currVal;
end
end
if(alphaThresh>0)
    PlotTypeRange = 2;
    flipCmap=1;
    surfval(surfLabels==0) = 1;
    data(~isnan(data))=double(data(~isnan(data))<alphaThresh);
else
    PlotTypeRange = 1;
    flipCmap=0;
end

CoMDisp = CoM;

if(~exist('setViews','var') || isempty(setViews))
    setViews=1:6;
    ViewLoc=1:6;
end


if(~exist('fillFlag','var') || isempty(fillFlag))
    fillFlag=3;
end




for i = PlotTypeRange
    for vSel = 1:length(setViews) %6 views L,Lmedial,R,Rmedial,Dorsal,Ventral
%        valDisp = val;

        v=setViews(vSel);
        if(exist('surfDisp','var') && ~isempty(surfDisp))
            valDisp = surfval;
        else
            valDisp = zeros(length(surfval),1);        
        end
        regionSel = (data~=0);
        
        switch v 
            case 1
                viewA=-90;
                viewE=15;
                valDisp(RightIdx) = nan;
                regionSel(CoM(:,1)>(XMid-15))=0;
                lightStr=2;
            case 2
                viewA=90;
                viewE=35;%
                valDisp(RightIdx) = nan;
                regionSel((CoM(:,1)>XMid)|(CoM(:,1)<=-15))=0;
                lightStr=2;
            case 3
                viewA=90;
                viewE=15;
                valDisp(LeftIdx) = nan;
                regionSel(CoM(:,1)<(XMid+15))=0;
                lightStr=2;
            case 4
                viewA=-90;
                viewE=35;
                valDisp(LeftIdx) = nan;
                regionSel((CoM(:,1)<XMid)|(CoM(:,1)>=15))=0;
                lightStr=2;
            case 5
                viewA=0;
                viewE=90;
                regionSel(CoM(:,3)<(ZMid))=0;
                lightStr=1;
            case 6
                viewA=-180;
                viewE=-90;
                regionSel(CoM(:,3)>=(ZMid))=0;
                lightStr=1;
                
        end
    if (i==1) %if normal map
        plotType=[];
        cmapType = 1; %jet
        CRange=[0 1];
    else %if threshold map
        valDisp(~isnan(valDisp))=double(valDisp(~isnan(valDisp))>alphaThresh);
        plotType='thresh';
        CRange=[0 1];
        cmapType = 2; %gray to red
        flipCmap=0;
    end
    
    if (exist('cmapTypeForce','var'))
        cmapType=cmapTypeForce;
    end

    if (~exist('pointSize','var'))
       pointSize = 50;
    end
    
%     H=figure(1)
%     clf
%     plotMesh(mesh, valDisp, viewA,viewE,CRange,flipCmap)
%     colorbar
%     setTightFig
%     
%     if ( ~exist(saveDir, 'dir'))
%         mkdir(saveDir)
%     end
%     saveas(H,fullfile(saveDir,[ saveName '_' plotNames{v} '_' plotType '.png']))
    figure(2)
    if(pSelect_row>0)
    panelAll(pSelect_row,ViewLoc(vSel)).select();
%    plotMesh(mesh, valDisp, viewA,viewE,CRange,flipCmap)
    
 
    switch cmapType
        case 1
            cmap=brewermap(100,'Greens');
            cmap(1,:)=[1 1 1];
            CRange=[0 1];
        case 2
            cmap = flipud(brewermap(100,'RdBu'));
            cmap(50,:)=[1 1 1];
            CRange=[-.25 .25];
        case 3
            cmap = flipud(brewermap(100,'PiYG'));
            cmap(50,:)=[1 1 1];
            CRange=[-.25 .25];
        case 4
            CRange=[0 6];
            cmap = zeros(6,3);
            cmap(1,:) = [0 0 0];
            cmap(2,:) = [1 0 0];
            cmap(3,:) = [0.3281 0.0859 0.7031];
            cmap(4,:) = [0 .5 0];
            cmap(5,:) = [1 .5 0];
            cmap(6,:) = [1, 1, 0];
            cmap(7,:) = [240, 163, 255]/255;
        case 5
            cmap=brewermap(100,'Greens');
            %cmap=brewermap(100,'Blues');
            %cmap = flipud(brewermap(100,'RdBu'));
            cmap(size(cmap,1),:)=[1 1 1];
            cmap = flipud(cmap);
            flipCmap=1;
        case 6
            CRange=[-1 1];
            cmap = zeros(3,3);
            cmap(1,:) = [45/255 115/255 178/255];
            cmap(2,:) = [106/255 26/255 74/255];
            cmap(3,:) = [1 1 1];
            cmap(4,:) = [0/255 82/255 33/255];%[106/255 26/255 74/255];%[90/255 39/255 41/255];
            cmap(5,:) = [228/255 63/255 41/255];
        case 7 
            CRange=[0 6];
            cmap=colormap('hot');
    end
    
    if (fillFlag == 1)
    hold on
    scatter3(CoMDisp(regionSel,1),CoMDisp(regionSel,2),CoMDisp(regionSel,3),pointSize,data(regionSel)','filled');
    hold off
    elseif (fillFlag == 0)
    scatter3(CoMDisp(regionSel,1),CoMDisp(regionSel,2),CoMDisp(regionSel,3),pointSize,data(regionSel)');
    colormap(gca,cmap);
    hold on
    t=plotgiiSurf(GIImesh.giiSurface_Both, valDisp, viewA,viewE,CRange,flipCmap,cmapType,lightStr,cmap);
    %alpha(t,.05)
    alpha(t,.2)
    hold off
    elseif (fillFlag == 3)
        scatter3(CoMDisp(regionSel,1),CoMDisp(regionSel,2),CoMDisp(regionSel,3),pointSize,data(regionSel)','filled');
        colormap(gca,cmap);
        hold on
        t=plotgiiSurf(GIImesh.giiSurface_Both, valDisp, viewA,viewE,CRange,flipCmap,cmapType,lightStr,cmap);
        %alpha(t,.2)
        alpha(t,.2)
        hold off
    end
    
    set(gca,'visible','off')
    end
    end
    %H2=figure(2)
    %saveas(H2,fullfile(saveDir,[ saveName '_' plotNames{v} '_' plotType 'thumbnail.png']))
end
