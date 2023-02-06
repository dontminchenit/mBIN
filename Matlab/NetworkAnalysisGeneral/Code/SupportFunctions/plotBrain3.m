function plotBrain3(GIImesh,data,CRange,saveDir,saveName,alpha,panelAll,pSelect,cmapTypeForce)


N = length(data);
LeftIdx = GIImesh.R0_L1_Index;
RightIdx = ~GIImesh.R0_L1_Index;
surfLabels = GIImesh.cdata;
val=GIImesh.cdata;

for i = 1:N
    currLabel = GIImesh.LabelLUT{i,1};
    currVal = data(i);
    val(surfLabels==currLabel) = currVal;
end
if(alpha>0)
    PlotTypeRange = 2;
    flipCmap=1;
    val(surfLabels==0) = 1;
else
    PlotTypeRange = 1;
    flipCmap=0;
end


for i = PlotTypeRange
    for v = 1:6 %5 views L,Lmedial,R,Rmedial,Dorsal,Ventral
        valDisp = val;

        switch v 
            case 1
                viewA=-90;
                viewE=15;
                valDisp(RightIdx) = nan;
                lightStr=2;
            case 2
                viewA=90;
                viewE=35;%
                valDisp(RightIdx) = nan;
                lightStr=2;
            case 3
                viewA=90;
                viewE=15;
                valDisp(LeftIdx) = nan;
                lightStr=2;
            case 4
                viewA=-90;
                viewE=35;
                valDisp(LeftIdx) = nan;
                lightStr=2;
            case 5
                viewA=0;
                viewE=90;
                lightStr=1;
            case 6
                viewA=-180;
                viewE=-90;
                lightStr=1;
                
        end
    if (i==1) %if normal map
        plotType=[];
        cmapType = 1; %jet
    else %if threshold map
        valDisp(~isnan(valDisp))=double(valDisp(~isnan(valDisp))<alpha);
        plotType='thresh';
        CRange=[0 1];
        cmapType = 2; %gray to red
        flipCmap=0;
    end
    
    if (exist('cmapTypeForce','var'))
        cmapType=cmapTypeForce;
    end

%     H=figure(1)
%     clf
%     plotMesh(mesh, valDisp, viewA,viewE,CRange,flipCmap)
%     colorbar
%     setTightFig
%     
     if ( ~exist(saveDir, 'dir'))
         mkdir(saveDir)
     end
%     saveas(H,fullfile(saveDir,[ saveName '_' plotNames{v} '_' plotType '.png']))
    figure(2)
    if(pSelect>0)
    panelAll(pSelect,v).select();
%    plotMesh(mesh, valDisp, viewA,viewE,CRange,flipCmap)
    plotgiiSurf(GIImesh.giiSurface_Both, valDisp, viewA,viewE,CRange,flipCmap,cmapType,lightStr)
    
    set(gca,'visible','off')
    end
    end
    %H2=figure(2)
    %saveas(H2,fullfile(saveDir,[ saveName '_' plotNames{v} '_' plotType 'thumbnail.png']))
end
