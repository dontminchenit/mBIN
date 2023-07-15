function t=plotgiiSurf(giiSurf, valDisp, viewA,viewE,CRange,flipCmap,cmapType,lightStr,cmapIn)

t=trisurf(giiSurf.faces,giiSurf.vertices(:,1),giiSurf.vertices(:,2),giiSurf.vertices(:,3),valDisp(:),'EdgeColor','none');

lighting flat;

for i = 1:lightStr
    a4 = camlight(0,90,'infinite');
    a4 = camlight(0,-90,'infinite'); 
end

axis equal

view(viewA,viewE)

if(~exist('cmapIn','var'))
    if(cmapType==1)
        c = flipud(brewermap(100,'RdBu'));
    elseif(cmapType==2)
        c = brewermap(100,'Greens');
    elseif(cmapType==3)
        c = flipud(brewermap(100,'PuOr'));
    elseif(cmapType==4)
        c = brewermap(100,'Oranges');
    elseif(cmapType==5)
        c= colormap('jet');
        c(1,:)=[1 1 1];
    end
else
    c=cmapIn;
end

if(flipCmap)
    c = flipud(c);
    c(1,:)=[1 1 1];
end
clim(CRange)
