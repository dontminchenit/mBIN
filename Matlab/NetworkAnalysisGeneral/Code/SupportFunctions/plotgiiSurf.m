function t=plotgiiSurf(giiSurf, valDisp, viewA,viewE,CRange,flipCmap,cmapType,lightStr,cmapIn)

t=trisurf(giiSurf.faces,giiSurf.vertices(:,1),giiSurf.vertices(:,2),giiSurf.vertices(:,3),valDisp(:),'EdgeColor','none');


%set(gca, 'color', 'k');
%lighting gouraud;
lighting flat;
%    a1 = camlight('right');
%    a2 = camlight('left');
%    a3 = camlight('headlight');

for i = 1:lightStr
    a4 = camlight(0,90,'infinite');
    a4 = camlight(0,-90,'infinite');
    %    a4 = camlight(90,0,'infinite');
    %    a4 = camlight(-90,0,'infinite');
    
end


%camlight;
%camlight(-80,-10);
axis equal

view(viewA,viewE)
%a4 = camlight('headlight','infinite');
%a4 = camlight('headlight','infinite');
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
colormap(c);

%     xlim([min(giiSurf.vertices(:,1)) max(giiSurf.vertices(:,1))])
%     ylim([min(giiSurf.vertices(:,2)) max(giiSurf.vertices(:,2))])
%     zlim([min(giiSurf.vertices(:,3)) max(giiSurf.vertices(:,3))])
