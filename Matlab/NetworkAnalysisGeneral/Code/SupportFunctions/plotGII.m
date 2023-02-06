function plotGII(mesh, valDisp, viewA,viewE,CRange,flipCmap)

  trisurf(mesh.triangles,mesh.points(:,1),mesh.points(:,2),mesh.points(:,3),valDisp(:),'EdgeColor','none');

    
    set(gca, 'color', 'k');
    lighting gouraud;
    a1 = camlight('right');
    a2 = camlight('left');
    a3 = camlight('headlight');
    a4 = camlight(0,-90);
    axis equal
    
    view(viewA,viewE)

    c = jet;
    if(flipCmap)
        c = flipud(c);
    end
    caxis(CRange)
    colormap(c);