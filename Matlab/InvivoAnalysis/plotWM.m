H2=figure(2)
clf
panelAll = panel();
panelAll.pack(2,6);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];

cRange = [0 .2];
 MakerVec=ones(1,400);
 colorVec=ones(1,400);
LabelNames=[];
 pSelect=1;
 
 
plotNetwork3(schaefer400x17_gfa_end_connectivity, schaefer400x7SurfaceMeshFromGII, CoM,LabelNames,cRange,MakerVec,panelAll,pSelect,colorVec,0)
