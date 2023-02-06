inLabelL = 'D:\Min\Dropbox (Personal)\Research\Projects\PPA\Data\lausanne250\lh.lausanne250.label.gii';
inSurfL = 'D:\Min\Dropbox (Personal)\Research\Projects\PPA\Data\lausanne250\lh.white.surf.gii';
inLabelR = 'D:\Min\Dropbox (Personal)\Research\Projects\PPA\Data\lausanne250\rh.lausanne250.label.gii';
inSurfR = 'D:\Min\Dropbox (Personal)\Research\Projects\PPA\Data\lausanne250\rh.white.surf.gii';


%first we go through each label can convert 
giiLabel_L=gifti(inLabelL);
giiSurface_L=gifti(inSurfL);
oldCdata_L = giiLabel_L.cdata;
newCData_L = nan(size(giiLabel_L.cdata));

%Left
for i = 1:length(oldCdata_L)
    oldVal = oldCdata_L(i);
    if(oldVal==0 || oldVal==4) %Unknown and corpuscallosum
        newCData_L(i) = 0;    
    else
    oldIndex = find(giiLabel_L.labels.key==oldVal);
    oldName = giiLabel_L.labels.name{oldIndex};
    newIndex = find(strcmpi(LabelLUT.Label_Name,['left-' oldName]));
    newVal = LabelLUT.Label_ID(newIndex);
    newCData_L(i) = newVal;
    end
end

giiLabel_R=gifti(inLabelR);
giiSurface_R=gifti(inSurfR);
oldCdata_R = giiLabel_R.cdata;
newCData_R = nan(size(giiLabel_R.cdata));

%Right
for i = 1:length(oldCdata_R)
    oldVal = oldCdata_R(i);
    if(oldVal==0 || oldVal==4) %Unknown and corpuscallosum
        newCData_R(i) = 0;    
    else
    oldIndex = find(giiLabel_R.labels.key==oldVal);
    oldName = giiLabel_R.labels.name{oldIndex};
    newIndex = find(strcmpi(LabelLUT.Label_Name,['right-' oldName]));
    newVal = LabelLUT.Label_ID(newIndex);
    newCData_R(i) = newVal;
    end
end

giiSurface_Both = struct;
offset=size(giiSurface_R.vertices,1);
giiSurface_Both.faces = [giiSurface_R.faces; (giiSurface_L.faces+offset)];
giiSurface_Both.vertices = [giiSurface_R.vertices; (giiSurface_L.vertices)];

R0_L1_Index = ones(size(giiSurface_Both.vertices,1),1);
R0_L1_Index(1:offset) = 0;
newCData_Both = [newCData_R; newCData_L];

Lausanne250SurfaceMeshFromGII=struct;
Lausanne250SurfaceMeshFromGII.giiSurface_Both=giiSurface_Both;
Lausanne250SurfaceMeshFromGII.R0_L1_Index = logical(R0_L1_Index);
Lausanne250SurfaceMeshFromGII.cdata = newCData_Both;
Lausanne250SurfaceMeshFromGII.LabelLUT = LabelLUT;

save('./Lausanne250SurfaceMeshFromGII.mat','Lausanne250SurfaceMeshFromGII');

%next we combine into one surface

% figure(3)
% data = zeros(numLab,1);
% data(314)=1;
% plotBrain2(labelmesh,data,LabelLUT,[0 1],0,[],0,[],1)
