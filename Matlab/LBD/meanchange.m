networkNames={'DorsAttn','VentAttn','Limbic','SomMot','Default','Cont','_Vis_','All'};

meanAll=nan(2,8);
stdAll=nan(2,8);

for ii= 1:8
    if(ii==8)
        currNetwork=true(1,length(NetworkDataGeneral.Schaefer400x7.LUT.Label_Name2));
    else
        currNetwork=contains(NetworkDataGeneral.Schaefer400x7.LUT.Label_Name2,networkNames{ii})';
    end


networkNames(ii)
yesdata=thickyesChange(:,currNetwork);
yesdata([6 8 20],:)=[];
nodata=thicknoChange(:,currNetwork);
nodata([15 19 20],:)=[];
meanAll(1,ii)=mean(yesdata(:),'omitnan');
meanAll(2,ii)=mean(nodata(:),'omitnan');

stdAll(1,ii)=std(yesdata(:),'omitnan');
stdAll(2,ii)=std(nodata(:),'omitnan');

end
