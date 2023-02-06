
TauCrossCov=nan(N,1);
TDPCrossCov=nan(N,1);

for i = 1:N

    
            NTau=sum(~isnan(Tau_W(:,i)) & ~isnan(pathTau(:,i)));
        NTDP=sum(~isnan(TDP_W(:,i)) & ~isnan(pathTDP(:,i)));

        if(NTau > 3)
            TauCrossCov(i)=corr(Tau_W(:,i),pathTau(:,i),'Rows','complete');
        end
        if(NTDP > 3)
            TDPCrossCov(i)=corr(TDP_W(:,i),pathTDP(:,i),'Rows','complete');
        end
    
    
end

CrossData= [TauCrossCov TDPCrossCov];
H=figure(1)
clf
imagesc(CrossData)
colorbar
set(gca,'ytick',[1:N],'yticklabel',regNames)
set(gca,'TickLabelInterpreter','none','FontSize',5)
print(H,fullfile(baseSaveDir,[ saveName '_crossCov.png']),'-dpng','-r400');

%create pannels
H2=figure(2)
clf
panelAll = panel();
panelAll.pack(2,6);

panelAll.de.margin = 1;
panelAll.marginbottom = 1;
panelAll.marginleft = -10;
panelAll.margin = [0 0 0 0];


plotBrain3Points(currCoM,schaefer400x7SurfaceMeshFromGII,TauCrossCov,0,panelAll,1,2,10)
plotBrain3Points(currCoM,schaefer400x7SurfaceMeshFromGII,TDPCrossCov,0,panelAll,2,2,10)
print(H2,fullfile(baseSaveDir,[ saveName '_crossCovBrain.png']),'-dpng','-r400');