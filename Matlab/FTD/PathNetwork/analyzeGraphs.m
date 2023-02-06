function [TauNetworkBurdenCorr,TDPNetworkBurdenCorr, MarkerNames,MarkerSizeTau,MarkerSizeTDP]=analyzeGraphs(inTable,TableIdxCore,removeRegions,baseSaveDir,GMWMOffset,ExistingMarkers)

plotGraphsFlag = 1;

%set Tau/TDP data
Table_Tau=inTable(inTable.Tau2TDP1==2 & inTable.Include ==1 ,TableIdxCore);
Table_TDP=inTable(inTable.Tau2TDP1==1 & inTable.Include ==1 ,TableIdxCore);


%do this for both pathology
for Tau1TDP2 = 1:2
    if(Tau1TDP2 == 1)
        TableCurrent = removevars(Table_Tau,removeRegions);
        TableCurrentLogVals=log(TableCurrent{:,:}+.00015);%Log Transform
    elseif(Tau1TDP2==2)
        TableCurrent = removevars(Table_TDP,removeRegions);
        TableCurrentLogVals=log(TableCurrent{:,:}+.00015);%Log Transform
    end
    
    N = size(TableCurrentLogVals,2);
    
    CorrAdjMtx = nan(N,N);
    CorradjMtxLowConf = nan(N,N);
    CorrAdjMtxCount = nan(N,N);
    CorrAdjMtxPval = nan(N,N);
    for i = 1:N
        for j = (i+1):N
            A= TableCurrentLogVals(:,i);
            B= TableCurrentLogVals(:,j);
            [R,P,RLO,RUP]= corrcoef(A, B, 'alpha', 0.05,'rows','complete');
            CorrAdjMtxPval(i,j) = P(1,2);
            CorrAdjMtx(i,j) = R(1,2);
            CorradjMtxLowConf(i,j) = RLO(1,2);
            CorrAdjMtxCount(i,j) = sum(~isnan(A) & ~isnan(B));
        end
    end
    
    
    SN=length(TableCurrent.Properties.VariableNames);%Number of Regions
    baseNames = cell(1,SN);%BaseNames of Regions (remove GM/WM from name)
    
    for j = 1:SN
        nameSplit = strsplit(TableCurrent.Properties.VariableNames{j},'_');
        baseNames{j} = nameSplit{:,3};
    end
    
    %Plot CorrAdjMtx
    if(plotGraphsFlag)
        H=figure(2);
        clf
        imagesc(CorrAdjMtx')
        set(gca,'xtick',[1:N],'xticklabel',baseNames)
        set(gca,'ytick',[1:N],'yticklabel',baseNames)
        %    set(gca,'xtick',[1:N],'xticklabel',gmT_Curr.Properties.VariableNames)
        %    set(gca,'ytick',[1:N],'yticklabel',gmT_Curr.Properties.VariableNames)
        set(gca,'TickLabelInterpreter','none')
        colormap('jet')
        %colorbar
        caxis([0 1])
        xtickangle(90)
        set(gca,'FontSize',12)
        
        if(Tau1TDP2==1)
            saveName = 'tau_MTX.png';
        else
            saveName = 'tdp_MTX.png';
        end
        print(H,fullfile(baseSaveDir,saveName),'-dpng','-r400');%save
    end
    
    MarkerNames=TableCurrent.Properties.VariableNames;
    %plot number of samples available for each correlation
    if(plotGraphsFlag)
        
        H=figure(2)
        clf
        imagesc(CorrAdjMtxCount')
        set(gca,'xtick',[1:N],'xticklabel',TableCurrent.Properties.VariableNames)
        set(gca,'ytick',[1:N],'yticklabel',TableCurrent.Properties.VariableNames)
        set(gca,'TickLabelInterpreter','none')
        colorbar
        xtickangle(90)
        if(Tau1TDP2==1)
            saveName = 'tau_Count.png';
            savePrefix ='tau_'
        else
            saveName = 'tdp_Count.png';
            savePrefix ='tdp_'
        end

        saveas(H,fullfile(baseSaveDir,saveName));
    end

    NetworkBurdenCorr=relateConnectivityAndBurden(CorrAdjMtx,nanmean(TableCurrentLogVals,1),MarkerNames);
    
    if(Tau1TDP2 == 1)
        TauMtx = CorrAdjMtx;
        TauCount = CorrAdjMtxCount;
        TauNetworkBurdenCorr = NetworkBurdenCorr;
        if(~exist('ExistingMarkers','var'))
            MarkerSizeTau = normalize(nanmean(TableCurrentLogVals,1),'range')+1;
            MarkerSizeTau=MarkerSizeTau.^2;
            MarkerSizeTau(isnan(MarkerSizeTau)) = .1;
        else
            [matched matchIdx] = ismember(MarkerNames,ExistingMarkers.MarkerNames);
            MarkerSizeTau=ExistingMarkers.MarkerSizeTau(matchIdx);
        end
    elseif(Tau1TDP2==2)
        TDPMtx = CorrAdjMtx;
        TDPCount = CorrAdjMtxCount;
        TDPNetworkBurdenCorr = NetworkBurdenCorr;
        if(~exist('ExistingMarkers','var'))
            MarkerSizeTDP = normalize(nanmean(TableCurrentLogVals,1),'range')+1;
            MarkerSizeTDP=MarkerSizeTDP.^2;
            MarkerSizeTDP(isnan(MarkerSizeTDP)) = .1;
        else
            [matched matchIdx] = ismember(MarkerNames,ExistingMarkers.MarkerNames);
            MarkerSizeTDP=ExistingMarkers.MarkerSizeTDP(matchIdx);
        end
        
    end
    
    if(Tau1TDP2==1)
        MarkerSize=MarkerSizeTau;
        saveName = 'tau_graph';
    else
        MarkerSize=MarkerSizeTDP;
        saveName = 'TDP_graph';
    end
    
    plotMtx = CorrAdjMtx;
    plotMtx(isnan(plotMtx)) = 0;
    if(plotGraphsFlag)
        plotGraphs(plotMtx,TableCurrent,baseSaveDir,saveName,GMWMOffset,MarkerSize,Tau1TDP2)
    end
end

%find places where TAU/TDP is significantly lower

%valid = TauCount>3 & TDPCount>3;
valid = TauCount>3 & TDPCount>3;
zTau=atanh(TauMtx(valid));
zTDP=atanh(TDPMtx(valid));
NTau=TauCount(valid);
NTDP=TDPCount(valid);
zObs = (zTau - zTDP)./ ((1 ./ (NTau - 3)) + (1 ./ (NTDP - 3))).^(.5);
Probs = normcdf(zObs);

TauSigLower=zeros(N,N);
TauSigLower(valid) = Probs <= .05;

TDPSigLower=zeros(N,N);
TDPSigLower(valid) = Probs >= .95;

if(plotGraphsFlag)
    
    for TauTDPcmp = 1:4
        if TauTDPcmp==1
            plotMtx = TauSigLower;
            MarkerSize=MarkerSizeTau;
            saveName = 'tauLowerThanTDP_comp';
            Tau1TDP2=2;
        elseif TauTDPcmp==2
            plotMtx = TDPSigLower;
            MarkerSize=MarkerSizeTDP;
            saveName = 'TDPLowerThanTau_comp';
            Tau1TDP2=2;
        elseif TauTDPcmp==3
            plotMtx = TDPSigLower;
            MarkerSize=MarkerSizeTau;
            saveName = 'tauGreaterThanTDP_comp';
            Tau1TDP2=1;
        elseif TauTDPcmp==4
            plotMtx = TauSigLower;
            MarkerSize=MarkerSizeTDP;
            saveName = 'TDPGreaterThanTau_comp';
            Tau1TDP2=2;
        end
        
        
        
        plotGraphs(plotMtx,TableCurrent,baseSaveDir,saveName,GMWMOffset,MarkerSize,Tau1TDP2)
    end
end
