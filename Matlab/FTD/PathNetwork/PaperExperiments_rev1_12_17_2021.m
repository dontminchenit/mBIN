%inT = readtable('D:/Min/Dropbox (Personal)/Research/Projects/FTDAnalysis/Data/Final_Network_Path_Sheet_1_14_22_wCore_wHip.xlsx');
inT = readtable('D:/Min/Dropbox (Personal)/Research/Projects/FTD/FTDAnalysis/Data/Final_Network_Path_Sheet_12_17_21_wCore_wHip.xlsx');
%inT = readtable('/Users/min/Dropbox (Personal)/Research/Projects/FTDAnalysis/Data/Final_Network_Path_Sheet_1_14_22_wCore_wHip.xlsx');

tidxKeyCore = contains(inT.Properties.VariableNames,"KeyCore");
removeGMRegions = {'Core_GM_ATC',...
    'Core_GM_M1',...
    'Core_GM_mPFC',...
    'Core_GM_PC',...
    'Core_GM_pSTC',...
    'Core_GM_S1',...
    'Core_GM_V1',...
    };

removeWMRegions = {'Core_WM_ATC',...
    'Core_WM_M1',...
    'Core_WM_mPFC',...
    'Core_WM_PC',...
    'Core_WM_pSTC',...
    'Core_WM_S1',...
    'Core_WM_V1',...
    };


% clear some persistent Figures
H1=figure(10)    
clf    

H2=figure(20)
clf

baseSaveDir = './PaperExpCore_Rev1_2_11_2022';
mkdir(baseSaveDir)


%GM+WM Graph
SaveDir = fullfile(baseSaveDir,'GM+WM');
mkdir(SaveDir);
tidxCoreAll = contains(inT.Properties.VariableNames,"Core_");
tidxCoreAll = tidxCoreAll & ~tidxKeyCore;
removeRegions=[removeGMRegions, removeWMRegions];
[BothTauNetworkBurdenCorr,BothTDPNetworkBurdenCorr,MarkerNames,MarkerSizeTauAll,MarkerSizeTDPAll]=analyzeGraphs(inT,tidxCoreAll,removeRegions,SaveDir,1);

ExistingMarkers.MarkerNames=MarkerNames;
ExistingMarkers.MarkerSizeTau=MarkerSizeTauAll;
ExistingMarkers.MarkerSizeTDP=MarkerSizeTDPAll;

%GM Graph
SaveDir = fullfile(baseSaveDir,'GM/');
mkdir(SaveDir);
tidxCoreGM = contains(inT.Properties.VariableNames,"Core_GM");
tidxCoreGM = tidxCoreGM & ~tidxKeyCore;
removeRegions=removeGMRegions;
[GMTauNetworkBurdenCorr,GMTDPNetworkBurdenCorr]=analyzeGraphs(inT,tidxCoreGM,removeRegions,SaveDir,0,ExistingMarkers);


%WM Graph
SaveDir = fullfile(baseSaveDir,'WM/');
mkdir(SaveDir);
tidxCoreWM = contains(inT.Properties.VariableNames,"Core_WM");
tidxCoreWM = tidxCoreWM & ~tidxKeyCore;
removeRegions=removeWMRegions;
[WMTauNetworkBurdenCorr,WMTDPNetworkBurdenCorr]=analyzeGraphs(inT,tidxCoreWM,removeRegions,SaveDir,0,ExistingMarkers);



mName = {'Degree','ClusterCoeff','BetweenCen',};
mFullName = {'Degree', 'Clustering\nCoefficient', 'Betweenness\nCentrality'};

CmpCorr=nan(3,3);
IndCorrP=nan(3,6);
IndCorr=nan(3,6);
IndRegP=nan(3,6);
CmpCorr1UB=nan(3,3);
CmpCorr1LB=nan(3,3);
CmpCorr2UB=nan(3,3);
CmpCorr2LB=nan(3,3);

for j = 1:3
    
        switch j
            
            case 1
                CurrTauMetrics = GMTauNetworkBurdenCorr;
                CurrTDPMetrics = GMTDPNetworkBurdenCorr;
                
            case 2
                CurrTauMetrics = WMTauNetworkBurdenCorr;
                CurrTDPMetrics = WMTDPNetworkBurdenCorr;
                
            case 3
                CurrTauMetrics = BothTauNetworkBurdenCorr;
                CurrTDPMetrics = BothTDPNetworkBurdenCorr;
        end
        
    for m =1:3
        Corr1=CurrTauMetrics.([mName{m} '_corr']);
        Corr2=CurrTDPMetrics.([mName{m} '_corr']);
        Corr1LB=CurrTauMetrics.([mName{m} '_corrLB']);
        Corr1UB=CurrTauMetrics.([mName{m} '_corrUB']);
        Corr2LB=CurrTDPMetrics.([mName{m} '_corrLB']);
        Corr2UB=CurrTDPMetrics.([mName{m} '_corrUB']);
        
        NCorr1=CurrTauMetrics.([mName{m} '_corrN']);
        NCorr2=CurrTDPMetrics.([mName{m} '_corrN']);
        CmpCorr(m,j)=TestCorrSig(Corr1,Corr2,NCorr1,NCorr2);

        CmpCorr1UB(m,j)=Corr1+Corr1UB <= Corr2;
        CmpCorr1LB(m,j)=Corr1+Corr1LB >= Corr2;
        CmpCorr2UB(m,j)=Corr2+Corr1UB <= Corr1;
        CmpCorr2LB(m,j)=Corr2+Corr1LB >= Corr1;
    
         
        IndCorr(m,(j-1)*2 + 1) =CurrTauMetrics.([mName{m} '_corr']);
        IndCorr(m,(j-1)*2 + 2) =CurrTDPMetrics.([mName{m} '_corr']);
        IndCorrP(m,(j-1)*2 + 1)= CurrTauMetrics.([mName{m} '_corrPval']);
        IndCorrP(m,(j-1)*2 + 2)= CurrTDPMetrics.([mName{m} '_corrPval']);
        
        IndRegP(m,(j-1)*2 + 1)= CurrTauMetrics.([mName{m} '_regpval']);
        IndRegP(m,(j-1)*2 + 2)= CurrTDPMetrics.([mName{m} '_regpval']);
    end
end

%Correlate networkMeasures

H=figure(10)
clf
%tiledlayout(3,6, 'Padding', 'none', 'TileSpacing', 'compact'); 

for m = 1
for ColID = 1:6
    
    switch ColID
        case 1
            CurrMetrics = GMTauNetworkBurdenCorr;
            xlabelName='GM ln(%AO)';
            plotxLoc=1;
            plotyLoc=1;
        case 2
            CurrMetrics = GMTDPNetworkBurdenCorr;
            xlabelName='GM ln(%AO)';
            plotxLoc=1;
            plotyLoc=2;
        case 3
            CurrMetrics = WMTauNetworkBurdenCorr;
            xlabelName='WM ln(%AO)';
            plotxLoc=2;
            plotyLoc=1;
        case 4
            CurrMetrics = WMTDPNetworkBurdenCorr;
            xlabelName='WM ln(%AO)';
            plotxLoc=2;
            plotyLoc=2;
        case 5
            CurrMetrics = BothTauNetworkBurdenCorr;
            xlabelName='GM & WM ln(%AO)';
            plotxLoc=3;
            plotyLoc=1;
        case 6
            CurrMetrics = BothTDPNetworkBurdenCorr;
            xlabelName='GM & WM ln(%AO)';
            plotxLoc=3;
            plotyLoc=2;
    end
    
    subplot(3,2,ColID);
    %subplot(1,6,ColID + 6*(m-1));
%    nexttile
    X=CurrMetrics.Burden;
    Y=CurrMetrics.(mName{m});
    valid=~isnan(X') & ~isnan(Y);
    
    scatter(X,Y,15)
    
    %if (ColID==1 || ColID==4)
       % ylabel(sprintf(mFullName{m}))
    %end
    
%    if(ColID>4)
        switch ColID
            
            case 1
                ylabel('GM Node Strength')
            case 2
                %xlabel('GM ln(%AO)')
                
            case 3
                ylabel('WM Network Degree')
            case 4
                %xlabel('WM ln(%AO)')
                
            case 5
                ylabel('GM & WM Network Degree')
                xlabel('FTLD-Tau ln(%AO)')
            case 6
                xlabel('FTLD-TDP ln(%AO)')
        end
%    end
    
    lsline
    corr=CurrMetrics.([mName{m} '_corr']);
    
    if(IndCorrP(m,ColID)<=0.05)
        str=sprintf('r= %1.2f*',corr);    
    else
        str=sprintf('r= %1.2f',corr);
    end
    
    
    
    
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(gca,'fontsize',8)

end
end
SaveDir = baseSaveDir;
saveName = 'NetworkMetricCmp'
print(H,fullfile(SaveDir,[ saveName '.png']),'-dpng','-r400'); 

% H=figure(20)
% saveName = 'TDPMetricCmp'
% print(H,fullfile(baseSaveDir,[ saveName '.png']),'-dpng','-r400'); 



%CompareMetrics