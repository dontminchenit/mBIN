inT = readtable('D:/Min/Dropbox (Personal)/Research/Projects/FTDAnalysis/Data/Final_Network_Path_Sheet_5_7_21_wCore.xlsx');
%inT = readtable('/Users/min/Dropbox (Personal)/Research/Projects/FTDAnalysis/Data/Final_Network_Path_Sheet_5_7_21_wCore.xlsx');


%tidxKeyCore = contains(inT.Properties.VariableNames,"KeyCore");
removeGMRegions = {'L_GM_ATC',...
    'L_GM_HIP',...
    'L_GM_M1',...
    'L_GM_mPFC',...
    'L_GM_PC',...
    'L_GM_pSTC',...
    'L_GM_S1',...
    'L_GM_V1',...
    'R_GM_ATC',...
    'R_GM_HIP',...
    'R_GM_M1',...
    'R_GM_mPFC',...
    'R_GM_PC',...
    'R_GM_PFC',...
    'R_GM_pSTC',...
    'R_GM_S1',...
    'R_GM_V1',...
    };

removeWMRegions = {'L_WM_ATC',...
    'L_WM_HIP',...
    'L_WM_M1',...
    'L_WM_mPFC',...
    'L_WM_PC',...
    'L_WM_pSTC',...
    'L_WM_S1',...
    'L_WM_V1',...
    'R_WM_ATC',...
    'R_WM_HIP',...
    'R_WM_M1',...
    'R_WM_mPFC',...
    'R_WM_PC',...
    'R_WM_PFC',...
    'R_WM_pSTC',...
    'R_WM_S1',...
    'R_WM_V1',...
    };


% clear some persistent Figures
H1=figure(10)    
clf    

H2=figure(20)
clf

%GM+WM Graph
baseSaveDir = './PaperExpBilateral_6_25_2021/GM+WM/';
mkdir(baseSaveDir);
%tidxCoreAll = contains(inT.Properties.VariableNames,"Core_");
%tidxCoreAll = tidxCoreAll & ~tidxKeyCore;
tidxCoreAll = contains(inT.Properties.VariableNames,"L_") | contains(inT.Properties.VariableNames,"R_");

removeRegions=[removeGMRegions, removeWMRegions];
[BothTauNetworkBurdenCorr,BothTDPNetworkBurdenCorr,MarkerNames,MarkerSizeTauAll,MarkerSizeTDPAll]=analyzeGraphs(inT,tidxCoreAll,removeRegions,baseSaveDir,1);

ExistingMarkers.MarkerNames=MarkerNames;
ExistingMarkers.MarkerSizeTau=MarkerSizeTauAll;
ExistingMarkers.MarkerSizeTDP=MarkerSizeTDPAll;

%GM Graph
baseSaveDir = './PaperExpBilateral_6_25_2021/GM/';
mkdir(baseSaveDir);
%tidxCoreGM = contains(inT.Properties.VariableNames,"Core_GM");
%tidxCoreGM = tidxCoreGM & ~tidxKeyCore;
tidxCoreGM = contains(inT.Properties.VariableNames,"L_GM") | contains(inT.Properties.VariableNames,"R_GM");

removeRegions=removeGMRegions;
[GMTauNetworkBurdenCorr,GMTDPNetworkBurdenCorr]=analyzeGraphs(inT,tidxCoreGM,removeRegions,baseSaveDir,0,ExistingMarkers);


%WM Graph
baseSaveDir = './PaperExpBilateral_6_25_2021/WM/';
mkdir(baseSaveDir);
%tidxCoreWM = contains(inT.Properties.VariableNames,"Core_WM");
%tidxCoreWM = tidxCoreWM & ~tidxKeyCore;
tidxCoreWM = contains(inT.Properties.VariableNames,"L_WM") | contains(inT.Properties.VariableNames,"R_WM");
removeRegions=removeWMRegions;
[WMTauNetworkBurdenCorr,WMTDPNetworkBurdenCorr]=analyzeGraphs(inT,tidxCoreWM,removeRegions,baseSaveDir,0,ExistingMarkers);



mName = {'Degree','ClusterCoeff','BetweenCen',};
mFullName = {'Degree', 'Clustering\nCoefficient', 'Betweenness\nCentrality'};

CmpCorr=nan(3,3);
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
    
        
        
        IndCorr(m,(j-1)*2 + 1)= CurrTauMetrics.([mName{m} '_corrPval']);
        IndCorr(m,(j-1)*2 + 2)= CurrTDPMetrics.([mName{m} '_corrPval']);
        
        IndRegP(m,(j-1)*2 + 1)= CurrTauMetrics.([mName{m} '_regpval']);
        IndRegP(m,(j-1)*2 + 2)= CurrTDPMetrics.([mName{m} '_regpval']);
    end
end

%Correlate networkMeasures

H=figure(10)
clf
%tiledlayout(3,6, 'Padding', 'none', 'TileSpacing', 'compact'); 

for m = 1:3
for ColID = 1:6
    
        switch ColID
            
            case 1
                CurrMetrics = GMTauNetworkBurdenCorr;
                xlabelName='GM ln(%AO)';
            case 2
                CurrMetrics = GMTDPNetworkBurdenCorr;
                xlabelName='GM ln(%AO)';
                
            case 3
                CurrMetrics = WMTauNetworkBurdenCorr;
                xlabelName='WM ln(%AO)';
            case 4
                CurrMetrics = WMTDPNetworkBurdenCorr;
                xlabelName='WM ln(%AO)';
                
            case 5
                CurrMetrics = BothTauNetworkBurdenCorr;
                xlabelName='GM & WM ln(%AO)';
            case 6
                CurrMetrics = BothTDPNetworkBurdenCorr;
                xlabelName='GM & WM ln(%AO)';
        end
    

    subplot(3,6,ColID + 6*(m-1));
    %nexttile
    X=CurrMetrics.Burden;
    Y=CurrMetrics.(mName{m});
    valid=~isnan(X') & ~isnan(Y);
    
    scatter(X,Y,15)
    xlabel(xlabelName);
    if (ColID==1)
        ylabel(sprintf(mFullName{m}))
    end
    
    if(m==3)
        switch ColID
            
            case 1
                xlabel('GM ln(%AO)')
            case 2
                xlabel('GM ln(%AO)')
                
            case 3
                xlabel('WM ln(%AO)')
            case 4
                xlabel('WM ln(%AO)')
                
            case 5
                xlabel('GM & WM ln(%AO)')
            case 6
                xlabel('GM & WM ln(%AO)')
        end
    end
    
    lsline
    corr=CurrMetrics.([mName{m} '_corr']);
    
    if(IndCorr(m,ColID)<=0.05)
        str=sprintf('r= %1.2f*',corr);    
    else
        str=sprintf('r= %1.2f',corr);
    end
    
    
    
    
    
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    set(gca,'fontsize',8)

end
end
baseSaveDir = './PaperExpBilateral_6_25_2021/';
saveName = 'NetworkMetricCmp'
print(H,fullfile(baseSaveDir,[ saveName '.png']),'-dpng','-r400'); 

% H=figure(20)
% saveName = 'TDPMetricCmp'
% print(H,fullfile(baseSaveDir,[ saveName '.png']),'-dpng','-r400'); 



%CompareMetrics