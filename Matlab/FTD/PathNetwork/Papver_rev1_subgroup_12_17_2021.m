inT = readtable('D:/Min/Dropbox (Personal)/Research/Projects/FTDAnalysis/Data/Final_Network_Path_Sheet_12_17_21_wCore_wHip.xlsx');
%inT = readtable('/Users/min/Dropbox (Personal)/Research/Projects/FTDAnalysis/Data/Final_Network_Path_Sheet_12_17_21_wCore_wHip.xlsx');

tidxKeyCore = contains(inT.Properties.VariableNames,"KeyCore");
removeGMRegions=[];
removeWMRegions=[];
removeGMRegions = {'Core_GM_ATC',...
    'Core_GM_HIP',...
    'Core_GM_M1',...
    'Core_GM_mPFC',...
    'Core_GM_PC',...
    'Core_GM_pSTC',...
    'Core_GM_S1',...
    'Core_GM_V1',...
    };

removeWMRegions = {'Core_WM_ATC',...
    'Core_WM_HIP',...
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

baseSaveDir = './PaperExpCore_Rev1_Subgroup_1_21_2022/';
mkdir(baseSaveDir)

%TauSubTypes ={'CBD','PID','PSP','PSP','AllTau'};
%TDPSubTypes ={'A','B','C','E','C'};

TauSubTypes ={'CBD','PID','PSP','All','All','All','All'};
TDPSubTypes ={'All', 'All', 'All', 'A','B','C','All'};


for q = 1:7
%GM+WM Graph
SaveDir = fullfile(baseSaveDir,'GM+WM');
mkdir(SaveDir);
tidxCoreAll = contains(inT.Properties.VariableNames,"Core_");
tidxCoreAll = tidxCoreAll & ~tidxKeyCore;
%tidxCoreAll = tidxKeyCore;
removeRegions=[removeGMRegions, removeWMRegions];
[BothTauNetworkBurdenCorr,BothTDPNetworkBurdenCorr,MarkerNames,MarkerSizeTauAll,MarkerSizeTDPAll]=analyzeGraphsSubGroup(inT,tidxCoreAll,removeRegions,SaveDir,1,[],TauSubTypes{q},TDPSubTypes{q});
%[BothTauNetworkBurdenCorr,BothTDPNetworkBurdenCorr,MarkerNames,MarkerSizeTauAll,MarkerSizeTDPAll]=analyzeGraphs(inT,tidxCoreAll,removeRegions,baseSaveDir,1);
end