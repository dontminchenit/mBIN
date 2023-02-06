inT = readtable('D:\Min\Dropbox (Personal)\Research\Projects\FTDAnalysis\Data\Working Network Path Sheet 1_29_21.xlsx');

baseSaveDir = './PaperExpCore/WM/';
mkdir(baseSaveDir);

tidxL = contains(inT.Properties.VariableNames,"L_WM");

for j = 1:length(tidxL)
    
    if(tidxL(j))
        LName = inT.Properties.VariableNames{j};
%        NewName = currName;
%        NewName(1) = '_';
%        NewName=['Core' NewName];
        
%        inT.(NewName) = nan(height(inT),1);
        
        for i = 1:height(inT)
            
            CoreName = LName;
            CoreName(1) = inT.CoreHemisphere{i};
            
            if(~isnan(inT.(CoreName)(i)))
                inT.(LName)(i)=inT.(CoreName)(i);
            end
        end
        
    end
    
end
% removeRegions = {'L_GM_aCING',...
% 'L_GM_M1', ...
% 'L_GM_mPFC', ...
% 'L_GM_PC', ...
% 'L_GM_pSTC', ...
% 'L_GM_S1', ...
% 'L_GM_V1', ...
% 'L_GM_ATC', ...
% 'L_GM_pCING', ...
% 'L_GM_HIP', ...
% 'L_GM_IFC', ...
% 'L_GM_aITC' ...
% 'L_GM_SPC' ...
% 'L_WM_SPC' ...
% 'L_WM_ATC',...
% 'L_WM_aINS',...
% 'L_WM_mPFC', ...
% 'L_WM_PC', ...
% 'L_WM_pCING', ...
% 'L_WM_pSTC', ...
% 'L_WM_S1', ...    %'L_WM_PFC',...
% 'L_WM_V1',...
% 'L_WM_M1'
% };


%gmT_tau=inT(inT.Tau2TDP1==2 & strcmpi(inT.ClinicalCAT,'bvFTD') ,tidx);
%gmT_tdp=inT(inT.Tau2TDP1==1 & strcmpi(inT.ClinicalCAT,'bvFTD') ,tidx);

gmT_tau=inT(inT.Tau2TDP1==2 & inT.Include ==1 ,tidxL);
gmT_tdp=inT(inT.Tau2TDP1==1 & inT.Include ==1 ,tidxL);



for k = 1:2
    if(k == 1)
%        gmT_Curr = removevars(gmT_tau,removeRegions);
        gmT_Curr = gmT_tau
        gmwmT_AllTau =gmT_Curr;
        writetable(gmT_Curr,fullfile(baseSaveDir,'TauGMComb.xlsx'));
        gmTvals=log(gmT_Curr{:,:});    
    elseif(k==2)
%        gmT_Curr = removevars(gmT_tdp,removeRegions);
        gmT_Curr = gmT_tdp
        gmwmT_AllTDP =gmT_Curr;
        writetable(gmT_Curr,fullfile(baseSaveDir,'TDPGMComb.xlsx'));
        gmTvals=log(gmT_Curr{:,:}+.00015);    
    end


% table_data=gmT_Curr{:,:};
% figure(10)
% boxplot(table_data,'Labels',gmT_Curr.Properties.VariableNames);
% set(gca,'FontSize',10,'XTickLabelRotation',90)    
    
N = size(gmTvals,2);

adjMtx = zeros(N,N);
adjMtxLB = zeros(N,N);
adjMtxCount = zeros(N,N);
for i = 1:N
    for j = (i+1):N
        A= gmTvals(:,i);
        B= gmTvals(:,j);
        [R,P,RLO,RUP]= corrcoef(A, B, 'alpha', 0.05,'rows','complete');
        
        adjMtx(i,j) = R(1,2);
        adjMtxLB(i,j) = RLO(1,2);
        adjMtxCount(i,j) = sum(~isnan(A) & ~isnan(B));
    end
end

H=figure(k)
imagesc(adjMtx')
set(gca,'xtick',[1:N],'xticklabel',gmT_Curr.Properties.VariableNames)
set(gca,'ytick',[1:N],'yticklabel',gmT_Curr.Properties.VariableNames)
set(gca,'TickLabelInterpreter','none')
colormap('jet')
colorbar
caxis([0 1])
xtickangle(90)


if(k==1)
saveName = 'tau_MTX.png';    
else
saveName = 'tdp_MTX.png';    
end

saveas(H,fullfile(baseSaveDir,saveName));

 H=figure(2)
 imagesc(adjMtxCount')
 set(gca,'xtick',[1:N],'xticklabel',gmT_Curr.Properties.VariableNames)
 set(gca,'ytick',[1:N],'yticklabel',gmT_Curr.Properties.VariableNames)
 set(gca,'TickLabelInterpreter','none')
 colorbar
 xtickangle(90)

if(k==1)
saveName = 'tau_Count.png';    
else
saveName = 'tdp_Count.png';    
end

saveas(H,fullfile(baseSaveDir,saveName));

 
 

    if(k == 1)
        TauWMMtx = adjMtx;
        TauMtx = adjMtx;
        TauLBmtx = adjMtxLB;
    elseif(k==2)
        TDPWMMtx = adjMtx;
        TDPMtx = adjMtx;
        TDPLBmtx = adjMtxLB;
    end

if(k==1)
saveName = 'tau_graph';    
else
saveName = 'tdp_graph';    
end
    
plotMtx = adjMtx;
plotGraphs



end

CompareGraphs
%CompareMetrics