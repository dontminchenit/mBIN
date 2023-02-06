inT = readtable('D:\Min\Dropbox (Personal)\Research\Projects\FTDAnalysis\Data\Working Network Path Sheet 1_29_21.xlsx');

baseSaveDir = './PaperExpCore/RatioGMWM/';
mkdir(baseSaveDir);

tidxL = contains(inT.Properties.VariableNames,"L_");
tidxGM = contains(inT.Properties.VariableNames,"L_GM");
tidxWM = contains(inT.Properties.VariableNames,"L_WM");

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
removeRegionsGM = {'L_GM_aCING',...
'L_GM_pCING',...
'L_GM_aINS',...
'L_GM_M1', ...
'L_GM_mPFC', ...
'L_GM_PC', ...
'L_GM_pSTC', ...
'L_GM_S1', ...
'L_GM_V1', ...
'L_GM_ATC', ...
};
removeRegionsWM = {'L_WM_aCING',...
'L_WM_pCING',...
'L_WM_aINS',...
'L_WM_M1', ...
'L_WM_mPFC', ...
'L_WM_PC', ...
'L_WM_pSTC', ...
'L_WM_S1', ...
'L_WM_V1', ...
'L_WM_ATC', ...
};


%gmT_tau=inT(inT.Tau2TDP1==2 & strcmpi(inT.ClinicalCAT,'bvFTD') ,tidx);
%gmT_tdp=inT(inT.Tau2TDP1==1 & strcmpi(inT.ClinicalCAT,'bvFTD') ,tidx);

gmT_tau=inT(inT.Tau2TDP1==2 & inT.Include ==1 ,tidxGM);
gmT_tdp=inT(inT.Tau2TDP1==1 & inT.Include ==1 ,tidxGM);

wmT_tau=inT(inT.Tau2TDP1==2 & inT.Include ==1 ,tidxWM);
wmT_tdp=inT(inT.Tau2TDP1==1 & inT.Include ==1 ,tidxWM);



for k = 1:2
    if(k == 1)
        gmT_Curr = removevars(gmT_tau,removeRegionsGM);
        wmT_Curr = removevars(wmT_tau,removeRegionsWM);
 %       gmT_Curr = gmT_tau
 %       wmT_Curr = wmT_tau
        gmwmT_AllTau =gmT_Curr;
        writetable(gmT_Curr,fullfile(baseSaveDir,'TauGMComb.xlsx'));
        gmTvals=log(gmT_Curr{:,:}+.00015);
        wmTvals=log(wmT_Curr{:,:}+.00015);
    elseif(k==2)
        gmT_Curr = removevars(gmT_tdp,removeRegionsGM);
        wmT_Curr = removevars(wmT_tdp,removeRegionsWM);
%        gmT_Curr = gmT_tdp
%        wmT_Curr = wmT_tdp
        gmwmT_AllTDP =gmT_Curr;
        writetable(gmT_Curr,fullfile(baseSaveDir,'TDPGMComb.xlsx'));
        gmTvals=log(gmT_Curr{:,:}+.00015);    
        wmTvals=log(wmT_Curr{:,:}+.00015);    
        
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
        
        A= gmTvals(:,i)./wmTvals(:,i);
        B= gmTvals(:,j)./wmTvals(:,j)
        [R,P,RLO,RUP]= corrcoef(A, B, 'alpha', 0.05,'rows','complete');
        
        
        adjMtxCount(i,j) = sum(~isnan(A) & ~isnan(B));
        
        if(adjMtxCount(i,j) <= 5)
            adjMtx(i,j) = nan;
            adjMtxLB(i,j) = nan;
        else
            adjMtx(i,j) = R(1,2);
            adjMtxLB(i,j) = RLO(1,2);
        end
        
        
        
%         Aw= wmTvals(:,i)
%         Ag= gmTvals(:,i);
%         Bw= wmTvals(:,j)
%         Bg= gmTvals(:,j);
%         [Rw,P,RLOw,RUP]= corrcoef(Aw, Bw, 'alpha', 0.05,'rows','complete');
%         [Rg,P,RLOg,RUP]= corrcoef(Ag, Bg, 'alpha', 0.05,'rows','complete');
%         
%         adjMtx(i,j) = Rw(1,2)/Rg(1,2);
%         adjMtxLB(i,j) = RLOw(1,2)/RLOg(1,2);
%         adjMtxCount(i,j) = sum(~isnan(Aw) & ~isnan(Bw) & ~isnan(Ag) & ~isnan(Bg));
    end
end

H=figure(k)
imagesc(adjMtx')
set(gca,'xtick',[1:N],'xticklabel',gmT_Curr.Properties.VariableNames)
set(gca,'ytick',[1:N],'yticklabel',gmT_Curr.Properties.VariableNames)
set(gca,'TickLabelInterpreter','none')
colormap('jet')
colorbar
caxis([-1 1])
xtickangle(90)


if(k==1)
saveName = 'tau_MTX.png';    
else
saveName = 'tdp_MTX.png';    
end

saveas(H,fullfile(baseSaveDir,saveName));

 H=figure(2)
 imagesc(adjMtxCount' > 5)
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
        TauAllMtx = adjMtx;
        TauMtx = adjMtx;
        TauLBmtx = adjMtxLB;
    elseif(k==2)
        TDPAllMtx = adjMtx;
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