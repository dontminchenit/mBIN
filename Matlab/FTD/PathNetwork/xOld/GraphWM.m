inT = readtable('./Working Network Path spreadhseet 1_20_2021.xlsx');



tidx = contains(inT.Properties.VariableNames,"WM");
inT.('Tau1_TDP2') = zeros(height(inT),1);

removeRegions = {'L_WM_ATC',...
    'L_WM_aINS',...
    'L_WM_mPFC', ...
    'L_WM_PC', ...
    'L_WM_pCING', ...
    'L_WM_pSTC', ...
    'L_WM_S1', ...    %'L_WM_PFC',...
    'L_WM_V1',...
    'R_WM_V1',...
    'L_WM_M1',...
    'R_WM_M1',...
    'R_WM_ATC',...
    'R_WM_aINS',...
    'R_WM_mPFC', ...
    'R_WM_PC', ...
    'R_WM_pCING', ...
    'R_WM_pSTC', ...
    'R_WM_PFC',...
    'R_WM_S1'
    
};

tidxL = contains(inT.Properties.VariableNames,"L_WM");

for j = 1:length(tidxL)
    
    if(tidxL(j))
        currName = inT.Properties.VariableNames{j};
        RName = currName;
        RName(1) = 'R';
        if(contains(inT.Properties.VariableNames,RName))
            inT.(currName) = nanmax([inT.(currName) inT.(RName)],[],2);
        end
    end
    
end
removeRegions = {'L_WM_ATC',...
    'L_WM_aINS',...
    'L_WM_mPFC', ...
    'L_WM_PC', ...
    'L_WM_pCING', ...
    'L_WM_pSTC', ...
    'L_WM_S1', ...    %'L_WM_PFC',...
    'L_WM_V1',...
    'L_WM_M1' ...
    };

tidx=tidxL


%gmT_tau=inT(inT.Tau2TDP1==2 & strcmpi(inT.ClinicalCAT,'bvFTD') ,tidx);
%gmT_tdp=inT(inT.Tau2TDP1==1 & strcmpi(inT.ClinicalCAT,'bvFTD') ,tidx);

gmT_tau=inT(inT.Tau2TDP1==2 ,tidx);
gmT_tdp=inT(inT.Tau2TDP1==1 ,tidx);



for k = 1:2
    if(k == 1)
        gmT_Curr = removevars(gmT_tau,removeRegions);
        %gmT_Curr = gmT_tau;
    elseif(k==2)
        gmT_Curr = removevars(gmT_tdp,removeRegions);
        %gmT_Curr = gmT_tdp;
    end

gmTvals=log(gmT_Curr{:,:} + .00015);    
%     
%table_data=gmTvals;
% table_data=gmT_Curr{:,:};
% figure(10)
% boxplot(table_data,'Labels',gmT_Curr.Properties.VariableNames);
% set(gca,'FontSize',10,'XTickLabelRotation',90)    
%     
% %     



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

saveDir = './WMResults'
mkdir(saveDir);
if(k==1)
saveName = 'tau_MTX.png';    
else
saveName = 'tdp_MTX.png';    
end

saveas(H,fullfile(saveDir,saveName));

H=figure(2)
imagesc(adjMtxCount'>5)
set(gca,'xtick',[1:N],'xticklabel',gmT_Curr.Properties.VariableNames)
set(gca,'ytick',[1:N],'yticklabel',gmT_Curr.Properties.VariableNames)
set(gca,'TickLabelInterpreter','none')
colorbar
xtickangle(90)


    if(k == 1)
        TauWMMtx = adjMtx;
        TauLBmtx = adjMtxLB;
    elseif(k==2)
        TDPWMMtx = adjMtx;
        TDPLBmtx = adjMtxLB;
    end

if(k==1)
saveName = 'tau_graph';    
else
saveName = 'tdp_graph';    
end
    
%     plotMtx = adjMtx;
% plotGraphs
% 
% 

end

% CompareGraphs