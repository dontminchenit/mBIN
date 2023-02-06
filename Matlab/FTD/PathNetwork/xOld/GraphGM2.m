inT = readtable('./Working Network Path spreadhseet 1_20_2021.xlsx');



tidx = contains(inT.Properties.VariableNames,"GM");
inT.('Tau1_TDP2') = zeros(height(inT),1);

removeRegions = {'L_GM_aCING',...
    'L_GM_M1', ...
'L_GM_mPFC', ...
'L_GM_PC', ...
'L_GM_pSTC', ...
'L_GM_S1', ...
'L_GM_V1', ...
'R_GM_M1', ...
'R_GM_mPFC', ...
'R_GM_PC', ...
'R_GM_PFC', ...
'R_GM_pSTC', ...
'R_GM_S1', ...
'R_GM_V1', ...
'R_GM_ATC', ...
'L_GM_ATC', ...
'R_GM_pCING', ...
'L_GM_pCING', ...
'R_GM_HIP', ...
'L_GM_HIP', ...
'R_GM_ANG', ...
'L_GM_IFC', ...
'L_GM_aITC' ...
'R_GM_IFC', ...
'R_GM_aITC' ...
};

tidxL = contains(inT.Properties.VariableNames,"L_GM");

for j = 1:length(tidxL)
    
    if(tidxL(j))
        currName = inT.Properties.VariableNames{j};
        RName = currName;
        RName(1) = 'R';
        if(contains(inT.Properties.VariableNames,RName))
            inT.(currName) = nanmax([inT.(currName) inT.(RName)],2);
        end
    end
    
end
removeRegions = {'L_GM_aCING',...
'L_GM_M1', ...
'L_GM_mPFC', ...
'L_GM_PC', ...
'L_GM_pSTC', ...
'L_GM_S1', ...
'L_GM_V1', ...
'L_GM_ATC', ...
'L_GM_pCING', ...
'L_GM_HIP', ...
'L_GM_IFC', ...
'L_GM_aITC' ...
'L_GM_SPC' ...
};


%gmT_tau=inT(inT.Tau2TDP1==2 & strcmpi(inT.ClinicalCAT,'bvFTD') ,tidx);
%gmT_tdp=inT(inT.Tau2TDP1==1 & strcmpi(inT.ClinicalCAT,'bvFTD') ,tidx);

gmT_tau=inT(inT.Tau2TDP1==2 ,tidxL);
gmT_tdp=inT(inT.Tau2TDP1==1 ,tidxL);



for k = 1:2
    if(k == 1)
        gmT_Curr = removevars(gmT_tau,removeRegions);
        writetable(gmT_Curr,'./TauGMComb.xlsx');
    elseif(k==2)
        gmT_Curr = removevars(gmT_tdp,removeRegions);
        writetable(gmT_Curr,'./TDPGMComb.xlsx');
    end

gmTvals=log(gmT_Curr{:,:}+.00015);    
%     
%table_data=gmTvals;
table_data=gmT_Curr{:,:};
figure(10)
boxplot(table_data,'Labels',gmT_Curr.Properties.VariableNames);
set(gca,'FontSize',10,'XTickLabelRotation',90)    
    
%     



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

saveDir = './GMResults_2_5'
mkdir(saveDir);
if(k==1)
saveName = 'tau_MTX.png';    
else
saveName = 'tdp_MTX.png';    
end

saveas(H,fullfile(saveDir,saveName));

 H=figure(2)
 imagesc(adjMtxCount'>10)
 set(gca,'xtick',[1:N],'xticklabel',gmT_Curr.Properties.VariableNames)
 set(gca,'ytick',[1:N],'yticklabel',gmT_Curr.Properties.VariableNames)
 set(gca,'TickLabelInterpreter','none')
 colorbar
 xtickangle(90)


    if(k == 1)
        TauMtx = adjMtx;
        TauLBmtx = adjMtxLB;
    elseif(k==2)
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