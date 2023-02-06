%convert Dan's data to same table format
T=readtable('D:\Min\Dropbox (Personal)\Research\Projects\FTDAnalysis\Data\DanData_Layer_project_master_1_22_2021.xlsx');

INDDID = unique(T.INDDID);
NS=length(IDList);

RegionList = unique(T.Region);
NR=length(RegionList);

Tnew = table(INDDID);

for r = 1:NR
    Tnew.(['L_GM_' RegionList{r}])=nan(NS,1);
    Tnew.(['R_GM_' RegionList{r}])=nan(NS,1);
end

Tnew.Tau1TDP2=nan(NS,1);
Tnew.CDx1_updated=cell(NS,1);


for t = 1:height(T)

newIdx = find(Tnew.INDDID == T.INDDID(t));
currRegion = [T.Hemi_updated{t} '_GM_' T.Region{t}];
Tnew.(currRegion)(newIdx) =  T.ln_layers_avg(t);
Tnew.Tau1TDP2(newIdx)=double(strcmpi(T.Cohort(t),'TDP'))+1;
Tnew.CDx1_updated{newIdx}=T.CDx1_updated{t};
end

%we have the new Table

inT=Tnew;

tidx = contains(inT.Properties.VariableNames,"GM");


cond=strcmpi(Tnew.CDx1_updated,'PSPS');

gmT_tau=inT(inT.Tau1TDP2==1 & cond,tidx);
gmT_tdp=inT(inT.Tau1TDP2==2 & cond,tidx);



for k = 1:2
    if(k == 1)
        gmT_Curr = gmT_tau;
    elseif(k==2)
        gmT_Curr = gmT_tdp;
    end

gmTvals=gmT_Curr{:,:};%log(gmT_Curr{:,:}*.01);    
%     
%table_data=gmTvals;
%table_data=gmT_Curr{:,:};
%figure(10)
%boxplot(table_data,'Labels',gmT_Curr.Properties.VariableNames);
%set(gca,'FontSize',10,'XTickLabelRotation',90)    
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

saveDir = './testResultsDan_PSPS';
mkdir(saveDir);
if(k==1)
saveName = 'tau_MTX.png';    
else
saveName = 'tdp_MTX.png';    
end

saveas(H,fullfile(saveDir,saveName));

% H=figure(2)
% imagesc(adjMtxCount'>5)
% set(gca,'xtick',[1:N],'xticklabel',gmT_Curr.Properties.VariableNames)
% set(gca,'ytick',[1:N],'yticklabel',gmT_Curr.Properties.VariableNames)
% set(gca,'TickLabelInterpreter','none')
% colorbar
% xtickangle(90)


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