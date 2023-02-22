

MN = 5;
SN = 6;
allpvalsLeft = zeros(MN,SN);
allpvalsRight = zeros(MN,SN);


for m = 1:5
    
    currMetric = allMetricData{m};

    for s = 1:6

        
    vals = squeeze(currMetric(s,1,:));
    vals(isnan(vals)) = [];
    vals2 = squeeze(currMetric(s,2,:));
    vals2(isnan(vals2)) = [];
    
    [h,p] = ttest2(vals,vals2,'Tail','left');
    
    allpvalsLeft(m,s) = p;

    [h,p] = ttest2(vals,vals2,'Tail','right');

    allpvalsRight(m,s) = p;

    
    end
    
end