
data = path_yesAD;

%X=data(:,5); %stmc
%Y=data(:,3); %m1
%M=data(:,1); %ACing

X=data(:,6); %stri
Y=data(:,1); %aCING
M=data(:,5); %STMC


[paths, stats] = mediation(X, Y, M, 'plots', 'verbose','boottop');