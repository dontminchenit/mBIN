function [rescaledDataGM, rescaledDataWM] = rescalePathData(DataInGM,DataInWM)


rescaledDataGM=DataInGM;
rescaledDataWM=DataInWM;

offset=0.00015;

maxVal = max(max(DataInGM(:)),max(DataInWM(:)));
minVal = min(min(DataInGM(:)),min(DataInWM(:)));

%This rescale Tau and TDP data to be max/min of their GM+WM Combined

rescaledDataGM=(DataInGM-minVal)/maxVal;
rescaledDataWM=(DataInWM-minVal)/maxVal;
rescaledDataGM=log(rescaledDataGM+offset);
rescaledDataWM=log(rescaledDataWM+offset);

