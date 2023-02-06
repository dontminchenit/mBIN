
Tinfo = 'D:\Min\Dropbox (Personal)\Research\Projects\FTDAnalysis\Data\Working copy path network sheet streamlines_finalized empiric algorithm detection data2_24_2020.xlsx';
Tinfo2 = 'Finalized_FTLD-Tau_TDP_CoreExtended_Threshold_PrelimLibLongWideCleaned_10_29_20_long.xlsx';

TTest = 'D:\Min\Dropbox (Personal)\Research\Projects\FTDAnalysis\Data\Finalized_FTLD-Tau_TDP_Threshold_PrelimLibLongWideCleaned_subset_pickedCases_hemiMatch_10_16_20_added demographics for MICE.xlsx';

TIn = readtable('./Finalized_FTLD-Tau_TDP_CoreExtended_Threshold_PrelimLibLongWideCleaned_10_29_20.xlsx');
T1=readtable(Tinfo);
T2=readtable(Tinfo2);

%TOut = agumentTableWith(TIn,Tinfo,'INDDID','CNDRDx',1);
%TOut = agumentTableWith(TOut,Tinfo,'INDDID','FTDDx',1);
%TOut = agumentTableWith(TOut,Tinfo2,'INDDID',T2.Properties.VariableNames{12},2);


TOut = agumentTableWith(TIn,TTest,'INDDID','MFC_GM',2);
TOut = agumentTableWith(TOut,TTest,'INDDID','ANG_GM',2);
TOut = agumentTableWith(TOut,TTest,'INDDID','STC_GM',2);
TOut = agumentTableWith(TOut,TTest,'INDDID','OFC_GM',2);
TOut = agumentTableWith(TOut,TTest,'INDDID','Tau1_TDP2',2);
TOut = agumentTableWith(TOut,TTest,'INDDID','Hemisphere',1);

writetable(TOut,'temp.xlsx')