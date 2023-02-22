%winBase = fullfile('D:','Min','Dropbox (Personal)','Research','Projects','PhasingAnalysis');
%macBase = fullfile('/Users','Min','Dropbox (Personal)','Research','Projects','PhasingAnalysis');
%currBase = winBase;

%inInfo = fullfile(currBase,'Data','subj95_baseline.xlsx');
%T = readtable(inInfo);

%inInfo = fullfile(currBase,'Data','dti_subs_5subRemoved.csv');
%TPhases = readtable(inInfo);


inInfoMRI = fullfile(currBase,'PhaseData','demog_withDuration_phaseNetworks.csv');

TMRI = readtable(inInfoMRI);

%inInfo = 'D:\Min\Dropbox (Personal)\Research\Projects\FTDAnalysis\Data\Final Network Path Sheet 5_7_21.csv';
inInfoMRI = fullfile(currBase,'PhaseData','phase_network_demographics.csv');
TAO = readtable(inInfoMRI);

SarahDiseaseDurAO=nan(length(TPhases.INDDID),1);
SarahDiseaseDurMR=nan(length(TPhases.INDDID),1);
for i = 1:length(TPhases.INDDID)

    currID = TPhases.INDDID(i)
    ind=find(TMRI.INDDID==currID)
    if(~isempty(ind))
    %SarahDiseaseDurMR(i) = TDemo.DurationAtMRI(ind);
    end
    
    ind=find(TAO.INDDID==currID);
    if(~isempty(ind))
    SarahDiseaseDurAO(i) = TAO.GlobalDisDur_toAut(ind);
    end
    
    
    
end

SarahDiseaseDur = SarahDiseaseDurAO;