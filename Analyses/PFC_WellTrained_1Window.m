fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

binSize = 200;
dsRate = 50;
trlWindow = [-800 1200];
alignment = 'PokeIn';
lfpWindow = [16 32];

for ani = 1:length(fileDirs)
    %% Create initial object & extract trial info
    mlb = MLB_SM(fileDirs{ani});
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    [trlSpikes, trlTimeVect] = mlb.PP_TrialMatrix_Spiking(trlWindow, alignment);
%     [trlLFPphase, trlLFPpower] = mlb.PP_TrialMatrix_LFP(lfpWindow, trlWindow, alignment);
    mlb.PP_IdentifyFISCseqs;
    mlb.SummarizeSessionBehavior;
    
    %% Decode FISC via Leave-1-Out
    %% Decode ISC via FISC 
    %% Decode ISC via Sub-Sampling

