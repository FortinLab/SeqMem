% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
% 
fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% 
% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
binSize = 200;
dsRate = 50;
trlWindow = {[-1000 2000]};
alignment = {'PokeIn'};
% trlWindow = {[-2000 800]};
% alignment = {'PokeOut'};
lfpWindow = [16 32];
numPerms = 10;
ssProportion = 0.4;
ssType = 1; % 0 = use all ISC for decoding; 1 = use subsampled ISC types
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

postCLim = [0 0.05];
decodeCLim = [0 0.2];

%% Create Output Variables
% Behavior Variables
fiscPokeOutLat = cell(length(fileDirs),1);
fiscRwdDelivLat = cell(length(fileDirs),1);
smi = nan(length(fileDirs),1);
dPrm = nan(length(fileDirs),1);
ri = nan(length(fileDirs),1);
smiByOP = nan(length(fileDirs),4,2);
dPrmByOP = nan(length(fileDirs),4,2);
riByOP = nan(length(fileDirs),4,2);
% Neural Variables
posts = cell(size(fileDirs));
posDecodes = cell(size(fileDirs));
trlLFP = cell(3,1,length(fileDirs));
distMtx = cell(1,1,length(fileDirs));
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.windows = trlWindow;
    mlb.alignments = alignment;
    mlb.numPerms = numPerms;
    mlb.ssProportion = ssProportion;
    mlb.ssType = ssType;
    mlb.bayesType = bayesType;
    %% Extract Behavioral Variables
%     fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).PokeOutIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
%     fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
%     smi(ani) = mlb.smi;
%     dPrm(ani) = mlb.dPrime;
%     ri(ani) = mlb.ri;
%     for op = 1:mlb.seqLength
%         smiByOP(ani,:,1) = mlb.smiByPos;
%         smiByOP(ani,:,2) = mlb.smiByOdr;
%         dPrmByOP(ani,:,1) = mlb.dPrimeByPos;
%         dPrmByOP(ani,:,2) = mlb.dPrimeByOdr;
%         riByOP(ani,:,1) = mlb.riByPos;
%         riByOP(ani,:,2) = mlb.riByOdr;
%     end
%     nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(fiscPokeOutLat{ani}),1,'last'));
%     nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(fiscRwdDelivLat{ani}),1,'last'));
    %% Extract LFP
%     [betaPhase, betaPower] = mlb.PP_TrialMatrix_LFP([16 32], trlWindow{1}, alignment{1});
%     [thetaPhase, thetaPower] = mlb.PP_TrialMatrix_LFP([4 12], trlWindow{1}, alignment{1});
    %% Decode ISC via Sub-Sampling
%     mlb.bayesType = 1;  % Comment In to decode using Gaussian rather than Poisson
    mlb.SetLikes_SubSample;
    distMtx{ani} = mlb.PFC_validPCA;
end