% PFC_WellTrained_FullTrial_Group
%%
fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
setupSeqLength = 4; 

% fileDirs = [{'D:\WorkBigDataFiles\HC\1. Well-Trained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Stella'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Mitt'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Barat'}];
% tets = [1,22,17,18,17];
% setupSeqLength = 5;

% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% setupSeqLength = 4; 
% 
% fileDirs = [
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
% setupSeqLength = 4;  

binSize = 200;
dsRate = 50;
trlWindow = {[-1000 2000]};
alignment = {'PokeIn'};
% trlWindow = {[-2000 800]};
% alignment = {'PokeOut'};
lfpWindow = [16 32];
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
smiByOP = nan(length(fileDirs),setupSeqLength,2);
dPrmByOP = nan(length(fileDirs),setupSeqLength,2);
riByOP = nan(length(fileDirs),setupSeqLength,2);
% fiscLOO analysis
fiscL1O_Posts = cell(setupSeqLength,1,length(fileDirs));
fiscL1O_TimePosDecode = cell(setupSeqLength,length(fileDirs));
fiscL1O_MargDecodeTime = cell(setupSeqLength,length(fileDirs));
fiscL1O_MargDecodePos = cell(setupSeqLength,length(fileDirs));
% fiscISC analysis
fiscISC_Posts = cell(setupSeqLength,1,length(fileDirs));
fiscISC_TimePosDecode = cell(setupSeqLength,length(fileDirs));
fiscISC_MargDecodeTime = cell(setupSeqLength,length(fileDirs));
fiscISC_MargDecodePos = cell(setupSeqLength,length(fileDirs));
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
%     uniInfo = mlb.unitInfo;
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05; % only MODULATED cells
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))>0.05; % only NON-MODULATED cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.windows = trlWindow;
    mlb.alignments = alignment;
    mlb.bayesType = bayesType;
    
    %% Extract Behavioral Variables
    fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).PokeOutIndex] - [mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).PokeInIndex])'/1000;
    fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).PokeInIndex])'/1000;
    smi(ani) = mlb.smi;
    dPrm(ani) = mlb.dPrime;
    ri(ani) = mlb.ri;
    for op = 1:mlb.seqLength
        smiByOP(ani,:,1) = mlb.smiByPos;
        smiByOP(ani,:,2) = mlb.smiByOdr;
        dPrmByOP(ani,:,1) = mlb.dPrimeByPos;
        dPrmByOP(ani,:,2) = mlb.dPrimeByOdr;
        riByOP(ani,:,1) = mlb.riByPos;
        riByOP(ani,:,2) = mlb.riByOdr;
    end
    %% Decode FISC via Leave-1-Out
    mlb.SetLikes_FISC;
%     mlb.bayesType = 3; % Comment In to decode using Gaussian rather than Poisson
    mlb.Process_LikelyL1O;
    % Extract Posteriors
    fiscL1O_Posts(:,:,ani) = mlb.post;
    % Decode Index
    fiscL1O_TimePosDecode(:,ani) = mlb.DecodeBayesPost(mlb.post, (1:size(mlb.decodeIDvects,1))');
    % Decode Time
    fiscL1O_MargDecodeTime(:,ani) = mlb.DecodeBayesPost(mlb.post, mlb.decodeIDvects(:,1));
    % Decode Position
    fiscL1O_MargDecodePos(:,ani) = mlb.DecodeBayesPost(mlb.post, mlb.decodeIDvects(:,3)); 
      
    %% Decode ISC via FISC 
    mlb.Process_Observes;
    % Decode Index
    ndxDecode = mlb.DecodeBayesPost(mlb.post, (1:size(mlb.decodeIDvects,1))');
    % Decode Time
    timeDecode = mlb.DecodeBayesPost(mlb.post, mlb.decodeIDvects(:,1));
    % Decode Position
    posDecode = mlb.DecodeBayesPost(mlb.post, mlb.decodeIDvects(:,3));
    
    trlOdrVect = [mlb.trialInfo(mlb.postTrlIDs).Odor];
    trlPosVect = [mlb.trialInfo(mlb.postTrlIDs).Position];
    trlPerfVect = [mlb.trialInfo(mlb.postTrlIDs).Performance];
    tempPosts = cell(mlb.seqLength,1);
    tempTimePosDecode = cell(mlb.seqLength,1);
    tempMargTime = cell(mlb.seqLength,1);
    tempMargPos = cell(mlb.seqLength,1);
    for pos = 1:mlb.seqLength
        curISClog = trlOdrVect==pos & trlPosVect==pos & trlPerfVect==1;
        tempPosts{pos} = mlb.post(:,:,curISClog);
        tempTimePosDecode{pos} = ndxDecode(:,curISClog);
        tempMargTime{pos} = timeDecode(:,curISClog);
        tempMargPos{pos} = posDecode(:,curISClog);
    end
    fiscISC_Posts(:,:,ani) = tempPosts;
    fiscISC_TimePosDecode(:,ani) = tempTimePosDecode;
    fiscISC_MargDecodeTime(:,ani) = tempMargTime;
    fiscISC_MargDecodePos(:,ani) = tempMargPos;
end

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));

%% Collapse across rats and calculate decodings
%% FISC-L1O
fiscPosts = mean(cell2mat(fiscL1O_Posts),3);
fiscTP = cell2mat(fiscL1O_TimePosDecode);
fiscTPacc = nan(size(mlb.decodeIDvects,1));
for t = 1:size(fiscTP,1)
    curTime = fiscTP(t,:);
    for ndx = 1:size(mlb.decodeIDvects,1)
        fiscTPacc(t,ndx) = mean(curTime==ndx);
    end
end
fiscTime = cell2mat(fiscL1O_MargDecodeTime)
timeVect = unique(mlb.decodeIDvects(:,1))
fiscL1O_MargDecodePos
%% FISC-ISC
iscPosts = cell(mlb.seqLength,1);
for pos = 1:mlb.seqLength
    iscPosts{pos} = mean(cell2mat(fiscISC_Posts(pos,:,:)),3);
end

fiscISC_TimePosDecode
fiscISC_MargDecodeTime
fiscISC_MargDecodePos
