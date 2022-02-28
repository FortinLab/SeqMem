% PFC_WellTrained_Group_OsTAOs_PositionChance
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
trlWindow = {[-800 500]; [-500 800]};
alignment = [{'PokeIn'}, {'PokeOut'}];
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

numChancePerms = 100;

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
% Posteriors
% May need to include regular decoding to show things aren't different....
realPost_TransMat = cell(setupSeqLength, setupSeqLength, length(fileDirs));
chancePost_TransMat = cell(setupSeqLength, setupSeqLength, numChancePerms, length(fileDirs));
realPost_TAO = cell(setupSeqLength, setupSeqLength, length(fileDirs));
chancePost_TAO = cell(setupSeqLength, setupSeqLength, numChancePerms, length(fileDirs));

%%
initTic = tic;
for ani = 1:length(fileDirs)
    tic;
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
    %% Analyze things
    mlb.SetLikes_FISC;
    mlb.Process_Observes;
    trlOdrVect = [mlb.trialInfo(mlb.postTrlIDs).Odor];
    trlPosVect = [mlb.trialInfo(mlb.postTrlIDs).Position];
    trlPerfVect = [mlb.trialInfo(mlb.postTrlIDs).Performance];
    tempTmatPosts = cell(mlb.seqLength);
    tempTaoPosts = cell(mlb.seqLength);
    for odr = 1:mlb.seqLength
        for pos = 1:mlb.seqLength
            curTrlLog = trlOdrVect==odr & trlPosVect==pos & trlPerfVect==1;
            tempTmatPosts{odr,pos} = mlb.post(:,:,curTrlLog);
            if pos~=1 && odr==pos
                posTAO = find(curTrlLog);
                taoLog = trlOdrVect(posTAO-1)~=trlPosVect(posTAO-1);
                tempTaoPosts{odr,pos} = mlb.post(:,:,posTAO(taoLog));
            end
        end
    end
    realPost_TransMat(:,:,ani) = tempTmatPosts;
    realPost_TAO(:,:,ani) = tempTaoPosts;
    %% Calculate Chance
    for perm = 1:numChancePerms
        chancePerm = permute(reshape(randperm(numel(mlb.fiscTrials)), size(mlb.fiscTrials)), [1,3,2]);
        mlb.likeTrlSpikes = mlb.likeTrlSpikes(chancePerm);
        mlb.Process_Observes;
        trlOdrVect = [mlb.trialInfo(mlb.postTrlIDs).Odor];
        trlPosVect = [mlb.trialInfo(mlb.postTrlIDs).Position];
        trlPerfVect = [mlb.trialInfo(mlb.postTrlIDs).Performance];
        tempTmatPosts = cell(mlb.seqLength,1);
        tempTaoPosts = cell(mlb.seqLength);
        for odr = 1:mlb.seqLength
            for pos = 1:mlb.seqLength
                curTrlLog = trlOdrVect==odr & trlPosVect==pos & trlPerfVect==1;
                tempTmatPosts{odr,pos} = mlb.post(:,:,curTrlLog);
                if pos~=1 && odr==pos
                    posTAO = find(curTrlLog);
                    taoLog = trlOdrVect(posTAO-1)~=trlPosVect(posTAO-1);
                    tempTaoPosts{odr,pos} = mlb.post(:,:,posTAO(taoLog));
                end
            end
        end
        chancePost_TransMat(:,:,perm,ani) = tempTmatPosts;
        chancePost_TAO(:,:,perm,ani) = tempTaoPosts;
    end
    toc;
end
toc(initTic);

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));
%% Collapse Trials Across Animals
groupPost_TransMat = cell(size(realPost_TransMat,1), size(realPost_TransMat,2));
groupPost_TAO = cell(size(realPost_TAO,1), size(realPost_TAO,2));
groupChancePost_TransMat = cell(size(chancePost_TransMat,1), size(chancePost_TransMat,2), size(chancePost_TransMat,3));
groupChancePost_TAO = cell(size(chancePost_TAO,1), size(chancePost_TAO,2), size(chancePost_TAO,3));
for odr = 1:mlb.seqLength
    for pos = 1:mlb.seqLength
        groupPost_TransMat{odr,pos} = cell2mat(realPost_TransMat(odr,pos,:));
        groupPost_TAO{odr,pos} = cell2mat(realPost_TAO(odr,pos,:));
        for perm = 1:numChancePerms
            groupChancePost_TransMat{odr,pos,perm} = cell2mat(permute(chancePost_TransMat(odr,pos,perm,:), [1,2,4,3]));
            groupChancePost_TAO{odr,pos,perm} = cell2mat(permute(chancePost_TAO(odr,pos,perm,:), [1,2,4,3]));
        end
    end
end
clear realPost_TransMat realPost_TAO chancePost_TransMat chancePost_TAO
save('PFC_WT_OsTAOs_Group_Posts.mat', 'groupPost_TransMat', 'groupPost_TAO', 'groupChancePost_TransMat', 'groupChancePost_TAO', 'mlb', '-v7.3');