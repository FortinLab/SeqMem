% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];
% setupSeqLength = 4; 

% % CA1 Data
% fileDirs = [{'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Stella'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Mitt'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\Barat'}];
% tets = [1,22,17,18,17]; % Lateral/Distal
% % tets = [7,3,1,5,5]; % Medial/Proximal
% setupSeqLength = 5;
% 
fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];
setupSeqLength = 4; 
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
ssProportion = 0.5;
numPerms = 100;
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
smiByOP = nan(length(fileDirs),setupSeqLength,2);
dPrmByOP = nan(length(fileDirs),setupSeqLength,2);
riByOP = nan(length(fileDirs),setupSeqLength,2);
% Neural Variables
grpDprmDecodePosts = cell(1,length(fileDirs));
grpDecodes = cell(1,length(fileDirs));
grpDecodeTrlNums = cell(1,length(fileDirs));
grpDecodeTrlPosIDs = cell(1,length(fileDirs));
grpDecodePermIDs = cell(1,length(fileDirs));
grpLFP = cell(1,length(fileDirs));
%%
for ani = 1:length(fileDirs)
    %% Create initial object & set object parameters
    mlb = MLB_SM(fileDirs{ani});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENT IN TO ONLY RUN ON BETA MODULATED CELLS & COMMENT OUT TO RUN ON ALL CELLS
%     uniInfo = mlb.unitInfo;
%     betaModCells(ani) = mean(cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05);
%     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))<0.05; % only MODULATED cells
% %     mlb.popVectIncludeLog = cell2mat(cellfun(@(a){a.Beta.R_Test(1)}, arrayfun(@(a){a.Spike_Phase_Relations}, mlb.unitInfo)))>0.05; % only NON-MODULATED cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mlb.binSize = binSize;
    mlb.dsRate = dsRate;
    mlb.windows = trlWindow;
    mlb.alignments = alignment;
    mlb.numPerms = 1;
    mlb.ssProportion = ssProportion;
    mlb.ssType = ssType;
    mlb.bayesType = bayesType;
    %% Extract Behavioral Variables
    if strcmp(alignment{1}, 'PokeIn')
        fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).PokeOutIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
        fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeInIndex])'/1000;
    elseif strcmp(alignment{1}, 'PokeOut')
        fiscPokeOutLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).PokeInIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeOutIndex])'/1000;
        fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials).PokeOutIndex])'/1000;
    end
       
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
    %% Process Neural Data
    [betaPhase, betaPower] = mlb.PP_TrialMatrix_LFP([16 32], trlWindow{1}, alignment{1});
    [thetaPhase, thetaPower] = mlb.PP_TrialMatrix_LFP([4 12], trlWindow{1}, alignment{1});
    
    dPrmPerm = cell(1,numPerms);
    tempDecodes = cell(1,1,numPerms);
    tempTrlNums = cell(1,numPerms);
    tempTrlIDs = cell(1,numPerms);
    tempPermIDs = cell(1,numPerms);
    for perm = 1:numPerms
        fprintf('Iteration #%i', perm);
        mlb.SetLikes_SubSample;
        mlb.Process_IterativeObserves;
        [decodes, ~] = mlb.DecodeBayesPost(mlb.post{1});
        dPrmPerm{perm} = mlb.CalcDprmMtxFromDecode(decodes, [mlb.trialInfo(mlb.obsvTrlIDs{1}).Position]);
        tempDecodes{perm} = decodes;
        tempTrlNums{perm} = squeeze(mlb.postTrlIDs{1})';
        tempTrlIDs{perm} = [mlb.trialInfo(mlb.postTrlIDs{1}).Position];
        tempPermIDs{perm} = squeeze(ones(size(mlb.postTrlIDs{1})))'*perm;
        fprintf(' complete\n');
    end
    grpDprmDecodePosts{ani} = dPrmPerm;
    grpDecodes{ani}	= cell2mat(tempDecodes);
    grpDecodeTrlNums{ani} = cell2mat(tempTrlNums);
    grpDecodeTrlPosIDs{ani} = cell2mat(tempTrlIDs);
    grpDecodePermIDs{ani} = cell2mat(tempPermIDs);
    grpLFP{ani} = {[betaPower, thetaPower]};
    
    %%
    figure; 
    for t = 1:4
        subplot(2,2,t); 
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect, mean(cell2mat(permute(cellfun(@(a){a(:,:,t)},dPrmPerm), [3,1,2])),3), [0 1]); 
        set(gca, 'ydir', 'normal'); 
        hold on;
        plot([0 0], get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), [0 0], '-k');
        plot(repmat(mean(fiscPokeOutLat{ani}), [1,2]), get(gca, 'ylim'), '-k');
        plot(get(gca, 'xlim'), repmat(mean(mean(fiscPokeOutLat{ani})), [1,2]), '-k');
        plot(repmat(mean(fiscRwdDelivLat{ani}), [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 2);
        plot(get(gca, 'xlim'), repmat(mean(mean(fiscRwdDelivLat{ani})), [1,2]), ':k', 'linewidth', 2);
        title(t);
    end
    annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
        'String', mlb.pathDir,...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    drawnow;
end    

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));