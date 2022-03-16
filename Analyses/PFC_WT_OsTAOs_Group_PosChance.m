% PFC_WellTrained_Group_OsTAOs_PositionChance
%%
% fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\HC\1. Well-Trained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Stella'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Mitt'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Barat'}];
% tets = [1,22,17,18,17];

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
trlWindow = {[-800 500]; [-500 800]};
alignment = [{'PokeIn'}, {'PokeOut'}];
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

numChancePerms = 100;

postCLim = [0 0.05];
decodeCLim = [0 0.2];


%%
for ani = 1:length(fileDirs)
    %% Create & setup initial object and data variables (if initial file)
    mlb = MLB_SM(fileDirs{ani});
    % Create Analysis Variables
    if ani == 1
        % Behavior Variables
        fiscRwdDelivLat = cell(length(fileDirs),1);
        smi = nan(length(fileDirs),1);
        dPrm = nan(length(fileDirs),1);
        ri = nan(length(fileDirs),1);
        smiByOP = nan(length(fileDirs),mlb.seqLength,2);
        dPrmByOP = nan(length(fileDirs),mlb.seqLength,2);
        riByOP = nan(length(fileDirs),mlb.seqLength,2);
        % Posteriors
        realPost_TransMat = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
        chancePost_TransMat = cell(mlb.seqLength, mlb.seqLength, numChancePerms, length(fileDirs));
        realPost_TAO = cell(mlb.seqLength, mlb.seqLength, length(fileDirs));
        chancePost_TAO = cell(mlb.seqLength, mlb.seqLength, numChancePerms, length(fileDirs));
    end
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
    fiscRwdDelivLat{ani} = ([mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).RewardIndex] - [mlb.trialInfo(mlb.fiscTrials(~isnan(mlb.fiscTrials))).PokeOutIndex])'/1000;
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
                for prevOdr = 1:mlb.seqLength
                    if prevOdr~=pos-1 && prevOdr~=pos                                   % This second condition removes skip dupes from the TAO analysis... may be losing significant power?
                        posTAO = find(curTrlLog);
                        prevOdrLog = trlOdrVect(posTAO-1)==prevOdr;
                        tempTaoPosts{prevOdr,pos-1} = mlb.post(:,:,posTAO(prevOdrLog));
                    end
                end
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
                    for prevOdr = 1:mlb.seqLength
                        if prevOdr~=pos-1
                            posTAO = find(curTrlLog);
                            prevOdrLog = trlOdrVect(posTAO-1)==prevOdr;
                            tempTaoPosts{prevOdr,pos-1} = mlb.post(:,:,posTAO(prevOdrLog));
                        end
                    end
                end
            end
        end
        chancePost_TransMat(:,:,perm,ani) = tempTmatPosts;
        chancePost_TAO(:,:,perm,ani) = tempTaoPosts;
    end
end
alignBoundry = find(mlb.decodeIDvects(:,2)==2,1,'first');
alignBoundTime = mlb.decodeIDvects(alignBoundry,1);
nearestPOtime = mlb.decodeIDvects(alignBoundry+(abs(trlWindow{2}(1))/dsRate),1);
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<=(nearestPOtime+median(rwdDelivLat)),1,'last'));
%% Collapse Trials Across Animals
% Posteriors
groupPostTIP_TransMat = cell(mlb.seqLength);        chancePostTIP_TransMat = cell(mlb.seqLength, mlb.seqLength, 2);
groupPostTime_TransMat = cell(mlb.seqLength);       chancePostTime_TransMat = cell(mlb.seqLength, mlb.seqLength, 2);
groupPostPos_TransMat = cell(mlb.seqLength);        chancePostPos_TransMat = cell(mlb.seqLength, mlb.seqLength, 2);
groupPostOvP_OutSeq = cell(1,1,mlb.seqLength^2);    chancePostOvP_OutSeq = cell(1,2,2);
groupPostTIP_TAO = cell(mlb.seqLength);             chancePostTIP_TAO = cell(mlb.seqLength, mlb.seqLength, 2);
groupPostTime_TAO = cell(mlb.seqLength);            chancePostTime_TAO = cell(mlb.seqLength, mlb.seqLength, 2);
groupPostPos_TAO = cell(mlb.seqLength);             chancePostPos_TAO = cell(mlb.seqLength, mlb.seqLength, 2);
groupPostOvP_TAO = cell(1,1,mlb.seqLength^2);       chancePostOvP_TAO = cell(1,2,3);
% Decodings
groupDecodeTIP_TransMat = cell(mlb.seqLength);      chanceDecodeTIP_TransMat = cell(mlb.seqLength, mlb.seqLength, 2);
groupDecodeTime_TransMat = cell(mlb.seqLength);     chanceDecodeTime_TransMat = cell(mlb.seqLength, mlb.seqLength, 2);
groupDecodePos_TransMat = cell(mlb.seqLength);      chanceDecodePos_TransMat = cell(mlb.seqLength, mlb.seqLength, 2);
groupDecodeOvP_OutSeq = cell(1,2);                  chanceDecodeOvP_OutSeq = cell(1,2,2);
groupDecodeTIP_TAO = cell(mlb.seqLength);           chanceDecodeTIP_TAO = cell(mlb.seqLength, mlb.seqLength, 2);
groupDecodeTime_TAO = cell(mlb.seqLength);          chanceDecodeTime_TAO = cell(mlb.seqLength, mlb.seqLength, 2);
groupDecodePos_TAO = cell(mlb.seqLength);           chanceDecodePos_TAO = cell(mlb.seqLength, mlb.seqLength, 2);
groupDecodeOvP_TAO = cell(1,3);                     chanceDecodeOvP_TAO = cell(1,2,3);

for odr = 1:size(realPost_TransMat,1)
    for pos = 1:size(realPost_TransMat,2)
        %% Observations
        %% Time-In-Position
        % Posteriors
        groupPostTIP_TransMat{odr,pos} = cell2mat(realPost_TransMat(odr,pos,:));
        groupPostTIP_TAO{odr,pos} = cell2mat(realPost_TAO(odr,pos,:));
        % Decodings
        tempDecodeTM = mlb.DecodeBayesPost(groupPostTIP_TransMat{odr,pos}, 1:size(mlb.decodeIDvects,1));
        tempDecodeTAO = mlb.DecodeBayesPost(groupPostTIP_TAO{odr,pos}, 1:size(mlb.decodeIDvects,1));
        groupDecodeTIP_TransMat{odr,pos} = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect));
        groupDecodeTIP_TAO{odr,pos} = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect));
        for ot = 1:length(mlb.obsvTimeVect)
            for lt = 1:length(mlb.likeTimeVect)
                if ~isempty(tempDecodeTM)
                    groupDecodeTIP_TransMat{odr,pos}(ot,lt) = mean(tempDecodeTM(ot,:)==lt);
                end
                if ~isempty(tempDecodeTAO)
                    groupDecodeTIP_TAO{odr,pos}(ot,lt) = mean(tempDecodeTAO(ot,:)==lt);
                end
            end
        end
        %% Time
        % Posteriors
        groupPostTime_TransMat{odr,pos} = mlb.TabulateBayesPost(groupPostTIP_TransMat{odr,pos}, mlb.decodeIDvects(:,1));
        groupPostTime_TAO{odr,pos} = mlb.TabulateBayesPost(groupPostTIP_TAO{odr,pos}, mlb.decodeIDvects(:,1));
        % Decodings
        tempTimeDecodeTM = mlb.DecodeBayesPost(groupPostTIP_TransMat{odr,pos}, mlb.decodeIDvects(:,1));
        tempTimeDecodeTAO = mlb.DecodeBayesPost(groupPostTIP_TAO{odr,pos}, mlb.decodeIDvects(:,1));
        groupDecodeTime_TransMat{odr,pos} = nan(length(mlb.obsvTimeVect));
        groupDecodeTime_TAO{odr,pos} = nan(length(mlb.obsvTimeVect));
        for ot = 1:length(mlb.obsvTimeVect)
            for lt = 1:length(mlb.obsvTimeVect)
                if ~isempty(tempTimeDecodeTM)
                    groupDecodeTime_TransMat{odr,pos}(ot,lt) = mean(tempTimeDecodeTM(ot,:)==mlb.obsvTimeVect(lt));
                end
                if ~isempty(tempTimeDecodeTAO)
                    groupDecodeTime_TAO{odr,pos}(ot,lt) = mean(tempTimeDecodeTAO(ot,:)==mlb.obsvTimeVect(lt));
                end
            end
        end
        %% Position
        % Posteriors
        groupPostPos_TransMat{odr,pos} = mlb.TabulateBayesPost(groupPostTIP_TransMat{odr,pos}, mlb.decodeIDvects(:,3));
        groupPostPos_TAO{odr, pos} = mlb.TabulateBayesPost(groupPostTIP_TAO{odr,pos}, mlb.decodeIDvects(:,3));
        % Decodings
        tempPosDecodeTM = mlb.DecodeBayesPost(groupPostTIP_TransMat{odr,pos}, mlb.decodeIDvects(:,3));
        tempPosDecodeTAO = mlb.DecodeBayesPost(groupPostTIP_TAO{odr,pos}, mlb.decodeIDvects(:,3));
        groupDecodePos_TransMat{odr,pos} = nan(length(mlb.obsvTimeVect), mlb.seqLength);
        groupDecodePos_TAO{odr,pos} = nan(length(mlb.obsvTimeVect), mlb.seqLength);
        for ot = 1:length(mlb.obsvTimeVect)
            for lt = 1:mlb.seqLength
                if ~isempty(tempPosDecodeTM)
                    groupDecodePos_TransMat{odr,pos}(ot,lt) = mean(tempPosDecodeTM(ot,:)==lt);
                end
                if ~isempty(tempPosDecodeTAO)
                    groupDecodePos_TAO{odr,pos}(ot,lt) = mean(tempPosDecodeTAO(ot,:)==lt);
                end
            end
        end
        %% Odor Vs Position
        if odr~=pos
            osNdx = sub2ind([mlb.seqLength, mlb.seqLength], odr,pos);
            %% OutSeq
            % Posteriors
            if ~isempty(groupPostPos_TransMat{odr,pos}) && odr~=pos+1                                                               % Second conditional here removes cases where the OutSeq odor is the next position in the sequence (i.e. Repeat-1) 
                groupPostOvP_OutSeq{osNdx} = [groupPostPos_TransMat{odr,pos}(:,pos,:), groupPostPos_TransMat{odr,pos}(:,odr,:), nan(length(mlb.obsvTimeVect),1, size(groupPostPos_TransMat{odr,pos},3))];
                if pos~=4
                    groupPostOvP_OutSeq{osNdx}(:,3,:) = groupPostPos_TransMat{odr,pos}(:,pos+1,:);
                end
            end
            if ~isempty(groupPostPos_TAO{odr,pos})
                groupPostOvP_TAO{osNdx} = [groupPostPos_TAO{odr,pos}(:,pos,:), groupPostPos_TAO{odr,pos}(:,odr,:), nan(length(mlb.obsvTimeVect),1, size(groupPostPos_TAO{odr,pos},3))];
                if pos~=4
                    groupPostOvP_TAO{osNdx}(:,3,:) = groupPostPos_TAO{odr,pos,:}(:,pos+1,:);
                end
            end
        end
        %% Chance
        tempPostTIP_ChanceTM = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), numChancePerms);
        tempPostTIP_ChanceTAO = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), numChancePerms);
        tempDecodeTIP_ChanceTM = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), numChancePerms);
        tempDecodeTIP_ChanceTAO = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), numChancePerms);
        
        tempPostTime_ChanceTM = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), numChancePerms);
        tempPostTime_ChanceTAO = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), numChancePerms);
        tempDecodeTime_ChanceTM = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), numChancePerms);
        tempDecodeTime_ChanceTAO = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), numChancePerms);
        
        tempPostPos_ChanceTM = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
        tempPostPos_ChanceTAO = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
        tempDecodePos_ChanceTM = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
        tempDecodePos_ChanceTAO = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
        for perm = 1:numChancePerms
            %% Time-In-Position
            % Posteriors
            curChanceTM = cell2mat(permute(chancePost_TransMat(odr,pos,perm,:), [1,2,4,3]));
            curChanceTAO = cell2mat(permute(chancePost_TAO(odr,pos,perm,:), [1,2,4,3]));
            if ~isempty(curChanceTM)
                tempPostTIP_ChanceTM(:,:,perm) = mean(curChanceTM,3);
            end
            if ~isempty(curChanceTAO)
                tempPostTIP_ChanceTAO(:,:,perm) = mean(curChanceTAO,3);
            end
            % Decodings
            tempChanceDecodeTM = mlb.DecodeBayesPost(curChanceTM, 1:size(mlb.decodeIDvects,1));
            tempChanceDecodeTAO = mlb.DecodeBayesPost(curChanceTAO, 1:size(mlb.decodeIDvects,1));
            for ot = 1:length(mlb.obsvTimeVect)
                for lt = 1:length(mlb.likeTimeVect)
                    if ~isempty(curChanceTM)
                        tempDecodeTIP_ChanceTM(ot,lt,perm) = mean(tempChanceDecodeTM(ot,:)==lt);
                    end
                    if ~isempty(curChanceTAO)
                        tempDecodeTIP_ChanceTAO(ot,lt,perm) = mean(tempChanceDecodeTAO(ot,:)==lt);
                    end
                end
            end
            %% Time
            % Posteriors
            curChanceTimeTM = mlb.TabulateBayesPost(curChanceTM, mlb.decodeIDvects(:,1));
            curChanceTimeTAO = mlb.TabulateBayesPost(curChanceTAO, mlb.decodeIDvects(:,1));
            if ~isempty(curChanceTimeTM)
                tempPostTime_ChanceTM(:,:,perm) = mean(curChanceTimeTM,3);
            end
            if ~isempty(curChanceTimeTAO)
                tempPostTime_ChanceTAO(:,:,perm) = mean(curChanceTimeTAO,3);
            end
            % Decodings
            tempChanceDecodeTimeTM = mlb.DecodeBayesPost(curChanceTimeTM, mlb.decodeIDvects(:,1));
            tempChanceDecodeTimeTAO = mlb.DecodeBayesPost(curChanceTimeTAO, mlb.decodeIDvects(:,1));
            for ot = 1:length(mlb.obsvTimeVect)
                for lt = 1:length(mlb.obsvTimeVect)
                    if ~isempty(tempChanceDecodeTimeTM)
                        tempDecodeTime_ChanceTM(ot,lt,perm) = mean(tempChanceDecodeTimeTM(ot,:)==mlb.obsvTimeVect(lt));
                    end
                    if ~isempty(tempChanceDecodeTimeTAO)
                        tempDecodeTime_ChanceTAO(ot,lt,perm) = mean(tempChanceDecodeTimeTAO(ot,:)==mlb.obsvTimeVect(lt));
                    end
                end
            end
            %% Position
            % Posteriors
            curChancePosTM = mlb.TabulateBayesPost(curChanceTM, mlb.decodeIDvects(:,3));
            curChancePosTAO = mlb.TabulateBayesPost(curChanceTAO, mlb.decodeIDvects(:,3));
            if ~isempty(curChancePosTM)
                tempPostPos_ChanceTM(:,:,perm) = mean(curChancePosTM,3);
            end
            if ~isempty(curChancePosTAO)
                tempPostPos_ChanceTAO(:,:,perm) = mean(curChancePosTAO,3);
            end
            % Decodings
            tempChanceDecodePosTM = mlb.DecodeBayesPost(curChanceTM, mlb.decodeIDvects(:,3));
            tempChanceDecodePosTAO = mlb.DecodeBayesPost(curChanceTAO, mlb.decodeIDvects(:,3));
            for ot = 1:length(mlb.obsvTimeVect)
                for lt = 1:mlb.seqLength
                    if ~isempty(curChancePosTM)
                        tempDecodePos_ChanceTM(ot,lt,perm) = mean(tempChanceDecodePosTM(ot,:)==lt);
                    end
                    if ~isempty(curChancePosTAO)
                        tempDecodePos_ChanceTAO(ot,lt,perm) = mean(tempChanceDecodePosTAO(ot,:)==lt);
                    end
                end
            end
        end
        chancePostTIP_TransMat{odr,pos,1} = mean(tempPostTIP_ChanceTM,3,'omitnan');
        chancePostTIP_TransMat{odr,pos,2} = std(tempPostTIP_ChanceTM,0,3,'omitnan');
        chanceDecodeTIP_TransMat{odr,pos,1} = mean(tempDecodeTIP_ChanceTM,3,'omitnan');
        chanceDecodeTIP_TransMat{odr,pos,2} = std(tempDecodeTIP_ChanceTM,0,3,'omitnan');
        
        chancePostTime_TransMat{odr,pos,1} = mean(tempPostTime_ChanceTM,3,'omitnan');
        chancePostTime_TransMat{odr,pos,2} = std(tempPostTime_ChanceTM,0,3,'omitnan');
        chanceDecodeTime_TransMat{odr,pos,1} = mean(tempDecodeTime_ChanceTM,3,'omitnan');
        chanceDecodeTime_TransMat{odr,pos,2} = std(tempDecodeTime_ChanceTM,0,3,'omitnan');
        
        chancePostPos_TransMat{odr,pos,1} = mean(tempPostPos_ChanceTM,3,'omitnan');
        chancePostPos_TransMat{odr,pos,2} = std(tempPostPos_ChanceTM,0,3,'omitnan');
        chanceDecodePos_TransMat{odr,pos,1} = mean(tempDecodePos_ChanceTM,3,'omitnan');
        chanceDecodePos_TransMat{odr,pos,2} = std(tempDecodePos_ChanceTM,0,3,'omitnan');
        
        chancePostTIP_TAO{odr,pos,1} = mean(tempPostTIP_ChanceTAO,3,'omitnan');
        chancePostTIP_TAO{odr,pos,2} = std(tempPostTIP_ChanceTAO,0,3,'omitnan');
        chanceDecodeTIP_TAO{odr,pos,1} = mean(tempDecodeTIP_ChanceTAO,3,'omitnan');
        chanceDecodeTIP_TAO{odr,pos,2} = std(tempDecodeTIP_ChanceTAO,0,3,'omitnan');
        
        chancePostTime_TAO{odr,pos,1} = mean(tempPostTime_ChanceTAO,3,'omitnan');
        chancePostTime_TAO{odr,pos,2} = std(tempPostTime_ChanceTAO,0,3,'omitnan');
        chanceDecodeTime_TAO{odr,pos,1} = mean(tempDecodeTime_ChanceTAO,3,'omitnan');
        chanceDecodeTime_TAO{odr,pos,2} = std(tempDecodeTime_ChanceTAO,0,3,'omitnan');
        
        chancePostPos_TAO{odr,pos,1} = mean(tempPostPos_ChanceTAO,3,'omitnan');
        chancePostPos_TAO{odr,pos,2} = std(tempPostPos_ChanceTAO,0,3,'omitnan');
        chanceDecodePos_TAO{odr,pos,1} = mean(tempDecodePos_ChanceTAO,3,'omitnan');
        chanceDecodePos_TAO{odr,pos,2} = std(tempDecodePos_ChanceTAO,0,3,'omitnan');
    end
end
%% Calculate chance levels for Odor vs Position Collapsing across trial types
tempPostOvP_OS = nan(length(mlb.obsvTimeVect), 2, numChancePerms);
tempPostOvP_TAO = nan(length(mlb.obsvTimeVect), 3, numChancePerms);

tempDecodeOvP_OS = nan(length(mlb.obsvTimeVect), 2, numChancePerms);
tempDecodeOvP_TAO = nan(length(mlb.obsvTimeVect), 3, numChancePerms);
for perm = 1:numChancePerms
    tempOS = cell(1,3,mlb.seqLength^2);
    tempTAO = cell(1,3,mlb.seqLength^2);
    for odr = 1:mlb.seqLength
        for pos = 1:mlb.seqLength
            if odr~=pos
                osNdx = sub2ind([mlb.seqLength, mlb.seqLength], odr,pos);
                curChanceTM = cell2mat(permute(chancePost_TransMat(odr,pos,perm,:), [1,2,4,3]));
                if ~isempty(curChanceTM) && odr~=pos+1                                                               % Second conditional here removes cases where the OutSeq trial is the next trial in the sequence (i.e. Skip-1) 
                    chanceTMposts = mlb.TabulateBayesPost(curChanceTM, mlb.decodeIDvects(:,3));
                    tempOS{1,1,osNdx} = chanceTMposts(:,odr,:);
                    tempOS{1,2,osNdx} = chanceTMposts(:,pos,:);
                    if pos~=4
                        tempOS{1,3,osNdx} = chanceTMposts(:,pos,:);
                    end
                end
                curChanceTAO = cell2mat(permute(chancePost_TAO(odr,pos,perm,:), [1,2,4,3]));
                if ~isempty(curChanceTAO)
                    chanceTAOposts = mlb.TabulateBayesPost(curChanceTAO, mlb.decodeIDvects(:,3));
                    tempTAO{1,1,osNdx} = chanceTAOposts(:,odr,:);
                    tempTAO{1,2,osNdx} = chanceTAOposts(:,pos,:);
                    if pos~=4
                        tempTAO{1,3,osNdx} = chanceTAOposts(:,pos+1,:);
                    end
                end
            end
        end
    end
    tempPostOvP_OS(:,1,perm) = mean(cell2mat(tempOS(1,1,:)),3);
    tempPostOvP_OS(:,2,perm) = mean(cell2mat(tempOS(1,2,:)),3);
    tempPostOvP_OS(:,3,perm) = mean(cell2mat(tempOS(1,3,:)),3);
    tempPostOvP_TAO(:,1,perm) = mean(cell2mat(tempTAO(1,1,:)),3);
    tempPostOvP_TAO(:,2,perm) = mean(cell2mat(tempTAO(1,2,:)),3);
    tempPostOvP_TAO(:,3,perm) = mean(cell2mat(tempTAO(1,3,:)),3);
    
    % I'm feeling lazy, I'm not going to bother with this right now.
%     tempDecodeOvP_OS = nan(length(mlb.obsvTimeVect), 2, numChancePerms);
%     tempDecodeOvP_TAO = nan(length(mlb.obsvTimeVect), 3, numChancePerms);
end
chancePostOvP_OutSeq{1,1,1} = mean(tempPostOvP_OS(:,1,:),3);
chancePostOvP_OutSeq{1,2,1} = std(tempPostOvP_OS(:,1,:),0,3);
chancePostOvP_OutSeq{1,1,2} = mean(tempPostOvP_OS(:,2,:),3);
chancePostOvP_OutSeq{1,2,2} = std(tempPostOvP_OS(:,2,:),0,3);
chancePostOvP_OutSeq{1,1,3} = mean(tempPostOvP_OS(:,3,:),3);
chancePostOvP_OutSeq{1,2,3} = std(tempPostOvP_OS(:,3,:),0,3);

chancePostOvP_TAO{1,1,1} = mean(tempPostOvP_TAO(:,1,:),3);
chancePostOvP_TAO{1,2,1} = std(tempPostOvP_TAO(:,1,:),0,3);
chancePostOvP_TAO{1,1,2} = mean(tempPostOvP_TAO(:,2,:),3);
chancePostOvP_TAO{1,2,2} = std(tempPostOvP_TAO(:,2,:),0,3);
chancePostOvP_TAO{1,1,3} = mean(tempPostOvP_TAO(:,3,:),3);
chancePostOvP_TAO{1,2,3} = std(tempPostOvP_TAO(:,3,:),0,3);

%%
clear chancePost_TAO chancePost_TransMat

%% Plot Probability & Decodings relative to chance
piNdx = find(mlb.likeTimeVect==0)+0.5;
poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;

figure;
%%%%%%%%%%%% OutSeq
subplot(2,3,1:2)
groupTrlOvP_OS = cell2mat(groupPostOvP_OutSeq);
posMean = mean(groupTrlOvP_OS(:,1,:),3, 'omitnan');
posSEM = mlb.SEMcalc(groupTrlOvP_OS(:,1,:),0,3);
posCI = tinv(0.975, size(groupTrlOvP_OS,3)-1).*posSEM;
posPlot = plot(mlb.obsvTimeVect, posMean, 'color', 'b', 'linewidth', 1.5);
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(posMean+posSEM)', flipud(posMean-posSEM)'],...
    'linestyle', 'none', 'facecolor', 'b', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(posMean+posCI)', flipud(posMean-posCI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'b', 'facealpha', 0);

odrMean = mean(groupTrlOvP_OS(:,2,:),3, 'omitnan');
odrSEM = mlb.SEMcalc(groupTrlOvP_OS(:,2,:),0,3);
odrCI = tinv(0.975, size(groupTrlOvP_OS,3)-1).*odrSEM;
odrPlot = plot(mlb.obsvTimeVect, odrMean, 'color', 'r', 'linewidth', 1.5);
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(odrMean+odrSEM)', flipud(odrMean-odrSEM)'],...
    'linestyle', 'none', 'facecolor', 'r', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(odrMean+odrCI)', flipud(odrMean-odrCI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'r', 'facealpha', 0);

nxtPosMean = mean(groupTrlOvP_OS(:,3,:),3, 'omitnan');
nxtPosSEM = mlb.SEMcalc(groupTrlOvP_OS(:,3,:),0,3);
nxtPosCI = tinv(0.975, size(groupTrlOvP_OS,3)-1).*nxtPosSEM;
nxtPosPlot = plot(mlb.obsvTimeVect, nxtPosMean, 'color', 'g', 'linewidth', 1.5);
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(nxtPosMean+nxtPosSEM)', flipud(nxtPosMean-nxtPosSEM)'],...
    'linestyle', 'none', 'facecolor', 'g', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(nxtPosMean+nxtPosCI)', flipud(nxtPosMean-nxtPosCI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'g', 'facealpha', 0);

chanceTrlOvP_OS = cell2mat(chancePostOvP_OutSeq);
chanceMean = chanceTrlOvP_OS(:,1,1);
chanceSEM = chanceTrlOvP_OS(:,1,2)./sqrt(numChancePerms-1);
chanceCI = tinv(0.975, numChancePerms-1).*chanceSEM;
chncPlot = plot(mlb.obsvTimeVect, chanceMean, 'color', 'k', 'linewidth', 1.5);
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(chanceMean+chanceSEM)', flipud(chanceMean-chanceSEM)'],...
    'linestyle', 'none', 'facecolor', 'k', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(chanceMean+chanceCI)', flipud(chanceMean-chanceCI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facealpha', 0);
axis tight;
set(gca, 'ylim', [0 0.75]);
plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
plot(repmat(alignBoundTime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 0.5);
set(gca, 'xtick', [0, alignBoundTime, nearestPOtime, nearestRWDtime], 'xticklabels', [{'PI'}, {'PI/PO boundary'}, {'PO'}, {'RWD'}]);

preTrlLog = mlb.obsvTimeVect>=0 & mlb.obsvTimeVect<=nearestPOtime;
patch('XData', [mlb.obsvTimeVect(find(preTrlLog,1,'first')), mlb.obsvTimeVect(find(preTrlLog,1,'last')),mlb.obsvTimeVect(find(preTrlLog,1,'last')),mlb.obsvTimeVect(find(preTrlLog,1,'first'))],...
    'YData', [0,0,repmat(max(get(gca, 'ylim')),[1,2])],...
    'linestyle', 'none', 'facecolor', 'k', 'facealpha', 0.1);

legend([posPlot, odrPlot, nxtPosPlot, chncPlot], [{'Position'}, {'Odor'}, {'Next Position'}, {'Chance'}], 'location','northwest');

subplot(2,3,3)
hold on;
posPreTrl = mean(groupTrlOvP_OS(preTrlLog,1,:));
posPreTrlMean = mean(posPreTrl);
posPreTrlSEM = mlb.SEMcalc(posPreTrl,0,3);
posPreTrlCI = tinv(0.975, size(posPreTrl,3)-1).*posPreTrlSEM;
bar(1, posPreTrlMean, 'b');
errorbar(1,posPreTrlMean, posPreTrlSEM, posPreTrlSEM, 'color', 'k', 'capsize', 0);
scatter(normrnd(0,0.1, [numel(posPreTrl),1])+1, posPreTrl(:), 'ok', 'filled', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'markerfacealpha', 0.2);

odrPreTrl = mean(groupTrlOvP_OS(preTrlLog,2,:));
odrPreTrlMean = mean(odrPreTrl);
odrPreTrlSEM = mlb.SEMcalc(odrPreTrl,0,3);
odrPreTrlCI = tinv(0.975, size(odrPreTrl,3)-1).*odrPreTrlSEM;
bar(2, odrPreTrlMean, 'r');
errorbar(2,odrPreTrlMean, odrPreTrlSEM, odrPreTrlSEM, 'color', 'k', 'capsize', 0);
scatter(normrnd(0,0.1, [numel(odrPreTrl),1])+2, odrPreTrl(:), 'ok', 'filled', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markerfacealpha', 0.2);

nxtPosPreTrl = mean(groupTrlOvP_OS(preTrlLog,3,:));
nxtPosPreTrlMean = mean(nxtPosPreTrl, 'omitnan');
nxtPosPreTrlSEM = mlb.SEMcalc(nxtPosPreTrl,0,3);
nxtPosPreTrlCI = tinv(0.975, size(nxtPosPreTrl,3)-1).*nxtPosPreTrlSEM;
bar(3, nxtPosPreTrlMean, 'g');
errorbar(3,nxtPosPreTrlMean, nxtPosPreTrlSEM, nxtPosPreTrlSEM, 'color', 'k', 'capsize', 0);
scatter(normrnd(0,0.1, [numel(nxtPosPreTrl),1])+3, nxtPosPreTrl(:), 'ok', 'filled', 'markerfacecolor', 'g', 'markeredgecolor', 'k', 'markerfacealpha', 0.2);

[pOS,tblOS,statsOS] = anova1([posPreTrl(:);odrPreTrl(:);nxtPosPreTrl(:)], [ones(numel(posPreTrl),1); ones(numel(odrPreTrl),1)*2; ones(numel(nxtPosPreTrl),1)*3],'off');
[h,p,ci,stats] = ttest(posPreTrl, odrPreTrl);
% bar


%%%%%%%%%%%%% TAO
subplot(2,3,4:5)
groupTrlOvP_TAO = cell2mat(groupPostOvP_TAO);
posMean = mean(groupTrlOvP_TAO(:,1,:), 3, 'omitnan');
posSEM = mlb.SEMcalc(groupTrlOvP_TAO(:,1,:),0,3);
posCI = tinv(0.975, size(groupTrlOvP_TAO,3)-1).*posSEM;
posPlot = plot(mlb.obsvTimeVect, posMean, 'color', 'b', 'linewidth', 1.5);
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(posMean+posSEM)', flipud(posMean-posSEM)'],...
    'linestyle', 'none', 'facecolor', 'b', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(posMean+posCI)', flipud(posMean-posCI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'b', 'facealpha', 0);

odrMean = mean(groupTrlOvP_TAO(:,2,:),3, 'omitnan');
odrSEM = mlb.SEMcalc(groupTrlOvP_TAO(:,2,:),0,3);
odrCI = tinv(0.975, size(groupTrlOvP_TAO,3)-1).*odrSEM;
odrPlot = plot(mlb.obsvTimeVect, odrMean, 'color', 'r', 'linewidth', 1.5);
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(odrMean+odrSEM)', flipud(odrMean-odrSEM)'],...
    'linestyle', 'none', 'facecolor', 'r', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(odrMean+odrCI)', flipud(odrMean-odrCI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'r', 'facealpha', 0);

nxtPosMean = mean(groupTrlOvP_TAO(:,3,:),3, 'omitnan');
nxtPosSEM = mlb.SEMcalc(groupTrlOvP_TAO(:,3,:),0,3);
nxtPosCI = tinv(0.975, size(groupTrlOvP_TAO,3)-1).*nxtPosSEM;
nxtPosPlot = plot(mlb.obsvTimeVect, nxtPosMean, 'color', 'g', 'linewidth', 1.5);
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(nxtPosMean+nxtPosSEM)', flipud(nxtPosMean-nxtPosSEM)'],...
    'linestyle', 'none', 'facecolor', 'g', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(nxtPosMean+nxtPosCI)', flipud(nxtPosMean-nxtPosCI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'g', 'facealpha', 0);

chanceTrlOvP_TAO = cell2mat(chancePostOvP_TAO);
chanceMean = chanceTrlOvP_TAO(:,1,1);
chanceSEM = chanceTrlOvP_TAO(:,1,2)./sqrt(numChancePerms-1);
chanceCI = tinv(0.975, numChancePerms-1).*chanceSEM;
chncPlot = plot(mlb.obsvTimeVect, chanceMean, 'color', 'k', 'linewidth', 1.5);
hold on;
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(chanceMean+chanceSEM)', flipud(chanceMean-chanceSEM)'],...
    'linestyle', 'none', 'facecolor', 'k', 'facealpha', 0.25);
patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
    'YData', [(chanceMean+chanceCI)', flipud(chanceMean-chanceCI)'],...
    'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facealpha', 0);
axis tight;
set(gca, 'ylim', [0 0.75]);
plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
plot(repmat(alignBoundTime, [1,2]), get(gca, 'ylim'), ':k', 'linewidth', 0.5);
set(gca, 'xtick', [0, alignBoundTime, nearestPOtime, nearestRWDtime], 'xticklabels', [{'PI'}, {'PI/PO boundary'}, {'PO'}, {'RWD'}]);


preTrlLog = mlb.obsvTimeVect<=0;
patch('XData', [mlb.obsvTimeVect(find(preTrlLog,1,'first')), mlb.obsvTimeVect(find(preTrlLog,1,'last')),mlb.obsvTimeVect(find(preTrlLog,1,'last')),mlb.obsvTimeVect(find(preTrlLog,1,'first'))],...
    'YData', [0,0,repmat(max(get(gca, 'ylim')),[1,2])],...
    'linestyle', 'none', 'facecolor', 'k', 'facealpha', 0.1);

legend([posPlot, odrPlot, nxtPosPlot, chncPlot], [{'Prev-Position'}, {'Prev-Odor'}, {'Current Position'}, {'Chance'}], 'location','northwest');

subplot(2,3,6)
hold on;
posPreTrl = mean(groupTrlOvP_TAO(preTrlLog,1,:));
posPreTrlMean = mean(posPreTrl);
posPreTrlSEM = mlb.SEMcalc(posPreTrl,0,3);
posPreTrlCI = tinv(0.975, size(posPreTrl,3)-1).*posPreTrlSEM;
bar(1, posPreTrlMean, 'b');
errorbar(1,posPreTrlMean, posPreTrlSEM, posPreTrlSEM, 'color', 'k', 'capsize', 0);
scatter(normrnd(0,0.1, [numel(posPreTrl),1])+1, posPreTrl(:), 'ok', 'filled', 'markerfacecolor', 'b', 'markeredgecolor', 'k', 'markerfacealpha', 0.2);

odrPreTrl = mean(groupTrlOvP_TAO(preTrlLog,2,:));
odrPreTrlMean = mean(odrPreTrl);
odrPreTrlSEM = mlb.SEMcalc(odrPreTrl,0,3);
odrPreTrlCI = tinv(0.975, size(odrPreTrl,3)-1).*odrPreTrlSEM;
bar(2, odrPreTrlMean, 'r');
errorbar(2,odrPreTrlMean, odrPreTrlSEM, odrPreTrlSEM, 'color', 'k', 'capsize', 0);
scatter(normrnd(0,0.1, [numel(odrPreTrl),1])+2, odrPreTrl(:), 'ok', 'filled', 'markerfacecolor', 'r', 'markeredgecolor', 'k', 'markerfacealpha', 0.2);

nxtPosPreTrl = mean(groupTrlOvP_TAO(preTrlLog,3,:));
nxtPosPreTrlMean = mean(nxtPosPreTrl, 'omitnan');
nxtPosPreTrlSEM = mlb.SEMcalc(nxtPosPreTrl,0,3);
nxtPosPreTrlCI = tinv(0.975, size(nxtPosPreTrl,3)-1).*nxtPosPreTrlSEM;
bar(3, nxtPosPreTrlMean, 'g');
errorbar(3,nxtPosPreTrlMean, nxtPosPreTrlSEM, nxtPosPreTrlSEM, 'color', 'k', 'capsize', 0);
scatter(normrnd(0,0.1, [numel(nxtPosPreTrl),1])+3, nxtPosPreTrl(:), 'ok', 'filled', 'markerfacecolor', 'g', 'markeredgecolor', 'k', 'markerfacealpha', 0.2);
[pTAO,tblTAO,statsTAO] = anova1([posPreTrl(:);odrPreTrl(:);nxtPosPreTrl(:)], [ones(numel(posPreTrl),1); ones(numel(odrPreTrl),1)*2; ones(numel(nxtPosPreTrl),1)*3], 'off');

%%
clear chancePost_TransMat chancePost_TAO
save('PFC_OsTAOs_PosChance.mat', '-v7.3');
% save('PFC_WT_OsTAOs_Group_Posts.mat', 'groupPost_TransMat', 'groupPost_TAO', 'groupChancePost_TransMat', 'groupChancePost_TAO', 'mlb', '-v7.3');