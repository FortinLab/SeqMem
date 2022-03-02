% PFC_WellTrained_Group_PositionChance
%%
fileDirs = [{'D:\WorkBigDataFiles\PFC\Files To Process\GE11\GE11_Session132'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE13\GE13_Session083'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE14\GE14_Session123'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE17\GE17_Session095'},...
    {'D:\WorkBigDataFiles\PFC\Files To Process\GE24\Session096'}];

% fileDirs = [{'D:\WorkBigDataFiles\HC\1. Well-Trained session\SuperChris'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Stella'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Mitt'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Buchanan'},...
%     {'D:\WorkBigDataFiles\HC\1. Well-Trained session\Barat'}];
% tets = [1,22,17,18,17];
% 
% fileDirs = [{'D:\WorkBigDataFiles\PFC\GE11_Session132'},...
%     {'D:\WorkBigDataFiles\PFC\GE13_Session083'},...
%     {'D:\WorkBigDataFiles\PFC\GE14_Session123'},...
%     {'D:\WorkBigDataFiles\PFC\GE17_Session095'},...
%     {'D:\WorkBigDataFiles\PFC\GE24_Session096'}];

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
bayesType = 1; %1 = Poisson: use with raw spike counts; 2 = Bernoulli: use with binarized spike counts; 3 = Gaussian: Use with z-scored spike counts

numChancePerms = 100;

postCLim = [0 0.1];
decodeCLim = [0 0.2];

%%
for ani = 1:length(fileDirs)
    %% Create & setup initial object and data variables (if initial file)
    mlb = MLB_SM(fileDirs{ani});
    % Create Analysis Variables
    if ani == 1 
        % Behavior Variables
        fiscPokeOutLat = cell(length(fileDirs),1);
        fiscRwdDelivLat = cell(length(fileDirs),1);
        smi = nan(length(fileDirs),1);
        dPrm = nan(length(fileDirs),1);
        ri = nan(length(fileDirs),1);
        smiByOP = nan(length(fileDirs),mlb.seqLength,2);
        dPrmByOP = nan(length(fileDirs),mlb.seqLength,2);
        riByOP = nan(length(fileDirs),mlb.seqLength,2);
        % Posteriors
        realPost_L1O = cell(mlb.seqLength, 1, length(fileDirs));
        chancePost_L1O = cell(mlb.seqLength, numChancePerms, length(fileDirs));
        realPost_ISC = cell(mlb.seqLength, 1, length(fileDirs));
        chancePost_ISC= cell(mlb.seqLength, numChancePerms, length(fileDirs));
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
    %% Process Observations
    % FISC via Leave-1-Out
    mlb.SetLikes_FISC;
    mlb.Process_LikelyL1O;
    realPost_L1O(:,1,ani) = mlb.post;      
    % ISC via FISC 
    mlb.Process_Observes;
    trlOdrVect = [mlb.trialInfo(mlb.postTrlIDs).Odor];
    trlPosVect = [mlb.trialInfo(mlb.postTrlIDs).Position];
    trlPerfVect = [mlb.trialInfo(mlb.postTrlIDs).Performance];
    for pos = 1:mlb.seqLength
        curISClog = trlOdrVect==pos & trlPosVect==pos & trlPerfVect==1;
        realPost_ISC{pos,1,ani} = mlb.post(:,:,curISClog);
    end
    %% Process Chance
    for perm = 1:numChancePerms 
        chancePerm = permute(reshape(randperm(numel(mlb.fiscTrials)), size(mlb.fiscTrials)), [1,3,2]);
        mlb.likeTrlSpikes = mlb.likeTrlSpikes(chancePerm);
        mlb.Process_LikelyL1O;
        chancePost_L1O(:,perm,ani) = mlb.post;
        mlb.Process_Observes;
        trlOdrVect = [mlb.trialInfo(mlb.postTrlIDs).Odor];
        trlPosVect = [mlb.trialInfo(mlb.postTrlIDs).Position];
        trlPerfVect = [mlb.trialInfo(mlb.postTrlIDs).Performance];
        for pos = 1:mlb.seqLength
            curISClog = trlOdrVect==pos & trlPosVect==pos & trlPerfVect==1;
            chancePost_ISC{pos,perm,ani} = mlb.post(:,:,curISClog);
        end
    end
end

pokeOutLats = cell2mat(fiscPokeOutLat(:));
nearestPOtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(pokeOutLats),1,'last'));
rwdDelivLat = cell2mat(fiscRwdDelivLat(:));
nearestRWDtime = mlb.obsvTimeVect(find(mlb.obsvTimeVect<median(rwdDelivLat),1,'last'));

%% Collapse Trials Across Animals
groupPostTIP_RealL1O = cell(mlb.seqLength,1);           groupPostTIP_ChanceL1O = cell(mlb.seqLength,2);
groupDecodeTIP_RealL1O = cell(mlb.seqLength,1);         groupDecodeTIP_ChanceL1O = cell(mlb.seqLength,2);
groupPostTime_RealL1O = cell(mlb.seqLength,1);          groupPostTime_ChanceL1O = cell(mlb.seqLength,2);
groupDecodeTime_RealL1O = cell(mlb.seqLength,1);        groupDecodeTime_ChanceL1O = cell(mlb.seqLength,2);
groupPostPos_RealL1O = cell(mlb.seqLength,1);           groupPostPos_ChanceL1O = cell(mlb.seqLength,2);
groupDecodePos_RealL1O = cell(mlb.seqLength,1);         groupDecodePos_ChanceL1O = cell(mlb.seqLength,2);

groupPostTIP_RealISC = cell(mlb.seqLength,1);           groupPostTIP_ChanceISC = cell(mlb.seqLength,2);
groupDecodeTIP_RealISC = cell(mlb.seqLength,1);         groupDecodeTIP_ChanceISC = cell(mlb.seqLength,2);
groupPostTime_RealISC = cell(mlb.seqLength,1);          groupPostTime_ChanceISC = cell(mlb.seqLength,2);
groupDecodeTime_RealISC = cell(mlb.seqLength,1);        groupDecodeTime_ChanceISC = cell(mlb.seqLength,2);
groupPostPos_RealISC = cell(mlb.seqLength,1);           groupPostPos_ChanceISC = cell(mlb.seqLength,2);
groupDecodePos_RealISC = cell(mlb.seqLength,1);         groupDecodePos_ChanceISC = cell(mlb.seqLength,2);

for pos = 1:mlb.seqLength
%% Observations
    %% Time-In-Position
    % Posteriors
    groupPostTIP_RealL1O{pos} = cell2mat(realPost_L1O(pos,:,:));
    groupPostTIP_RealISC{pos} = cell2mat(realPost_ISC(pos,:,:));
    % Decodings
    tempDecodeL1O = mlb.DecodeBayesPost(groupPostTIP_RealL1O{pos}, 1:size(mlb.decodeIDvects,1));
    tempDecodeISC = mlb.DecodeBayesPost(groupPostTIP_RealISC{pos}, 1:size(mlb.decodeIDvects,1));
    tempTIP_L1O = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect));
    tempTIP_ISC = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect));
    for ot = 1:length(mlb.obsvTimeVect)
        for lt = 1:length(mlb.likeTimeVect)
            tempTIP_L1O(ot,lt) = mean(tempDecodeL1O(ot,:)==lt);
            tempTIP_ISC(ot,lt) = mean(tempDecodeISC(ot,:)==lt);
        end
    end
    groupDecodeTIP_RealL1O{pos} = tempTIP_L1O;
    groupDecodeTIP_RealISC{pos} = tempTIP_ISC;
    
    %% Time
    % Posteriors
    groupPostTime_RealL1O{pos} = mlb.TabulateBayesPost(groupPostTIP_RealL1O{pos}, mlb.decodeIDvects(:,1));
    groupPostTime_RealISC{pos} = mlb.TabulateBayesPost(groupPostTIP_RealISC{pos}, mlb.decodeIDvects(:,1));
    % Decodings
    tempTimeDecodeL1O = mlb.DecodeBayesPost(groupPostTIP_RealL1O{pos}, mlb.decodeIDvects(:,1));
    tempTimeDecodeISC = mlb.DecodeBayesPost(groupPostTIP_RealISC{pos}, mlb.decodeIDvects(:,1));
    tempTime_L1O = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    tempTime_ISC = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect));
    for ot = 1:length(mlb.obsvTimeVect)
        for lt = 1:length(mlb.obsvTimeVect)
            tempTime_L1O(ot,lt) = mean(tempTimeDecodeL1O(ot,:)==mlb.obsvTimeVect(lt));
            tempTime_ISC(ot,lt) = mean(tempTimeDecodeISC(ot,:)==mlb.obsvTimeVect(lt));
        end
    end
    groupDecodeTime_RealL1O{pos} = tempTime_L1O;
    groupDecodeTime_RealISC{pos} = tempTime_ISC;
            
    %% Position
    % Posteriors
    groupPostPos_RealL1O{pos} = mlb.TabulateBayesPost(groupPostTIP_RealL1O{pos}, mlb.decodeIDvects(:,3));
    groupPostPos_RealISC{pos} = mlb.TabulateBayesPost(groupPostTIP_RealISC{pos}, mlb.decodeIDvects(:,3));
    % Decodings
    tempPosDecodeL1O = mlb.DecodeBayesPost(groupPostTIP_RealL1O{pos}, mlb.decodeIDvects(:,3));
    tempPosDecodeISC = mlb.DecodeBayesPost(groupPostTIP_RealISC{pos}, mlb.decodeIDvects(:,3));
    tempPos_L1O = nan(length(mlb.obsvTimeVect), mlb.seqLength);
    tempPos_ISC = nan(length(mlb.obsvTimeVect), mlb.seqLength);
    for ot = 1:length(mlb.obsvTimeVect)
        for lt = 1:mlb.seqLength
            tempPos_L1O(ot,lt) = mean(tempPosDecodeL1O (ot,:)==lt);
            tempPos_ISC(ot,lt) = mean(tempPosDecodeISC(ot,:)==lt);
        end
    end
    groupDecodePos_RealL1O{pos} = tempPos_L1O;
    groupDecodePos_RealISC{pos} = tempPos_ISC;
    
%% Chance    
    tempPostTIP_ChanceL1O = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), numChancePerms);
    tempPostTIP_ChanceISC = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), numChancePerms);
    tempDecodeTIP_ChanceL1O = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), numChancePerms);
    tempDecodeTIP_ChanceISC = nan(length(mlb.obsvTimeVect), length(mlb.likeTimeVect), numChancePerms);
    
    tempPostTime_ChanceL1O = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), numChancePerms);
    tempPostTime_ChanceISC = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), numChancePerms);
    tempDecodeTime_ChanceL1O = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), numChancePerms);
    tempDecodeTime_ChanceISC = nan(length(mlb.obsvTimeVect), length(mlb.obsvTimeVect), numChancePerms);
    
    tempPostPos_ChanceL1O = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
    tempPostPos_ChanceISC = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
    tempDecodePos_ChanceL1O = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
    tempDecodePos_ChanceISC = nan(length(mlb.obsvTimeVect), mlb.seqLength, numChancePerms);
    for perm = 1:numChancePerms
        %% TIP
        % Posteriors
        curChanceL1O = cell2mat(chancePost_L1O(pos,perm,:));        
        curChanceISC = cell2mat(chancePost_ISC(pos,perm,:));
        tempPostTIP_ChanceL1O(:,:,perm) = mean(curChanceL1O,3);
        tempPostTIP_ChanceISC(:,:,perm) = mean(curChanceISC,3);    
        % Decodings
        tempDecodeL1O = mlb.DecodeBayesPost(curChanceL1O, 1:size(mlb.decodeIDvects,1));
        tempDecodeISC = mlb.DecodeBayesPost(curChanceISC, 1:size(mlb.decodeIDvects,1));
        for ot = 1:length(mlb.obsvTimeVect)
            for lt = 1:length(mlb.likeTimeVect)
                tempDecodeTIP_ChanceL1O(ot,lt,perm) = mean(tempDecodeL1O(ot,:)==lt);
                tempDecodeTIP_ChanceISC(ot,lt,perm) = mean(tempDecodeISC(ot,:)==lt);
            end
        end
        
        %% Time
        % Posteriors
        tempPostTime_ChanceL1O(:,:,perm) = mean(mlb.TabulateBayesPost(curChanceL1O, mlb.decodeIDvects(:,1)),3);
        tempPostTime_ChanceISC(:,:,perm) = mean(mlb.TabulateBayesPost(curChanceISC, mlb.decodeIDvects(:,1)),3);
        % Decodings
        tempDecodeTimeL1O = mlb.DecodeBayesPost(curChanceL1O, mlb.decodeIDvects(:,1));
        tempDecodeTimeISC = mlb.DecodeBayesPost(curChanceISC, mlb.decodeIDvects(:,1));
        for ot = 1:length(mlb.obsvTimeVect)
            for lt = 1:length(mlb.obsvTimeVect)
                tempDecodeTime_ChanceL1O(ot,lt,perm) = mean(tempDecodeTimeL1O(ot,:)==mlb.obsvTimeVect(lt));
                tempDecodeTime_ChanceISC(ot,lt,perm) = mean(tempDecodeTimeISC(ot,:)==mlb.obsvTimeVect(lt));
            end
        end
        
        %% Position
        % Posteriors
        tempPostPos_ChanceL1O(:,:,perm) = mean(mlb.TabulateBayesPost(curChanceL1O, mlb.decodeIDvects(:,3)),3);
        tempPostPos_ChanceISC(:,:,perm) = mean(mlb.TabulateBayesPost(curChanceISC, mlb.decodeIDvects(:,3)),3);
        % Decodings
        tempDecodePosL1O = mlb.DecodeBayesPost(curChanceL1O, mlb.decodeIDvects(:,3));
        tempDecodePosISC = mlb.DecodeBayesPost(curChanceISC, mlb.decodeIDvects(:,3));
        for ot = 1:length(mlb.obsvTimeVect)
            for lt = 1:mlb.seqLength
                tempDecodePos_ChanceL1O(ot,lt,perm) = mean(tempDecodePosL1O(ot,:)==lt);
                tempDecodePos_ChanceISC(ot,lt,perm) = mean(tempDecodePosL1O(ot,:)==lt);
            end
        end
    end
    % L1O
    groupPostTIP_ChanceL1O{pos,1} = mean(tempPostTIP_ChanceL1O, 3, 'omitnan');
    groupPostTIP_ChanceL1O{pos,2} = std(tempPostTIP_ChanceL1O, 0, 3, 'omitnan');
    groupDecodeTIP_ChanceL1O{pos,1} = mean(tempDecodeTIP_ChanceL1O, 3, 'omitnan');
    groupDecodeTIP_ChanceL1O{pos,2} = std(tempDecodeTIP_ChanceL1O, 0, 3, 'omitnan');
    
    groupPostTime_ChanceL1O{pos,1} = mean(tempPostTime_ChanceL1O, 3, 'omitnan');
    groupPostTime_ChanceL1O{pos,2} = std(tempPostTime_ChanceL1O, 0, 3, 'omitnan');
    groupDecodeTime_ChanceL1O{pos,1} = mean(tempDecodeTime_ChanceL1O, 3, 'omitnan');
    groupDecodeTime_ChanceL1O{pos,2} = std(tempDecodeTime_ChanceL1O, 0, 3, 'omitnan');
    
    groupPostPos_ChanceL1O{pos,1} = mean(tempPostPos_ChanceL1O, 3, 'omitnan');
    groupPostPos_ChanceL1O{pos,2} = std(tempPostPos_ChanceL1O, 0,  3, 'omitnan');
    groupDecodePos_ChanceL1O{pos,1} = mean(tempDecodePos_ChanceL1O, 3, 'omitnan');
    groupDecodePos_ChanceL1O{pos,2} = std(tempDecodePos_ChanceL1O, 0, 3, 'omitnan');
    
    %ISC
    groupPostTIP_ChanceISC{pos,1} = mean(tempPostTIP_ChanceISC, 3, 'omitnan');
    groupPostTIP_ChanceISC{pos,2} = std(tempPostTIP_ChanceISC, 0, 3, 'omitnan');
    groupDecodeTIP_ChanceISC{pos,1} = mean(tempDecodeTIP_ChanceISC, 3, 'omitnan');
    groupDecodeTIP_ChanceISC{pos,2} = std(tempDecodeTIP_ChanceISC, 0, 3, 'omitnan');
    
    groupPostTime_ChanceISC{pos,1} = mean(tempPostTime_ChanceISC, 3, 'omitnan');
    groupPostTime_ChanceISC{pos,2} = std(tempPostTime_ChanceISC, 0, 3, 'omitnan');
    groupDecodeTime_ChanceISC{pos,1} = mean(tempDecodeTime_ChanceISC, 3, 'omitnan');
    groupDecodeTime_ChanceISC{pos,2} = std(tempDecodeTime_ChanceISC, 0, 3, 'omitnan');
    
    groupPostPos_ChanceISC{pos,1} = mean(tempPostPos_ChanceISC, 3, 'omitnan');
    groupPostPos_ChanceISC{pos,2} = std(tempPostPos_ChanceISC, 0,  3, 'omitnan');
    groupDecodePos_ChanceISC{pos,1} = mean(tempDecodePos_ChanceISC, 3, 'omitnan');
    groupDecodePos_ChanceISC{pos,2} = std(tempDecodePos_ChanceISC, 0, 3, 'omitnan');
end
%%
clear realPost_L1O realPost_ISC chancePost_L1O chancePost_ISC
%% Plot Probability & Decodings relative to chance
piNdx = find(mlb.likeTimeVect==0)+0.5;
poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;

figure; 
chanceMean_L1O = cell(mlb.seqLength,1);
chanceCI_L1O = cell(mlb.seqLength,1);

chanceMean_ISC = cell(mlb.seqLength,1);
chanceCI_ISC = cell(mlb.seqLength,1);
for pos = 1:mlb.seqLength
    chanceMean_L1O{pos} = groupPostPos_ChanceL1O{pos,1}(:,pos);
    chanceCI_L1O{pos} = tinv(0.975, numChancePerms-1)*(groupPostPos_ChanceL1O{pos,2}(:,pos)./sqrt(numChancePerms));    
    tempPostMeanL1O = cell2mat(cellfun(@(a){mean(a(:,pos,:),3,'omitnan')}, groupPostPos_RealL1O));
    tempPostSEML1O = cell2mat(cellfun(@(a){std(a(:,pos,:),0,3)./sqrt(size(a,3))}, groupPostPos_RealL1O));
    tempPostCIL1O = cell2mat(cellfun(@(a){tinv(0.975, size(a,3)-1)*(std(a(:,pos,:),0,3)./sqrt(size(a,3)))}, groupPostPos_RealL1O));
    subplot(2,1,1);
    hold on;
    plot(tempPostMeanL1O, 'color', mlb.PositionColors(pos,:), 'linewidth', 1.5);
    patch('XData', [1:length(mlb.likeTimeVect), length(mlb.likeTimeVect):-1:1],...
        'YData', [(tempPostMeanL1O+tempPostSEML1O)', flipud(tempPostMeanL1O-tempPostSEML1O)'],...
        'linestyle', 'none', 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0.25);
    patch('XData', [1:length(mlb.likeTimeVect), length(mlb.likeTimeVect):-1:1],...
        'YData', [(tempPostMeanL1O+tempPostCIL1O)', flipud(tempPostMeanL1O-tempPostCIL1O)'],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0);
    
    chanceMean_ISC{pos} = groupPostPos_ChanceISC{pos,1}(:,pos);
    chanceCI_ISC{pos} = tinv(0.975, numChancePerms-1)*(groupPostPos_ChanceISC{pos,2}(:,pos)./sqrt(numChancePerms));
    tempPostMeanISC = cell2mat(cellfun(@(a){mean(a(:,pos,:),3,'omitnan')}, groupPostPos_RealISC));
    tempPostSEMISC = cell2mat(cellfun(@(a){std(a(:,pos,:),0,3)./sqrt(size(a,3))}, groupPostPos_RealISC));
    tempPostCIISC = cell2mat(cellfun(@(a){tinv(0.975, size(a,3)-1)*(std(a(:,pos,:),0,3)./sqrt(size(a,3)))}, groupPostPos_RealISC));
    subplot(2,1,2)
    hold on;
    plot(tempPostMeanISC, 'color', mlb.PositionColors(pos,:), 'linewidth', 1.5);
    patch('XData', [1:length(mlb.likeTimeVect), length(mlb.likeTimeVect):-1:1],...
        'YData', [(tempPostMeanISC+tempPostSEMISC)', flipud(tempPostMeanISC-tempPostSEMISC)'],...
        'linestyle', 'none', 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0.25);
    patch('XData', [1:length(mlb.likeTimeVect), length(mlb.likeTimeVect):-1:1],...
        'YData', [(tempPostMeanISC+tempPostCIISC)', flipud(tempPostMeanISC-tempPostCIISC)'],...
        'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(pos,:), 'facecolor', mlb.PositionColors(pos,:), 'facealpha', 0);
end
for sp = 1:2
    for ndx = 1:length(piNdx)
        subplot(2,1,sp)
        plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        
        if ndx<length(piNdx)
            plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
        end
    end
end
chanceMean_L1O = cell2mat(chanceMean_L1O);
chanceCI_L1O = cell2mat(chanceCI_L1O);
chanceMean_ISC = cell2mat(chanceMean_ISC);
chanceCI_ISC = cell2mat(chanceCI_ISC);
subplot(2,1,1)
plot(chanceMean_L1O, 'k');
patch('XData', [1:length(mlb.likeTimeVect), length(mlb.likeTimeVect):-1:1],...
    'YData', [(chanceMean_L1O+chanceCI_L1O)', flipud(chanceMean_L1O-chanceCI_L1O)'],...
    'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25, 'edgealpha', 0.5);
axis tight
title('Fully InSeq Correct: Leave-One-Out Decoding)');
ylabel('Probability');
set(gca, 'xticklabel', []);

subplot(2,1,2)
plot(chanceMean_ISC, 'k');
patch('XData', [1:length(mlb.likeTimeVect), length(mlb.likeTimeVect):-1:1],...
    'YData', [(chanceMean_ISC+chanceCI_ISC)', flipud(chanceMean_ISC-chanceCI_ISC)'],...
    'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25, 'edgealpha', 0.5);
axis tight
title('Fully InSeq Correct: Decode Remaining InSeq Correct)');
ylabel('Probability');
set(gca, 'xticklabel', []);

linkaxes;

%% Plot Probability summary figure
cMap = load('roma.mat'); % flip
% cMap = load('nuuk.mat');
% cMap = load('imola.mat');
% cMap = load('lapaz.mat'); %flip
cMap = cMap.(cell2mat(fieldnames(cMap)));
cMap = flipud(cMap);

piNdx = find(mlb.likeTimeVect==0)+0.5;
poNdx = find(mlb.likeTimeVect==nearestPOtime)+0.5;
rwdNdx = find(mlb.likeTimeVect==nearestRWDtime)+0.5;
posNdx = find(diff(mlb.likeTimeVect)<0)+0.5;
for dataToPlot = 1:2
    if dataToPlot == 1
        tipPostReal = groupPostTIP_RealL1O;
        tipPostChance = groupPostTIP_ChanceL1O;
        
        timePostReal = groupPostTime_RealL1O;
        timePostChance = groupPostTime_ChanceL1O;
        
        posPostReal = groupPostPos_RealL1O;
        posPostChance = groupPostPos_ChanceL1O;
        
    elseif dataToPlot == 2
        tipPostReal = groupPostTIP_RealISC;
        tipPostChance = groupPostTIP_ChanceISC;
        
        timePostReal = groupPostTime_RealISC;
        timePostChance = groupPostTime_ChanceISC;
        
        posPostReal = groupPostPos_RealISC;
        posPostChance = groupPostPos_ChanceISC;
    end
    
    figure;
    colormap(cMap);
    z = colormap;
    for pos = 1:mlb.seqLength
        subplot(5,5,sub2ind([mlb.seqLength+1,5],1,pos))
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect,mean(timePostReal{pos},3), [0 0.125]);
        hold on;
        plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(nearestPOtime, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(nearestRWDtime, [1,2]), ':k','linewidth', 2);
        if pos==mlb.seqLength
            cb = colorbar;
            cb.Location = 'westoutside';
            cb.Label.String = 'P(T|Spks,Pos)';
            cb.Label.Position(1) = 0;
            cb.Label.FontWeight = 'Bold';
            cb.Ticks = [0 0.125];
        end
        set(gca, 'xticklabel', [], 'yticklabel', []);
        title(sprintf('Time Info %i', pos));
        tempTimeReal = mean(timePostReal{pos},3)-(tinv(0.975,size(timePostReal{pos},3)-1).*mlb.SEMcalc(timePostReal{pos},0,3));
        tempTimeThresh = timePostChance{pos,1}+(tinv(0.975,numChancePerms-1).*(timePostChance{pos,1}./sqrt(numChancePerms-1)));
        abvThresh = tempTimeReal>tempTimeThresh;
        bounds = bwboundaries(abvThresh);
        for b = 1:length(bounds)
            plot(mlb.obsvTimeVect(bounds{b}(:,2)), mlb.obsvTimeVect(bounds{b}(:,1)), 'k', 'linewidth', 2);
        end
        
        subplot(5,5,sub2ind([mlb.seqLength+1,5],2:mlb.seqLength+1,repmat(pos, [1,mlb.seqLength])))
        imagesc(mean(tipPostReal{pos},3),[0 0.05]);
        hold on;
        plot(get(gca, 'xlim'),repmat(piNdx(1), [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(poNdx(1), [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(rwdNdx(1), [1,2]), ':k','linewidth', 2);
        for ndx = 1:length(piNdx)
            plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
            
            if ndx<length(piNdx)
                plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
            end
        end
        if pos==mlb.seqLength
            cb = colorbar;
            cb.Location = 'eastoutside';
            cb.Label.String = 'P(T,Pos|Spks)';
            cb.Label.Position(1) = 0;
            cb.Label.FontWeight = 'Bold';
            cb.Ticks = [0 0.05];
        end
        set(gca, 'xticklabel', [], 'yticklabel', []);
        title(sprintf('Time in Pos Info Position %i', pos));
        tempTIPReal = mean(tipPostReal{pos},3)-(tinv(0.975,size(tipPostReal{pos},3)-1).*mlb.SEMcalc(tipPostReal{pos},0,3));
        tempTIPThresh = tipPostChance{pos,1}+(tinv(0.975,numChancePerms-1).*(tipPostChance{pos,1}./sqrt(numChancePerms-1)));
        abvThresh = tempTIPReal>tempTIPThresh;
        bounds = bwboundaries(abvThresh);
        for b = 1:length(bounds)
            plot(bounds{b}(:,2), bounds{b}(:,1), 'k', 'linewidth', 2);
        end
        
        
        subplot(5,5,sub2ind([mlb.seqLength+1,5],pos+1,mlb.seqLength+1))
        tempChanceMean = posPostChance{pos,1}(:,pos);
        tempChanceSEM = posPostChance{pos,2}(:,pos)./sqrt(numChancePerms-1);
        tempChanceCI = tinv(0.975, numChancePerms-1).*tempChanceSEM;
        tempReal = posPostReal{pos};
        tempRealMean = mean(tempReal,3, 'omitnan');
        tempRealSEM = mlb.SEMcalc(tempReal,0,3);
        tempRealCI = tinv(0.975, size(tempReal,3)-1).*tempRealSEM;
        for p = 1:mlb.seqLength
            plot(mlb.obsvTimeVect, tempRealMean(:,p), 'color', mlb.PositionColors(p,:), 'linewidth', 1.5);
            hold on;
            patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
                'YData', [(tempRealMean(:,p)+tempRealSEM(:,p))', flipud(tempRealMean(:,p)-tempRealSEM(:,p))'],...
                'linestyle', 'none', 'edgecolor', mlb.PositionColors(p,:), 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0.25);
            patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
                'YData', [(tempRealMean(:,p)+tempRealCI(:,p))', flipud(tempRealMean(:,p)-tempRealCI(:,p))'],...
                'linestyle', ':', 'linewidth', 1.5, 'edgecolor', mlb.PositionColors(p,:), 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0);
            
            plot(mlb.obsvTimeVect, tempChanceMean, 'color', 'k', 'linewidth', 1.5);
            patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
                'YData', [(tempChanceMean+tempChanceSEM)', flipud(tempChanceMean-tempChanceSEM)'],...
                'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
            patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
                'YData', [(tempChanceMean+tempChanceCI)', flipud(tempChanceMean-tempChanceCI)'],...
                'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0);
            set(gca, 'ylim', [0 1]);
            
            plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
            if pos==1
                sortedYticks = sortrows([{0};{tempChanceMean(1)}; {0.5};{1}]);
                sortedYtickLabels = sortedYticks;
                sortedYtickLabels{cellfun(@(a)a==tempChanceMean(1),sortedYticks)} = 'Chance';
                set(gca, 'ytick', cell2mat(sortedYticks), 'yticklabel', sortedYtickLabels);
                ylabel('P(Pos|Spks,Pos)');
            else
                set(gca, 'yticklabel', []);
            end
            title(sprintf('Pos Info Position %i', pos));
            
        end
    end
    if dataToPlot == 1
        annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
            'String', sprintf("Train:FISC Test:FISC via Leave-One-Out; binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), BayesType = %.0f",...
            binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, bayesType),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
    else
        annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
            'String', sprintf("Train:FISC Test:ISC; binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), BayesType = %.0f",...
            binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, bayesType),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
    end
end

%% Plot Decoding summary figure
for dataToPlot = 1:2
    if dataToPlot == 1
        tipPostReal = groupDecodeTIP_RealL1O;
        tipPostChance = groupDecodeTIP_ChanceL1O;
        
        timePostReal = groupDecodeTime_RealL1O;
        timePostChance = groupDecodeTime_ChanceL1O;
        
        posPostReal = groupDecodePos_RealL1O;
        posPostChance = groupDecodePos_ChanceL1O;
        
    elseif dataToPlot == 2
        tipPostReal = groupDecodeTIP_RealISC;
        tipPostChance = groupDecodeTIP_ChanceISC;
        
        timePostReal = groupDecodeTime_RealISC;
        timePostChance = groupDecodeTime_ChanceISC;
        
        posPostReal = groupDecodePos_RealISC;
        posPostChance = groupDecodePos_ChanceISC;
    end
    
    figure;
    colormap(cMap);
    for pos = 1:mlb.seqLength
        subplot(5,5,sub2ind([mlb.seqLength+1,5],1,pos))
        imagesc(mlb.obsvTimeVect, mlb.obsvTimeVect,timePostReal{pos}, [0 0.125]);
        hold on;
        plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
        plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
        plot(get(gca, 'xlim'),zeros(1,2), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(nearestPOtime, [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(nearestRWDtime, [1,2]), ':k','linewidth', 2);
        if pos==mlb.seqLength
            cb = colorbar;
            cb.Location = 'westoutside';
            cb.Label.String = 'Accuracy';
            cb.Label.Position(1) = 0;
            cb.Label.FontWeight = 'Bold';
            cb.Ticks = [0 0.125];
        end
        set(gca, 'xticklabel', [], 'yticklabel', []);
        title(sprintf('Time Decoding %i', pos));
        tempTimeReal = timePostReal{pos};
        tempTimeThresh = timePostChance{pos,1}+(tinv(0.975,numChancePerms-1).*(timePostChance{pos,1}./sqrt(numChancePerms-1)));
        abvThresh = tempTimeReal>tempTimeThresh;
        bounds = bwboundaries(abvThresh);
        for b = 1:length(bounds)
            plot(mlb.obsvTimeVect(bounds{b}(:,2)), mlb.obsvTimeVect(bounds{b}(:,1)), 'k', 'linewidth', 2);
        end
        
        subplot(5,5,sub2ind([mlb.seqLength+1,5],2:mlb.seqLength+1,repmat(pos, [1,mlb.seqLength])))
        imagesc(tipPostReal{pos},[0 0.05]);
        hold on;
        plot(get(gca, 'xlim'),repmat(piNdx(1), [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(poNdx(1), [1,2]), '--k','linewidth', 2);
        plot(get(gca, 'xlim'),repmat(rwdNdx(1), [1,2]), ':k','linewidth', 2);
        for ndx = 1:length(piNdx)
            plot(repmat(piNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(poNdx(ndx), [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(rwdNdx(ndx), [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
            
            if ndx<length(piNdx)
                plot(repmat(posNdx(ndx), [1,2]), get(gca, 'ylim'), '-k','linewidth', 2);
            end
        end
        if pos==mlb.seqLength
            cb = colorbar;
            cb.Location = 'eastoutside';
            cb.Label.String = 'Accuracy';
            cb.Label.Position(1) = 0;
            cb.Label.FontWeight = 'Bold';
            cb.Ticks = [0 0.05];
        end
        set(gca, 'xticklabel', [], 'yticklabel', []);
        title(sprintf('Time in Pos Decoding Position %i', pos));
        tempTIPReal = tipPostReal{pos};
        tempTIPThresh = tipPostChance{pos,1}+(tinv(0.975,numChancePerms-1).*(tipPostChance{pos,1}./sqrt(numChancePerms-1)));
        abvThresh = tempTIPReal>tempTIPThresh;
        bounds = bwboundaries(abvThresh);
        for b = 1:length(bounds)
            plot(bounds{b}(:,2), bounds{b}(:,1), 'k', 'linewidth', 2);
        end
        
        
        subplot(5,5,sub2ind([mlb.seqLength+1,5],pos+1,mlb.seqLength+1))
        tempChanceMean = posPostChance{pos,1}(:,pos);
        tempChanceSEM = posPostChance{pos,2}(:,pos)./sqrt(numChancePerms-1);
        tempChanceCI = tinv(0.975, numChancePerms-1).*tempChanceSEM;
        for p = 1:mlb.seqLength
            plot(mlb.obsvTimeVect, posPostReal{pos}(:,p), 'color', mlb.PositionColors(p,:), 'linewidth', 1.5);
            hold on;
           
            plot(mlb.obsvTimeVect, tempChanceMean, 'color', 'k', 'linewidth', 1.5);
            patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
                'YData', [(tempChanceMean+tempChanceSEM)', flipud(tempChanceMean-tempChanceSEM)'],...
                'linestyle', 'none', 'edgecolor', 'k', 'facecolor', 'k', 'facealpha', 0.25);
            patch('XData', [mlb.obsvTimeVect; flipud(mlb.obsvTimeVect)],...
                'YData', [(tempChanceMean+tempChanceCI)', flipud(tempChanceMean-tempChanceCI)'],...
                'linestyle', ':', 'linewidth', 1.5, 'edgecolor', 'k', 'facecolor', mlb.PositionColors(p,:), 'facealpha', 0);
            set(gca, 'ylim', [0 1]);
            
            plot(zeros(1,2), get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(nearestPOtime, [1,2]), get(gca, 'ylim'), '--k','linewidth', 2);
            plot(repmat(nearestRWDtime, [1,2]), get(gca, 'ylim'), ':k','linewidth', 2);
            if pos==1
                sortedYticks = sortrows([{0};{tempChanceMean(1)}; {0.5};{1}]);
                sortedYtickLabels = sortedYticks;
                sortedYtickLabels{cellfun(@(a)a==tempChanceMean(1),sortedYticks)} = 'Chance';
                set(gca, 'ytick', cell2mat(sortedYticks), 'yticklabel', sortedYtickLabels);
                ylabel('Accuracy');
            else
                set(gca, 'yticklabel', []);
            end
            title(sprintf('Pos Decoding Position %i', pos));
            
        end
    end
    if dataToPlot == 1
        annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
            'String', sprintf("Train:FISC Decode:FISC via Leave-One-Out; binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), BayesType = %.0f",...
            binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, bayesType),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
    else
        annotation(gcf,'textbox', [0.1 0.95 0.9 0.05],...
            'String', sprintf("Train:FISC Decode:ISC; binSize = %.0fms, dsRate = %.0fms, Trial Window = (%.0fms:%.0fms to %s), BayesType = %.0f",...
            binSize, dsRate, trlWindow{1}(1), trlWindow{1}(2), alignment{1}, bayesType),...
            'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left');
    end
end
